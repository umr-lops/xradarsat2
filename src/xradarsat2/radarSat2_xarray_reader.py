# imports
import ast
import glob
import os
import re
import traceback
from datetime import datetime

import dask
import numpy as np
import rasterio
import rioxarray
import xarray as xr
import xmltodict
import yaml
from affine import Affine

xpath_dict = {
    "geolocation_grid": {
        "xpath": (
            "/product/imageAttributes/geographicInformation/"
            "geolocationGrid/imageTiePoint"
        )
    },
    "orbit_information": {
        "xpath": "/product/sourceAttributes/orbitAndAttitude/orbitInformation",
    },
    "attitude_information": {
        "xpath": "/product/sourceAttributes/orbitAndAttitude/attitudeInformation"
    },
    "doppler": {"xpath": "/product/imageGenerationParameters"},
    "radarParameters": {"xpath": "/product/sourceAttributes/radarParameters"},
}

radar_parameters_key_dict = {
    "ds_attributes": [
        "acquisitionType",
        "beams",
        "polarizations",
        "pulses",
        "radarCenterFrequency",
        "antennaPointing",
        "yawSteeringFlag",
        "geodeticFlag",
        "rawBitsPerSample",
        "pulseLength",
        "pulseBandwidth",
        "adcSamplingRate",
    ],
    "coord": {
        "pulsesReceivedPerDwell": ["beam"],
        "numberOfPulseIntervalsPerDwell": ["beam"],
        "rank": ["beam"],
        "settableGain": ["beam", "pole"],
        "pulseRepetitionFrequency": ["beam"],
        "samplesPerEchoLine": ["beam"],
        "incidenceAngleCorrection_Beta_Nought": [],
        "incidenceAngleCorrection_Sigma_Nought": [],
        "incidenceAngleCorrection_Gamma": [],
    },
}


def xpath_get(mydict, xpath):
    """
    Return a sub dictionary thanks to a xPath

    Parameters
    ----------
    mydict: dict
        Content of product.xml as a dictionary
    xpath: str
        xPath that shows the location of the dataset in the
        product.xml hierarchy

    Returns
    -------
    dict
        product.xml dataset information translated as a dictionary
    """

    elem = mydict
    try:
        for x in xpath.strip("/").split("/"):
            elem = elem.get(x)
    except (ImportError, NameError, KeyError, TypeError, ValueError, AttributeError):
        pass
    return elem


def create_dic_geolocation_grid(dictio):
    """
    Create a dictionary containing useful information
     of Geolocation Grid dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Useful information of Geolocation Grid as a dictionary
    """

    content_dict = xpath_get(dictio, xpath_dict["geolocation_grid"]["xpath"])
    final_dc = {
        "latitude": {"values": [], "units": ""},
        "longitude": {"values": [], "units": ""},
        "height": {"values": [], "units": ""},
        "coords": {"lines": [], "pixels": []},
        "attr": get_line_and_pix_info(fill_image_attribute(dictio)),
    }
    for element in content_dict:
        final_dc["coords"]["lines"].append(
            parse_value(element["imageCoordinate"]["line"])
        )
        final_dc["coords"]["pixels"].append(
            parse_value(element["imageCoordinate"]["pixel"])
        )
        final_dc["longitude"]["values"].append(
            parse_value(element["geodeticCoordinate"]["longitude"]["#text"])
        )
        final_dc["latitude"]["values"].append(
            parse_value(element["geodeticCoordinate"]["latitude"]["#text"])
        )
        final_dc["height"]["values"].append(
            parse_value(element["geodeticCoordinate"]["height"]["#text"])
        )
        final_dc["longitude"]["units"] = element["geodeticCoordinate"]["longitude"][
            "@units"
        ]
        final_dc["latitude"]["units"] = element["geodeticCoordinate"]["latitude"][
            "@units"
        ]
        final_dc["height"]["units"] = element["geodeticCoordinate"]["height"]["@units"]
    return final_dc


def create_dataset_geolocation_grid(dictio, folder_path):
    """
    Create a dataset for Geolocation Grid

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Geolocation Grid dataset
    """

    ds = xr.Dataset()
    lines = dictio["coords"]["lines"]
    pixs = dictio["coords"]["pixels"]
    lines = [int(float(lines[k])) for k in range(len(lines))]
    pixs = [int(float(pixs[k])) for k in range(len(pixs))]
    for key in dictio:
        if (key != "coords") and (key != "attr"):
            data = create_2d_matrix(lines, pixs, dictio[key]["values"])
            unit = dictio[key]["units"]
            if key == "longitude":
                xpath_suffix = os.path.join("geodeticCoordinate", "longitude")
            elif key == "latitude":
                xpath_suffix = os.path.join("geodeticCoordinate", "latitude")
            elif key == "height":
                xpath_suffix = os.path.join("geodeticCoordinate", "height")
            xpath = os.path.join(xpath_dict["geolocation_grid"]["xpath"], xpath_suffix)
            da = xr.DataArray(
                data=data,
                name=key,
                coords={
                    "line": np.unique(np.array(lines)),
                    "pixel": np.unique(np.array(pixs)),
                },
                dims=["line", "pixel"],
                attrs=(
                    {"units": unit, "xpath": xpath}
                    | generate_doc_vars(xpath, folder_path)
                ),
            )
            ds[key] = da
    try:
        ds.attrs = generate_doc_ds(xpath_dict["geolocation_grid"]["xpath"], folder_path)
    except Exception:
        pass
    ds["line"].attrs = dictio["attr"]["line"]
    ds["pixel"].attrs = dictio["attr"]["pixel"]
    return ds


def get_line_and_pix_info(dictio):
    """
    get line and pixel spacing information

    Parameters
    ----------
    dictio: dict
        Content of image/raster attributes as a dictionary

    Returns
    -------
    dict
        Dictionary of line and pixel spacing information
    """

    line = {}
    pixel = {}
    for key in dictio:
        if "PixelSpacing" in key:
            pixel[key] = dictio[key]
        if "LineSpacing" in key:
            line[key] = dictio[key]
    return {"line": line, "pixel": pixel}


def get_dic_orbit_information(dictio):
    """
    Create a dictionary containing useful information of Orbit Information dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Useful information of Orbit Information as a dictionary
    """

    content_dict = xpath_get(dictio, xpath_dict["orbit_information"]["xpath"])
    ds_attr = {}
    timestamp = []
    xPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "xPosition"
            ),
        },
    }
    yPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "yPosition"
            ),
        },
    }
    zPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "zPosition"
            ),
        },
    }
    xVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "xVelocity"
            ),
        },
    }
    yVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "yVelocity"
            ),
        },
    }
    zVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["geolocation_grid"]["xpath"], "stateVector", "zVelocity"
            ),
        },
    }
    for key in content_dict:
        if isinstance(content_dict[key], str):
            ds_attr[key] = content_dict[key]
        elif isinstance(content_dict[key], list):
            for value in content_dict[key]:
                dt_stamp = datetime.strptime(
                    value["timeStamp"], "%Y-%m-%dT%H:%M:%S.%fZ"
                )
                timestamp.append(np.array(dt_stamp).astype("datetime64[ns]"))
                xPosition["values"].append(float(value["xPosition"]["#text"]))
                xPosition["attr"]["units"] = value["xPosition"]["@units"]
                yPosition["values"].append(float(value["yPosition"]["#text"]))
                yPosition["attr"]["units"] = value["yPosition"]["@units"]
                zPosition["values"].append(float(value["zPosition"]["#text"]))
                zPosition["attr"]["units"] = value["zPosition"]["@units"]
                xVelocity["values"].append(float(value["xVelocity"]["#text"]))
                xVelocity["attr"]["units"] = value["xVelocity"]["@units"]
                yVelocity["values"].append(float(value["yVelocity"]["#text"]))
                yVelocity["attr"]["units"] = value["yVelocity"]["@units"]
                zVelocity["values"].append(float(value["zVelocity"]["#text"]))
                zVelocity["attr"]["units"] = value["zVelocity"]["@units"]
    return {
        "ds_attr": ds_attr,
        "timestamp": timestamp,
        "xPosition": xPosition,
        "yPosition": yPosition,
        "zPosition": zPosition,
        "xVelocity": xVelocity,
        "yVelocity": yVelocity,
        "zVelocity": zVelocity,
    }


def create_dataset_orbit_information(
    ds_attr, timestamp, xPos, yPos, zPos, xVel, yVel, zVel, folder_path
):
    """
    Create a dataset for Orbit Information

    Parameters
    ----------
    ds_attr: dict
        Dictionary of dataset attributes
    timestamp: list[datetime64]
        Timestamp list of values
    xPos: dict[str, list | dict[str, str]]
        xPosition list of values
    yPos: dict[str, list | dict[str, str]]
        yPosition list of values
    zPos: dict[str, list | dict[str, str]]
        zPosition list of values
    xVel: dict[str, list | dict[str, str]]
        xVelocity list of values
    yVel: dict[str, list | dict[str, str]]
        yVelocity list of values
    zVel: dict[str, list | dict[str, str]]
        zVelocity list of values
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Orbit Information dataset
    """

    ds = xr.Dataset()
    xpos_da = xr.DataArray(
        data=xPos["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(xPos["attr"] | generate_doc_vars(xPos["attr"]["xpath"], folder_path)),
    )
    ypos_da = xr.DataArray(
        data=yPos["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(yPos["attr"] | generate_doc_vars(yPos["attr"]["xpath"], folder_path)),
    )
    zpos_da = xr.DataArray(
        data=zPos["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(zPos["attr"] | generate_doc_vars(zPos["attr"]["xpath"], folder_path)),
    )
    xvel_da = xr.DataArray(
        data=xVel["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(xVel["attr"] | generate_doc_vars(xVel["attr"]["xpath"], folder_path)),
    )
    yvel_da = xr.DataArray(
        data=yVel["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(yVel["attr"] | generate_doc_vars(yVel["attr"]["xpath"], folder_path)),
    )
    zvel_da = xr.DataArray(
        data=zVel["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(zVel["attr"] | generate_doc_vars(zVel["attr"]["xpath"], folder_path)),
    )
    ds["xPosition"] = xpos_da
    ds["yPosition"] = ypos_da
    ds["zPosition"] = zpos_da
    ds["xVelocity"] = xvel_da
    ds["yVelocity"] = yvel_da
    ds["zVelocity"] = zvel_da
    try:
        ds.attrs = ds_attr | generate_doc_ds(
            xpath_dict["orbit_information"]["xpath"], folder_path
        )
        ds.attrs["Description"] += (
            ". " + generate_doc_ds("stateVector", folder_path)["Description"]
        )
    except Exception:
        pass
    return ds


def get_dic_attitude_info(dictio):
    """
    Create a dictionary containing useful information of Attitude Information dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Useful information of Attitude Information as a dictionary
    """

    content_dict = xpath_get(dictio, xpath_dict["attitude_information"]["xpath"])
    ds_attr = {}
    timestamp = []
    yaw = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["attitude_information"]["xpath"], "attitudeAngles", "yaw"
            ),
        },
    }
    roll = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["attitude_information"]["xpath"], "attitudeAngles", "roll"
            ),
        },
    }
    pitch = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": os.path.join(
                xpath_dict["attitude_information"]["xpath"], "attitudeAngles", "pitch"
            ),
        },
    }
    for key in content_dict:
        if isinstance(content_dict[key], str):
            ds_attr[key] = content_dict[key]
        elif isinstance(content_dict[key], list):
            for value in content_dict[key]:
                dt_stamp = datetime.strptime(
                    value["timeStamp"], "%Y-%m-%dT%H:%M:%S.%fZ"
                )
                timestamp.append(np.array(dt_stamp).astype("datetime64[ns]"))
                yaw["values"].append(float(value["yaw"]["#text"]))
                yaw["attr"]["units"] = value["yaw"]["@units"]
                roll["values"].append(float(value["roll"]["#text"]))
                roll["attr"]["units"] = value["roll"]["@units"]
                pitch["values"].append(float(value["pitch"]["#text"]))
                pitch["attr"]["units"] = value["pitch"]["@units"]
    return {
        "ds_attr": ds_attr,
        "timestamp": timestamp,
        "yaw": yaw,
        "roll": roll,
        "pitch": pitch,
    }


def create_dataset_attitude_information(
    ds_attr, timestamp, yaw, roll, pitch, folder_path
):
    """
    Create a dataset for Attitude Information

    Parameters
    ----------
    ds_attr: dict
        Dictionary of dataset attributes
    timestamp: list[datetime64]
        Timestamp list of values
    yaw: dict[str, list | dict[str, str]]
        Yaw list of values
    roll: dict[str, list | dict[str, str]]
        Roll list of values
    pitch: dict[str, list | dict[str, str]]
        Pitch list of values
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Attitude Information dataset
    """

    ds = xr.Dataset()
    interesting_files_yaw = list_xsd_files(yaw["attr"]["xpath"], folder_path)
    yaw_doc = find_doc_in_xsd_files(yaw["attr"]["xpath"], interesting_files_yaw)
    interesting_files_roll = list_xsd_files(roll["attr"]["xpath"], folder_path)
    roll_doc = find_doc_in_xsd_files(roll["attr"]["xpath"], interesting_files_roll)
    interesting_files_pitch = list_xsd_files(pitch["attr"]["xpath"], folder_path)
    pitch_doc = find_doc_in_xsd_files(pitch["attr"]["xpath"], interesting_files_pitch)
    yaw_da = xr.DataArray(
        data=yaw["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(yaw["attr"] | yaw_doc),
    )
    roll_da = xr.DataArray(
        data=roll["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(roll["attr"] | roll_doc),
    )
    pitch_da = xr.DataArray(
        data=pitch["values"],
        coords={"timeStamp": timestamp},
        dims="timeStamp",
        attrs=(pitch["attr"] | pitch_doc),
    )
    ds["yaw"] = yaw_da
    ds["roll"] = roll_da
    ds["pitch"] = pitch_da
    try:
        ds.attrs = ds_attr | generate_doc_ds(
            xpath_dict["attitude_information"]["xpath"], folder_path
        )
    except Exception:
        pass
    return ds


def get_dict_doppler_centroid(dictio):
    """
    Create a dictionary containing useful information of Doppler Centroid dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary
    Returns
    -------
    dict
        Useful information of Doppler Centroid as a dictionary
    """

    content_dict = xpath_get(dictio, xpath_dict["doppler"]["xpath"])
    ds_attr = {}
    times = []
    Ambiguity = {"values": [], "attr": {}}
    AmbiguityConfidence = {"values": [], "attr": {}}
    CentroidReferenceTime = {"values": [], "attr": {}}
    CentroidPolynomialPeriod = {"values": [], "attr": {}}
    CentroidCoefficients = {"values": [], "attr": {}}
    CentroidConfidence = {"values": [], "attr": {}}
    for key in content_dict:
        xpath = xpath_dict["doppler"]["xpath"]
        if key == "dopplerCentroid":
            xpath = os.path.join(xpath, key)
            for value in content_dict[key]:
                dt_dce = datetime.strptime(
                    value["timeOfDopplerCentroidEstimate"], "%Y-%m-%dT%H:%M:%S.%fZ"
                )
                times.append(np.array(dt_dce).astype("datetime64[ns]"))
                Ambiguity["values"].append(int(value["dopplerAmbiguity"]))
                Ambiguity["attr"]["xpath"] = os.path.join(xpath, "dopplerAmbiguity")
                AmbiguityConfidence["values"].append(
                    float(value["dopplerAmbiguityConfidence"])
                )
                AmbiguityConfidence["attr"]["xpath"] = os.path.join(
                    xpath, "dopplerAmbiguityConfidence"
                )
                CentroidReferenceTime["values"].append(
                    float(value["dopplerCentroidReferenceTime"]["#text"])
                )
                CentroidReferenceTime["attr"]["units"] = value[
                    "dopplerCentroidReferenceTime"
                ]["@units"]
                CentroidReferenceTime["attr"]["xpath"] = os.path.join(
                    xpath, "dopplerCentroidReferenceTime"
                )
                CentroidPolynomialPeriod["values"].append(
                    float(value["dopplerCentroidPolynomialPeriod"]["#text"])
                )
                CentroidPolynomialPeriod["attr"]["units"] = value[
                    "dopplerCentroidPolynomialPeriod"
                ]["@units"]
                CentroidPolynomialPeriod["attr"]["xpath"] = os.path.join(
                    xpath, "dopplerCentroidPolynomialPeriod"
                )
                CentroidCoefficients["values"].append(
                    [float(x) for x in value["dopplerCentroidCoefficients"].split(" ")]
                )
                CentroidCoefficients["attr"]["xpath"] = os.path.join(
                    xpath, "dopplerCentroidCoefficients"
                )
                CentroidConfidence["values"].append(
                    float(value["dopplerCentroidConfidence"])
                )
                CentroidConfidence["attr"]["xpath"] = os.path.join(
                    xpath, "dopplerCentroidConfidence"
                )
        elif len(times) != 0:
            # Doppler Centroid is the first key, so we want to skip the further keys (chirp...)
            break
    return {
        "ds_attr": ds_attr,
        "timeOfDopplerCentroidEstimate": times,
        "dopplerAmbiguity": Ambiguity,
        "dopplerAmbiguityConfidence": AmbiguityConfidence,
        "dopplerCentroidReferenceTime": CentroidReferenceTime,
        "dopplerCentroidPolynomialPeriod": CentroidPolynomialPeriod,
        "dopplerCentroidCoefficients": CentroidCoefficients,
        "dopplerCentroidConfidence": CentroidConfidence,
    }


def create_dataset_doppler_centroid(
    ds_attr,
    times,
    Ambiguity,
    AmbiguityConfidence,
    CentroidReferenceTime,
    CentroidPolynomialPeriod,
    CentroidCoefficients,
    CentroidConfidence,
    folder_path,
):
    """
    Create a dataset for Doppler Centroid

    Parameters
    ----------
    ds_attr: dict[str, str]
        Dictionary of dataset attributes
    times: list[datetime64]
        Times list of values
    Ambiguity: dict[str, dict | list | str]
        Ambiguity list of values
    AmbiguityConfidence: dict[str, dict | list | str]
        AmbiguityConfidence list of values
    CentroidReferenceTime: dict[str, dict | list | str]
        CentroidReferenceTime list of values
    CentroidPolynomialPeriod: dict[str, dict | list | str]
        CentroidPolynomialPeriod list of values
    CentroidCoefficients: dict[str, dict | list | str]
        CentroidCoefficients list of values
    CentroidConfidence: dict[str, dict | list | str]
        CentroidConfidence list of values
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Doppler Centroid dataset
    """

    ds = xr.Dataset()
    Ambiguity_doc = find_doc_in_xsd_files(
        Ambiguity["attr"]["xpath"],
        list_xsd_files(Ambiguity["attr"]["xpath"], folder_path),
    )
    AmbiguityConfidence_doc = find_doc_in_xsd_files(
        AmbiguityConfidence["attr"]["xpath"],
        list_xsd_files(AmbiguityConfidence["attr"]["xpath"], folder_path),
    )
    CentroidReferenceTime_doc = find_doc_in_xsd_files(
        CentroidReferenceTime["attr"]["xpath"],
        list_xsd_files(CentroidReferenceTime["attr"]["xpath"], folder_path),
    )
    CentroidPolynomialPeriod_doc = find_doc_in_xsd_files(
        CentroidPolynomialPeriod["attr"]["xpath"],
        list_xsd_files(CentroidPolynomialPeriod["attr"]["xpath"], folder_path),
    )
    CentroidCoefficients_doc = find_doc_in_xsd_files(
        CentroidCoefficients["attr"]["xpath"],
        list_xsd_files(CentroidCoefficients["attr"]["xpath"], folder_path),
    )
    CentroidConfidence_doc = find_doc_in_xsd_files(
        CentroidConfidence["attr"]["xpath"],
        list_xsd_files(CentroidConfidence["attr"]["xpath"], folder_path),
    )
    ambiguity_da = xr.DataArray(
        data=Ambiguity["values"],
        coords={"timeOfDopplerCentroidEstimate": times},
        dims=["timeOfDopplerCentroidEstimate"],
        attrs=(Ambiguity["attr"] | Ambiguity_doc),
    )
    ambiguityConfidence_da = xr.DataArray(
        data=AmbiguityConfidence["values"],
        coords={"timeOfDopplerCentroidEstimate": times},
        dims=["timeOfDopplerCentroidEstimate"],
        attrs=(AmbiguityConfidence["attr"] | AmbiguityConfidence_doc),
    )
    centroidReferenceTime_da = xr.DataArray(
        data=CentroidReferenceTime["values"],
        coords={"timeOfDopplerCentroidEstimate": times},
        dims=["timeOfDopplerCentroidEstimate"],
        attrs=(CentroidReferenceTime["attr"] | CentroidReferenceTime_doc),
    )
    centroidPolynomialPeriod_da = xr.DataArray(
        data=CentroidPolynomialPeriod["values"],
        coords={"timeOfDopplerCentroidEstimate": times},
        dims=["timeOfDopplerCentroidEstimate"],
        attrs=(CentroidPolynomialPeriod["attr"] | CentroidPolynomialPeriod_doc),
    )
    centroidCoefficients_da = xr.DataArray(
        data=np.array(CentroidCoefficients["values"]),
        coords={
            "timeOfDopplerCentroidEstimate": times,
            "n-Coefficients": [
                i for i in range(np.array(CentroidCoefficients["values"]).shape[1])
            ],
        },
        dims=["timeOfDopplerCentroidEstimate", "n-Coefficients"],
        attrs=(CentroidCoefficients["attr"] | CentroidCoefficients_doc),
    )
    centroidConfidence_da = xr.DataArray(
        data=CentroidConfidence["values"],
        coords={"timeOfDopplerCentroidEstimate": times},
        dims=["timeOfDopplerCentroidEstimate"],
        attrs=(CentroidConfidence["attr"] | CentroidConfidence_doc),
    )
    ds["dopplerAmbiguity"] = ambiguity_da
    ds["dopplerAmbiguityConfidence"] = ambiguityConfidence_da
    ds["dopplerCentroidReferenceTime"] = centroidReferenceTime_da
    ds["dopplerCentroidPolynomialPeriod"] = centroidPolynomialPeriod_da
    ds["dopplerCentroidCoefficients"] = centroidCoefficients_da
    ds["dopplerCentroidConfidence"] = centroidConfidence_da
    try:
        ds.attrs = ds_attr | generate_doc_ds(
            os.path.join(xpath_dict["doppler"]["xpath"], "dopplerCentroid"), folder_path
        )
    except Exception:
        pass
    return ds


def get_dic_doppler_rate_values(dictio):
    """
    Create a dictionary containing useful information of doppler rate values dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Useful information of doppler rate values as a dictionary
    """

    content_dict = xpath_get(dictio, xpath_dict["doppler"]["xpath"])
    ds_attr = {}
    RateReferenceTime = {"values": [], "attr": {}}
    RateValuesCoefficients = {"values": [], "attr": {}}
    for key in content_dict:
        xpath = xpath_dict["doppler"]["xpath"]
        if key == "dopplerRateValues":
            xpath = os.path.join(xpath, key)
            if isinstance(content_dict[key], dict):
                RateReferenceTime["values"].append(
                    float(content_dict[key]["dopplerRateReferenceTime"]["#text"])
                )
                RateReferenceTime["attr"]["RateReferenceTime units"] = content_dict[
                    key
                ]["dopplerRateReferenceTime"]["@units"]
                RateReferenceTime["attr"]["dopplerRateReferenceTime_xpath"] = (
                    os.path.join(xpath, "dopplerRateReferenceTime")
                )
                RateValuesCoefficients["values"].append(
                    [
                        float(x)
                        for x in content_dict[key][
                            "dopplerRateValuesCoefficients"
                        ].split(" ")
                    ]
                )
                RateValuesCoefficients["attr"][
                    "dopplerRateValuesCoefficients_xpath"
                ] = os.path.join(xpath, "dopplerRateValuesCoefficients")
            elif isinstance(content_dict[key], list):
                for value in content_dict[key]:
                    RateReferenceTime["values"].append(
                        float(content_dict[key]["dopplerRateReferenceTime"]["#text"])
                    )
                    RateReferenceTime["attr"]["units"] = content_dict[key][
                        "dopplerRateReferenceTime"
                    ]["@units"]
                    RateReferenceTime["attr"]["dopplerRateReferenceTime_xpath"] = (
                        os.path.join(xpath, "dopplerRateReferenceTime")
                    )
                    RateValuesCoefficients["values"].append(
                        [
                            float(x)
                            for x in content_dict[key][
                                "dopplerRateValuesCoefficients"
                            ].split(" ")
                        ]
                    )
                    RateValuesCoefficients["attr"][
                        "dopplerRateValuesCoefficients_xpath"
                    ] = os.path.join(xpath, "dopplerRateValuesCoefficients")
        elif len(RateReferenceTime["values"]) != 0:
            # Doppler Rate Values is in the first keys, so we want to skip the further keys (chirp...)
            break
    return {
        "ds_attr": ds_attr,
        "dopplerRateReferenceTime": RateReferenceTime,
        "dopplerRateValuesCoefficients": RateValuesCoefficients,
    }


def create_dataset_doppler_rate_values(
    ds_attr, rateTime, rateCoefficients, folder_path
):
    """
    Create a dataset for Doppler Rate Values

    Parameters
    ----------
    ds_attr: dict[str, str]
        Dictionary of dataset attributes
    rateTime: dict[str, dict | list | str]
        RateTime list of values
    rateCoefficients: dict[str, dict | list | str]
        RateCoefficients list of values
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Doppler Rate Values dataset
    """

    ds = xr.Dataset()
    rateTime_doc = find_doc_in_xsd_files(
        rateTime["attr"]["dopplerRateReferenceTime_xpath"],
        list_xsd_files(rateTime["attr"]["dopplerRateReferenceTime_xpath"], folder_path),
    )
    rateCoefficients_doc = find_doc_in_xsd_files(
        rateCoefficients["attr"]["dopplerRateValuesCoefficients_xpath"],
        list_xsd_files(
            rateCoefficients["attr"]["dopplerRateValuesCoefficients_xpath"], folder_path
        ),
    )
    rateCoefficients_da = xr.DataArray(
        data=np.array(rateCoefficients["values"]),
        coords={
            "dopplerRateReferenceTime": rateTime["values"],
            "n-RateValuesCoefficients": [
                i for i in range(np.array(rateCoefficients["values"]).shape[1])
            ],
        },
        dims=["dopplerRateReferenceTime", "n-RateValuesCoefficients"],
        attrs=(
            rateTime["attr"]
            | rateCoefficients["attr"]
            | rateTime_doc
            | rateCoefficients_doc
        ),
    )
    ds["dopplerRateValues"] = rateCoefficients_da
    try:
        ds.attrs = ds_attr | generate_doc_ds(
            os.path.join(xpath_dict["doppler"]["xpath"], "dopplerRateValues"),
            folder_path,
        )
    except Exception:
        pass
    return ds


def get_dict_chirp(dictio):
    """
    Create a dictionary containing useful information of chirp dataset

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Useful information of chirp as a dictionary
    """

    xpath = xpath_dict["doppler"]["xpath"]
    content_dict = xpath_get(dictio, xpath)
    pole = {"values": []}
    ds_attr = {}
    replicaQualityValid = {"values": [], "attr": {}}
    crossCorrelationWidth = {"values": [], "attr": {}}
    sideLobeLevel = {"values": [], "attr": {}}
    integratedSideLobeRatio = {"values": [], "attr": {}}
    crossCorrelationPeakLoc = {"values": [], "attr": {}}
    chirpPower = {"values": [], "attr": {}}
    amplitudeCoefficients = {"values": [], "attr": {}}
    phaseCoefficients = {"values": [], "attr": {}}
    for key in content_dict:
        if key == "chirp":
            xpath = os.path.join(xpath, key)
            if isinstance(content_dict[key], list):
                for value in content_dict[key]:
                    (
                        tmp_pole,
                        tmp_ds_attr,
                        tmp_replicaQualityValid,
                        tmp_crossCorrelationWidth,
                        tmp_sideLobeLevel,
                        tmp_integratedSideLobeRatio,
                        tmp_crossCorrelationPeakLoc,
                        tmp_chirpPower,
                        tmp_amplitudeCoefficients,
                        tmp_phaseCoefficients,
                    ) = chirp_dict_loop_processing(
                        value,
                        xpath,
                        pole,
                        ds_attr,
                        replicaQualityValid,
                        crossCorrelationWidth,
                        sideLobeLevel,
                        integratedSideLobeRatio,
                        crossCorrelationPeakLoc,
                        chirpPower,
                        amplitudeCoefficients,
                        phaseCoefficients,
                    )
            else:
                (
                    tmp_pole,
                    tmp_ds_attr,
                    tmp_replicaQualityValid,
                    tmp_crossCorrelationWidth,
                    tmp_sideLobeLevel,
                    tmp_integratedSideLobeRatio,
                    tmp_crossCorrelationPeakLoc,
                    tmp_chirpPower,
                    tmp_amplitudeCoefficients,
                    tmp_phaseCoefficients,
                ) = chirp_dict_loop_processing(
                    content_dict[key],
                    xpath,
                    pole,
                    ds_attr,
                    replicaQualityValid,
                    crossCorrelationWidth,
                    sideLobeLevel,
                    integratedSideLobeRatio,
                    crossCorrelationPeakLoc,
                    chirpPower,
                    amplitudeCoefficients,
                    phaseCoefficients,
                )
    new_ds_attr = {}
    for key in ds_attr:
        for intern_key in ds_attr[key]:
            value = parse_value(ds_attr[key][intern_key])
            new_ds_attr[f"{key}_{intern_key}"] = value

    return {
        "pole": pole,
        "ds_attr": new_ds_attr,
        "replicaQualityValid": replicaQualityValid,
        "crossCorrelationWidth": crossCorrelationWidth,
        "sideLobeLevel": sideLobeLevel,
        "integratedSideLobeRatio": integratedSideLobeRatio,
        "crossCorrelationPeakLoc": crossCorrelationPeakLoc,
        "chirpPower": chirpPower,
        "amplitudeCoefficients": amplitudeCoefficients,
        "phaseCoefficients": phaseCoefficients,
    }


def chirp_dict_loop_processing(
    dictio,
    xpath,
    pole,
    ds_attr,
    replicaQualityValid,
    crossCorrelationWidth,
    sideLobeLevel,
    integratedSideLobeRatio,
    crossCorrelationPeakLoc,
    chirpPower,
    amplitudeCoefficients,
    phaseCoefficients,
):
    """
    Processing of chirp intern loop to fill useful information in dictionaries

    Parameters
    ----------
    dictio: dict
        content of an intern dictionary key as a dictionary
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parents
    pole: dict
        Pole list of values
    ds_attr: dict
        Dictionary of dataset attributes
    replicaQualityValid: dict
        ReplicaQualityValid list of values
    crossCorrelationWidth: dict
        CrossCorrelationWidth list of values
    sideLobeLevel: dict
        SideLobeLevel list of values
    integratedSideLobeRatio: dict
        IntegratedSideLobeRatio list of values
    crossCorrelationPeakLoc: dict
        CrossCorrelationPeakLoc list of values
    chirpPower: dict
        ChirpPower list of values
    amplitudeCoefficients: dict
        AmplitudeCoefficients list of values
    phaseCoefficients: dict
        PhaseCoefficients list of values

    Returns
    -------
    dict
        Result dictionaries after one loop
    """

    pole["values"].append(dictio["@pole"])
    for value in dictio:
        if isinstance(dictio[value], str) and ("pole" not in value) and ("@" in value):
            if dictio["@pole"] not in list(ds_attr.keys()):
                ds_attr[dictio["@pole"]] = {}
            ds_attr[dictio["@pole"]][value.replace("@", "")] = dictio[value]
        elif (value == "amplitudeCoefficients") or (value == "phaseCoefficients"):
            eval(value)["values"].append([float(x) for x in dictio[value].split(" ")])
            eval(value)["attr"]["xpath"] = os.path.join(xpath, value)
        elif value == "chirpQuality":
            prefix_path = os.path.join(xpath, value)
            for var in dictio[value]:
                eval(var)["attr"]["xpath"] = os.path.join(prefix_path, var)
                if (
                    (var == "crossCorrelationPeakLoc")
                    or (var == "crossCorrelationWidth")
                    or (var == "replicaQualityValid")
                ):
                    eval(var)["values"].append(parse_value(dictio[value][var]))
                elif (var == "sideLobeLevel") or (var == "integratedSideLobeRatio"):
                    for intern_key in dictio[value][var]:
                        if "@" in intern_key:
                            eval(var)["attr"][intern_key.replace("@", "")] = dictio[
                                value
                            ][var][intern_key]
                        elif intern_key == "#text":
                            eval(var)["values"].append(
                                parse_value(dictio[value][var][intern_key])
                            )
        elif value == "chirpPower":
            eval(value)["attr"]["xpath"] = os.path.join(xpath, value)
            for intern_key in dictio[value]:
                if "@" in intern_key:
                    eval(value)["attr"][intern_key.replace("@", "")] = dictio[value][
                        intern_key
                    ]
                elif intern_key == "#text":
                    eval(value)["values"].append(parse_value(dictio[value][intern_key]))
    return (
        pole,
        ds_attr,
        replicaQualityValid,
        crossCorrelationWidth,
        sideLobeLevel,
        integratedSideLobeRatio,
        crossCorrelationPeakLoc,
        chirpPower,
        amplitudeCoefficients,
        phaseCoefficients,
    )


def create_dataset_chirp(
    pole,
    ds_attr,
    replicaQualityValid,
    crossCorrelationWidth,
    sideLobeLevel,
    integratedSideLobeRatio,
    crossCorrelationPeakLoc,
    chirpPower,
    amplitudeCoefficients,
    phaseCoefficients,
    folder_path,
):
    """
    Create a dataset for chirp thanks to its information

    Parameters
    ----------
    pole: dict[str, list]
        Pole list of values
    ds_attr: dict[str, Any]
        Dictionary of dataset attributes
    replicaQualityValid: dict[str, dict | list]
        ReplicaQualityValid list of values
    crossCorrelationWidth: dict[str, dict | list]
        CrossCorrelationWidth list of values
    sideLobeLevel: dict[str, dict | list]
        SideLobeLevel list of values
    integratedSideLobeRatio: dict[str, dict | list]
        IntegratedSideLobeRatio list of values
    crossCorrelationPeakLoc: dict[str, dict | list]
        CrossCorrelationPeakLoc list of values
    chirpPower: dict[str, dict | list]
        ChirpPower list of values
    amplitudeCoefficients: dict[str, dict | list]
        AmplitudeCoefficients list of values
    phaseCoefficients: dict[str, dict | list]
        PhaseCoefficients list of values
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        Chirp dataset
    """

    ds = xr.Dataset()
    try:
        ds.attrs = ds_attr | generate_doc_ds(
            os.path.join(xpath_dict["doppler"]["xpath"], "chirp"), folder_path
        )
    except Exception:
        pass
    replicaQualityValid_da = xr.DataArray(
        data=replicaQualityValid["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            replicaQualityValid["attr"]
            | generate_doc_vars(replicaQualityValid["attr"]["xpath"], folder_path)
        ),
    )
    crossCorrelationWidth_da = xr.DataArray(
        data=crossCorrelationWidth["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            crossCorrelationWidth["attr"]
            | generate_doc_vars(crossCorrelationWidth["attr"]["xpath"], folder_path)
        ),
    )
    sideLobeLevel_da = xr.DataArray(
        data=sideLobeLevel["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            sideLobeLevel["attr"]
            | generate_doc_vars(sideLobeLevel["attr"]["xpath"], folder_path)
        ),
    )
    integratedSideLobeRatio_da = xr.DataArray(
        data=integratedSideLobeRatio["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            integratedSideLobeRatio["attr"]
            | generate_doc_vars(integratedSideLobeRatio["attr"]["xpath"], folder_path)
        ),
    )
    crossCorrelationPeakLoc_da = xr.DataArray(
        data=crossCorrelationPeakLoc["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            crossCorrelationPeakLoc["attr"]
            | generate_doc_vars(crossCorrelationPeakLoc["attr"]["xpath"], folder_path)
        ),
    )
    chirpPower_da = xr.DataArray(
        data=chirpPower["values"],
        coords={"pole": pole["values"]},
        dims=["pole"],
        attrs=(
            chirpPower["attr"]
            | generate_doc_vars(chirpPower["attr"]["xpath"], folder_path)
        ),
    )
    amplitudeCoefficients_da = xr.DataArray(
        data=np.array(amplitudeCoefficients["values"]),
        coords={
            "pole": pole["values"],
            "n-amplitudeCoefficients": [
                i for i in range(np.array(amplitudeCoefficients["values"]).shape[1])
            ],
        },
        dims=["pole", "n-amplitudeCoefficients"],
        attrs=(
            amplitudeCoefficients["attr"]
            | generate_doc_vars(amplitudeCoefficients["attr"]["xpath"], folder_path)
        ),
    )
    phaseCoefficients_da = xr.DataArray(
        data=np.array(phaseCoefficients["values"]),
        coords={
            "pole": pole["values"],
            "n-phaseCoefficients": [
                i for i in range(np.array(phaseCoefficients["values"]).shape[1])
            ],
        },
        dims=["pole", "n-phaseCoefficients"],
        attrs=(
            phaseCoefficients["attr"]
            | generate_doc_vars(phaseCoefficients["attr"]["xpath"], folder_path)
        ),
    )
    ds["replicaQualityValid"] = replicaQualityValid_da
    ds["crossCorrelationWidth"] = crossCorrelationWidth_da
    ds["sideLobeLevel"] = sideLobeLevel_da
    ds["integratedSideLobeRatio"] = integratedSideLobeRatio_da
    ds["crossCorrelationPeakLoc"] = crossCorrelationPeakLoc_da
    ds["chirpPower"] = chirpPower_da
    ds["amplitudeCoefficients"] = amplitudeCoefficients_da
    ds["phaseCoefficients"] = phaseCoefficients_da
    ds["pole"].attrs = get_type_for_pole(folder_path)
    return ds


def get_dict_radar_parameters(dictio):
    """
    Get the dictionary containing the radar parameters useful information

    Parameters
    ----------
    dictio: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Dictionary with radar parameters information
    """

    xpath = xpath_dict["radarParameters"]["xpath"]
    content_dict = xpath_get(dictio, xpath)
    principal_dic = {"ds_attr": {}}
    ds_attr = list(radar_parameters_key_dict["ds_attributes"])
    vars = list(radar_parameters_key_dict["coord"].keys())
    for var in vars:
        if "incidenceAngleCorrection" in var:
            template_dic = {"noiseLevelValues": [], "coords": {}, "attr": {}}
        else:
            template_dic = {"values": [], "coords": {}, "attr": {}}
        for val in radar_parameters_key_dict["coord"][var]:
            template_dic["coords"][val] = []
        principal_dic[var] = template_dic
    for key in content_dict:
        if key in ds_attr:
            if (key == "polarizations") or (key == "beams"):
                principal_dic["ds_attr"][key] = content_dict[key].split(" ")
            elif isinstance(content_dict[key], dict):
                for intern_key in content_dict[key]:
                    if "@" in intern_key:
                        principal_dic["ds_attr"][
                            f"{key}_{intern_key.replace('@', '')}"
                        ] = parse_value(content_dict[key][intern_key])
                    else:
                        principal_dic["ds_attr"][key] = parse_value(
                            content_dict[key][intern_key]
                        )
            else:
                principal_dic["ds_attr"][key] = parse_value(content_dict[key])
        elif key in vars:
            prefix_path = os.path.join(xpath, key)
            if isinstance(content_dict[key], list):
                principal_dic[key]["attr"]["xpath"] = prefix_path
                for value in content_dict[key]:
                    for intern_key in value:
                        if (
                            intern_key.replace("@", "")
                            in radar_parameters_key_dict["coord"][key]
                        ):
                            principal_dic[key]["coords"][
                                intern_key.replace("@", "")
                            ].append(parse_value(value[intern_key]))
                        elif intern_key == "#text":
                            principal_dic[key]["values"].append(
                                parse_value(value[intern_key])
                            )
                        else:
                            principal_dic[key]["attr"][intern_key.replace("@", "")] = (
                                parse_value(value[intern_key])
                            )
        elif isinstance(content_dict[key], list):
            prefix_path = os.path.join(xpath, key)
            for value in content_dict[key]:
                var_name = ""
                # referenceNoiseLevel case
                for k in value:
                    if "@" in k:
                        var_name = f"{k.replace('@', '')}_{value[k].replace(' ', '_')}"
                    elif isinstance(value[k], str):
                        principal_dic[var_name]["attr"][k] = parse_value(value[k])
                    elif isinstance(value[k], dict):
                        principal_dic[var_name]["attr"]["xpath"] = os.path.join(
                            prefix_path, k
                        )
                        for intern_key in value[k]:
                            if "@" in intern_key:
                                principal_dic[var_name]["attr"][
                                    f"{k}_{intern_key.replace('@', '')}"
                                ] = value[k][intern_key]
                            else:
                                principal_dic[var_name]["noiseLevelValues"] = [
                                    parse_value(x)
                                    for x in value[k][intern_key].split(" ")
                                ]
                                principal_dic[var_name]["coords"][
                                    "NbOfNoiseLevelValues"
                                ] = np.arange(len(value[k][intern_key].split(" ")))

    return principal_dic


def create_dataset_radar_parameters(dictio, folder_path):
    """
    Create a dataset for the radar parameters

    ------------------------------------------------

    :rtype: xarray.Dataset
    :type dictio: dict[str, dict | dict[str, list | TypedDict]]
    :param dictio: Content of the useful dataset information as a dictionary
    :param folder_path: Folder path containing the level 1 files
    :return: Radar parameters dataset
    """
    general_ds = xr.Dataset()
    for key in dictio:
        if key == "ds_attr":
            general_ds.attrs = dictio[key]
        else:
            coords = {}
            data = []
            dims = list(dictio[key]["coords"])
            attr = dictio[key]["attr"]
            if "incidenceAngleCorrection" in key:
                data = dictio[key]["noiseLevelValues"]
                coords = dictio[key]["coords"]
                if "Beta" in key:
                    general_ds["noiseLevelValues_BetaNought"] = xr.DataArray(
                        data=data,
                        dims=dims,
                        attrs=(attr | generate_doc_vars(attr["xpath"], folder_path)),
                    )
                elif "Sigma" in key:
                    general_ds["noiseLevelValues_SigmaNought"] = xr.DataArray(
                        data=data,
                        dims=dims,
                        attrs=(attr | generate_doc_vars(attr["xpath"], folder_path)),
                    )
                elif "Gamma" in key:
                    general_ds["noiseLevelValues_Gamma"] = xr.DataArray(
                        data=data,
                        dims=dims,
                        attrs=(attr | generate_doc_vars(attr["xpath"], folder_path)),
                    )
            else:
                if len(dims) == 2:
                    data = create_2d_matrix(
                        dictio[key]["coords"][dims[0]],
                        dictio[key]["coords"][dims[1]],
                        dictio[key]["values"],
                    )
                    coords[dims[0]] = np.unique(dictio[key]["coords"][dims[0]])
                    coords[dims[1]] = np.unique(dictio[key]["coords"][dims[1]])
                elif len(dims) == 1:
                    data = dictio[key]["values"]
                    coords[dims[0]] = dictio[key]["coords"][dims[0]]
                general_ds[key] = xr.DataArray(
                    data=data,
                    dims=dims,
                    coords=coords,
                    attrs=(attr | generate_doc_vars(attr["xpath"], folder_path)),
                )
    try:
        general_ds.attrs |= generate_doc_ds(
            xpath_dict["radarParameters"]["xpath"], folder_path
        )
    except Exception:
        pass
    return general_ds


def create_2d_matrix(lines, cols, vals):
    """
    Create a matrix with sorted data when it depends on two coords

    Parameters
    ----------
    lines: list
        Lines coords
    cols: list
        Cols coords
    vals: list
        Values corresponding to the data

    Returns
    -------
    numpy.ndarray
        arranged matrix with values
    """

    height = len(np.unique(lines))
    width = len(np.unique(cols))
    tab = np.ones((height, width)) * np.nan
    unique_lines = np.unique(lines)
    unique_cols = np.unique(cols)
    indexs_lines = [np.where(li == unique_lines)[0][0] for li in lines]
    indexs_cols = [np.where(co == unique_cols)[0][0] for co in cols]
    for j in range(len(vals)):
        tab[indexs_lines[j], indexs_cols[j]] = vals[j]
    return tab


def fill_image_attribute(dictio):
    """
    Get raster attributes for image attributes

    Parameters
    ----------
    dictio: dict
        content of product.xml as a dictionary

    Returns
    -------
    dict[str, Any]
        dictionary containing the raster attributes
    """

    xpath = xpath_dict["geolocation_grid"]["xpath"].split("/geographicInformation")[0]
    content_dict = xpath_get(dictio, xpath)
    attr = {}
    for key in content_dict:
        if isinstance(content_dict[key], str):
            attr[key] = parse_value(content_dict[key])
        elif key == "rasterAttributes":
            for value in content_dict[key]:
                if isinstance(content_dict[key][value], str):
                    attr[f"{key}_{value}"] = parse_value(content_dict[key][value])
                elif isinstance(content_dict[key][value], dict):
                    dico_keys = list(content_dict[key][value].keys())
                    for k in dico_keys:
                        attr[
                            f"{key}_{value}_{k.replace('@', '').replace('#text', 'value')}"
                        ] = parse_value(content_dict[key][value][k])
    return attr


def parse_value(value):
    """
    Parse a value to return it in the appropriated type (except dates)

    Parameters
    ----------
    value: String
        value as a string to parse (type undefined)

    Returns
    -------
    Any
        typed value if the type is recognized, else the input string
    """

    try:
        return ast.literal_eval(value)
    except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
        return value


def list_xsd_files(xpath, folder_path):
    """
    Search for xsd file paths where names are linked with a dataset, thanks to xpaths

    Parameters
    ----------
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parents
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    list[list[str]]
        List of xsd file paths linked to a dataset
    """

    xsd_folder = os.path.join(folder_path, "schemas")
    var_name = xpath.split("/")[-1]
    parent_path = os.path.dirname(xpath)
    parent_var_name = parent_path.split("/")[-1]
    grand_parent_path = os.path.dirname(parent_path)
    grand_parent_var_name = grand_parent_path.split("/")[-1]
    list_files = glob.glob(os.path.join(xsd_folder, "*"))

    interesting_files = [
        [file for file in list_files if (var_name in os.path.basename(file))],
        [file for file in list_files if (parent_var_name in os.path.basename(file))],
        [
            file
            for file in list_files
            if (grand_parent_var_name in os.path.basename(file))
        ],
    ]
    return interesting_files


def find_doc_in_xsd_files(xpath, interesting_files):
    """
    Search dataset variable description in a list of xsd files

    Parameters
    ----------
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parents
    interesting_files: list[list[str]]

    Returns
    -------
    dict[str, str]
        Dictionary containing a dataset variable description
    """

    var_name = xpath.split("/")[-1]
    xsd_path = "xsd:schema/xsd:complexType/xsd:sequence/xsd:element"
    description = ""
    for element in interesting_files:
        # index = interesting_files.index(element)
        if len(element) != 0:
            for pathname in element:
                with open(pathname, "rb") as f:
                    xml_content = f.read()
                    dic = xmltodict.parse(xml_content)
                    f.close()
                content_dict = xpath_get(dic, xsd_path)
                for value in content_dict:
                    if value["@name"] == var_name:
                        description = (
                            value["xsd:annotation"]["xsd:documentation"]
                            .replace("  ", "")
                            .replace("\n", "")
                            .replace("\t", " ")
                        )
    return {"Description": description}


def get_type_for_pole(folder_path):
    """
    get the type for a polarization

    Parameters
    ----------
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    dict[str, str]
        Dictionary containing the xsd filename of the polarization description
    """

    pathname = os.path.join(folder_path, "schemas", "rs2prod_chirp.xsd")
    xpath = "/xsd:schema/xsd:complexType/xsd:attribute"
    with open(pathname, "rb") as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()
    content_dict = xpath_get(dic, xpath)
    for values in content_dict:
        if values["@name"] == "pole":
            return {"type": values["@type"]}


def find_doc_for_ds_in_xsd_files(xpath, interesting_files):
    """
    Search dataset general description in a list of xsd files

    Parameters
    ----------
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parents
    interesting_files: list[list[str]]
        List of files paths containing the dataset name, its parents name or grandparents names

    Returns
    -------
    dict[str, str]
        Dictionary containing the dataset general description
    """

    if "geolocationGrid" in xpath:
        interesting_files = interesting_files[1]
    else:
        interesting_files = interesting_files[0]
    xsd_path = "xsd:schema/xsd:annotation/xsd:documentation"
    for element in interesting_files:
        # index = interesting_files.index(element)
        if len(element) != 0:
            with open(element, "rb") as f:
                xml_content = f.read()
                dic = xmltodict.parse(xml_content)
                f.close()
            content_dict = (
                xpath_get(dic, xsd_path)
                .replace("  ", "")
                .replace("\n", "")
                .replace("\t", " ")
            )
    return {"Description": content_dict}


def generate_doc_vars(xpath, folder_path):
    """
    Generate a Dataset variable documentation found in xsd files from a folder

    Parameters
    ----------
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parents
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    dict[str, str]
        Dataset variable documentation
    """

    return find_doc_in_xsd_files(xpath, list_xsd_files(xpath, folder_path))


def generate_doc_ds(xpath, folder_path):
    """
    Generate a Dataset general documentation found in xsd files from a folder

    Parameters
    ----------
    xpath: str
        Xpath describing the hierarchy of a dataset to get its name and its parent
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    dict[str, str]
        Dataset general documentation
    """

    return find_doc_for_ds_in_xsd_files(xpath, list_xsd_files(xpath, folder_path))


def list_lut_files(folder_path):
    """
    Return a list of LookUpTable files path names present in a folder

    Parameters
    ----------
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    list[str]
        List of LookUpTable file paths located in a folder
    """

    return glob.glob(os.path.join(folder_path, "lut*xml"))


def create_data_array_lut(dictio, dt):
    """
    Return a dataArray for a single LookUpTable file

    Parameters
    ----------
    dictio: dict
        Content of a xml LookUpTable file described as a dictionary
    dt:xr.Datatree
        Contains every dataset of product.xml

    Returns
    -------
    xarray.DataArray
        DataArray describing a lookUpTable
    """

    final_lut_dict = {"attrs": {}}
    for value in dictio["lut"]:
        if "@" not in value:
            if value != "gains":
                final_lut_dict["attrs"][value] = parse_value(dictio["lut"][value])
            else:
                final_lut_dict[value] = [
                    parse_value(x) for x in dictio["lut"][value].split(" ")
                ]
                assert dt["geolocationGrid"].attrs[
                    "rasterAttributes_numberOfSamplesPerLine"
                ] == len(final_lut_dict[value])
    da = xr.DataArray(
        data=final_lut_dict["gains"],
        dims=["pixel"],
        coords={
            "pixel": np.arange(
                dt["geolocationGrid"].attrs["rasterAttributes_numberOfSamplesPerLine"]
            )
        },
    )
    da.attrs = final_lut_dict["attrs"]
    return da


def create_dataset_lut(files, dt, folder_path):
    """
    Return a dataset that contains LookUpTables

    Parameters
    ----------
    files: List[str]
        Path names of LookUpTables files in the current folder
    dt: xr.Datatree
        Contains every dataset of product.xml
    folder_path: str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Dataset
        LUT dataset
    """

    ds = xr.Dataset()
    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        with open(file, "rb") as f:
            xml_content = f.read()
            dic = xmltodict.parse(xml_content)
            f.close()
        ds[filename] = create_data_array_lut(dic, dt)
    try:
        ds.attrs = find_doc_for_ds_in_xsd_files(
            "lut", list_xsd_files("lut", folder_path)
        )
    except Exception:
        pass
    return ds


def sort_list_files_and_get_pols(list_tiff):
    """
    From a list of tiff files, sort it to get the co polarization tiff file as the first element, and extract pols

    Parameters
    ----------
    list_tiff: List[str]
        List of tiff files

    Returns
    -------
    (List[str], List[str])
        Tuple that contains the tiff files list sorted and the polarizations
    """
    pols = []
    if len(list_tiff) > 1:
        first_base = os.path.basename(list_tiff[0]).split(".")[0]
        first_pol = first_base[-2:]
        if first_pol[0] != first_pol[1]:
            list_tiff.reverse()
    for file in list_tiff:
        base = os.path.basename(file).split(".")[0]
        pol = base[-2:]
        pols.append(pol)
    return list_tiff, pols


def list_tiff_files(root_path):
    """
    Return a list that contains all tiff files paths

    Parameters
    ----------
    root_path: str
        Folder path containing the level 1 files

    Returns
    -------
    List[str]
        List of Tiff file paths located in a folder
    """

    return glob.glob(os.path.join(root_path, "*imagery*tif"))


def get_glob(strlist):
    """
    From list of str, replace diff by '?'

    Parameters
    ----------
    strlist: List[str]
        List of str

    Returns
    -------
    str
        replace diff between strings by `?`
    """

    def _get_glob(st):
        stglob = "".join(
            [
                "?" if len(charlist) > 1 else charlist[0]
                for charlist in [list(set(charset)) for charset in zip(*st)]
            ]
        )
        return re.sub(r"\?+", "*", stglob)

    strglob = _get_glob(strlist)
    if strglob.endswith("*"):
        strglob += _get_glob(s[::-1] for s in strlist)[::-1]
        strglob = strglob.replace("**", "*")

    return strglob


def load_digital_number(
    dt,
    resolution=None,
    chunks={"sample": 5000, "line": 5000},
    resampling=rasterio.enums.Resampling.rms,
):
    """
    Add to the initial datatree (that contains all datasets) a new dataset with digital numbers
    (name of the dataset = `digital_numbers`)

    Parameters
    ----------
    dt: xr.Datatree
        datatree containing every dataset
    resolution: str, dict[str, int], None or number
        Resampling dict like `{'line': 20, 'sample': 20}` where 20 is in pixels.
        if a number, dict will be constructed from `{'line': number, 'sample': number}`
        if str, it must end with 'm' (meters), like '100m'. dict will be computed from sensor pixel size
    chunks: dict[str, int]
        dict with keys ['pol','line','sample'] (dask chunks)
    resampling: rasterio.enums.Resampling
        Only used if `resolution` is not None.
        ` rasterio.enums.Resampling.rms` by default. `rasterio.enums.Resampling.nearest` (decimation) is fastest

    Returns
    -------
    xr.Datatree
        Initial datatree + dataset (possibly dual-pol), with basic coords/dims naming convention
    """

    # This func permit to calculate the resampling with Gaussian filter and xtiling (for the posting)
    """def resamp_signal(tiff_path, pol, resolution):
        im = Image.open(tiff_path)
        imarray = np.asarray(im)
        imarray = imarray[np.newaxis, :, :]
        res = scipy.ndimage.gaussian_filter(imarray, sigma=0.5 * resolution['sample'], mode='mirror', truncate=4)
        res_da = xr.DataArray(res, coords={'pol': [pol], 'line': np.arange(imarray.shape[1]),
                                           'sample': np.arange(imarray.shape[2])}, dims=['pol', 'line', 'sample'])
        bb_da = xr.DataArray(imarray, coords={'pol': [pol], 'line': np.arange(imarray.shape[1]),
                                              'sample': np.arange(imarray.shape[2])}, dims=['pol', 'line', 'sample'])
        nperseg = {'sample': int(resolution['sample']), 'line': int(resolution['line'])}
        posting = xtiling(bb_da, nperseg, noverlap=0, centering=False, side='left', prefix='')
        resampled_signal = res_da.isel(sample=posting['sample'].sample, line=posting['line'].line)
        return resampled_signal"""

    tiff_files = list_tiff_files(dt.attrs["product_path"])
    tiff_files, pols = sort_list_files_and_get_pols(tiff_files)
    map_dims = {"pol": "band", "line": "y", "sample": "x"}
    if resolution is not None:
        comment = f'resampled at "{resolution}" with {resampling.__module__}.{resampling.__class__.__name__}.{resampling.name}'
    else:
        comment = "read at full resolution"

    # arbitrary rio object, to get shape, etc ... (will not be used to read data)
    try:
        rio = rasterio.open(tiff_files[0])
    except rasterio.RasterioIOError:
        print(traceback.format_exc())
        raise rasterio.RasterioIOError(
            f"{os.path.basename(tiff_files[0])} couldn't be opened"
        )

    chunks["pol"] = 1
    # sort chunks keys like map_dims
    chunks = dict(
        sorted(chunks.items(), key=lambda pair: list(map_dims.keys()).index(pair[0]))
    )
    chunks_rio = {map_dims[d]: chunks[d] for d in map_dims.keys()}
    if resolution is None:
        # using tiff driver: need to read individual tiff and concat them
        # riofiles['rio'] is ordered like self.s1meta.manifest_attrs['polarizations']

        dn = xr.concat(
            [
                rioxarray.open_rasterio(f, chunks=chunks_rio, parse_coordinates=False)
                for f in tiff_files
            ],
            "band",
        ).assign_coords(band=np.arange(len(pols)) + 1)

        # set dimensions names
        dn = dn.rename(dict(zip(map_dims.values(), map_dims.keys())))

        # create coordinates from dimension index (because of parse_coordinates=False)
        dn = dn.assign_coords({"line": dn.line, "sample": dn.sample})
        dn = dn.drop_vars("spatial_ref", errors="ignore")
    else:
        if not isinstance(resolution, dict):
            if isinstance(resolution, str) and resolution.endswith("m"):
                resolution = float(resolution[:-1])
            resolution = dict(
                line=resolution
                / dt["geolocationGrid"]["line"].attrs[
                    "rasterAttributes_sampledLineSpacing_value"
                ],
                sample=resolution
                / dt["geolocationGrid"]["pixel"].attrs[
                    "rasterAttributes_sampledPixelSpacing_value"
                ],
            )

        # resample the DN at gdal level, before feeding it to the dataset
        out_shape = (
            int(rio.height / resolution["line"]),
            int(rio.width / resolution["sample"]),
        )
        out_shape_pol = (1,) + out_shape
        # read resampled array in one chunk, and rechunk
        # this doesn't optimize memory, but total size remain quite small

        if isinstance(resolution["line"], int):
            # legacy behaviour: winsize is the maximum full image size that can be divided  by resolution (int)
            winsize = (
                0,
                0,
                rio.width // resolution["sample"] * resolution["sample"],
                rio.height // resolution["line"] * resolution["line"],
            )
            window = rasterio.windows.Window(*winsize)
        else:
            window = None

        dn = xr.concat(
            [
                xr.DataArray(
                    dask.array.from_array(
                        rasterio.open(f).read(
                            out_shape=out_shape_pol,
                            resampling=resampling,
                            window=window,
                        ),
                        chunks=chunks_rio,
                    ),
                    dims=tuple(map_dims.keys()),
                    coords={"pol": [pol]},
                )
                for f, pol in zip(tiff_files, pols)
            ],
            "pol",
        ).chunk(chunks)
        # create coordinates at box center
        translate = Affine.translation(
            (resolution["sample"] - 1) / 2, (resolution["line"] - 1) / 2
        )
        scale = Affine.scale(
            rio.width // resolution["sample"] * resolution["sample"] / out_shape[1],
            rio.height // resolution["line"] * resolution["line"] / out_shape[0],
        )
        sample, _ = translate * scale * (dn.sample, 0)
        _, line = translate * scale * (0, dn.line)
        dn = dn.assign_coords({"line": line, "sample": sample})

    # for GTiff driver, pols are already ordered. just rename them
    dn = dn.assign_coords(pol=pols)
    # This commented part of code was used to resample without rasterio (gaussian filter)
    """ list_xr = []
        for f, pol in zip(tiff_files, pols):
            resampled_da = resamp_signal(f, pol, resolution)
            list_xr.append(xr.DataArray(
                dask.array.from_array(
                    resampled_da,
                    chunks=chunks_rio,
                ),
                dims=tuple(map_dims.keys()),
                coords={"pol": [pol], "line": resampled_da.line, "sample": resampled_da.sample},
            ))
        dn = xr.concat(list_xr, 'pol').chunk(chunks)
    dn = dn.assign_coords(pol=pols)"""
    descr = "not denoised"
    var_name = "digital_number"

    dn.attrs = {
        "comment": f"{descr} digital number, {comment}",
        "history": yaml.safe_dump(
            {
                var_name: get_glob(
                    [p.replace(dt.attrs["product_path"] + "/", "") for p in tiff_files]
                )
            }
        ),
    }
    ds = dn.to_dataset(name=var_name)
    dt["digital_numbers"] = xr.DataTree(ds)
    return dt


def get_product_attributes(dic):
    """
    Return a dictionary with every important general attributes from product.xml

    Parameters
    ----------
    dic: dict
        content of product.xml as a dictionary

    Returns
    -------
    dict
        General attributes from product.xml
    """

    dic = xpath_get(dic, "/product/sourceAttributes")
    final_dic = dict()
    useful_attributes = ["satellite", "inputDatasetId", "rawDataStartTime"]
    for key in dic:
        if key == "radarParameters":
            break
        else:
            if isinstance(dic[key], str) and key in useful_attributes:
                if "Time" in key:
                    dt = datetime.strptime(dic[key], "%Y-%m-%dT%H:%M:%S.%fZ")
                    final_dic[key] = np.array(dt).astype("datetime64[ns]")
                else:
                    final_dic[key] = dic[key]
    return final_dic


def get_satellite_height(dic):
    """
    Find information of satellite height in product.xml

    Parameters
    ----------
    dic: dict
        Content of product.xml as a dictionary

    Returns
    -------
    dict
        Contains value and unit of satellite height
    """
    xpath = (
        "/product/imageGenerationParameters/sarProcessingInformation/satelliteHeight"
    )
    dictio = xpath_get(dic, xpath)
    name = "satelliteHeight"
    content_dic = dict()
    for key in dictio:
        if "@" in key:
            content_dic[f"{name}_{key.replace('@', '')}"] = dictio[key]
        else:
            content_dic[name] = parse_value(dictio[key])
    return content_dic


def get_satellite_pass_direction(dic):
    """
    Get the pass direction of the satellite (ascending / descending)

    Parameters
    ----------
    dic : dict
        Content of product.xml as a dictionary
    Returns
    -------
    dict
        Contains the pass direction
    """
    xpath = "/product/sourceAttributes/orbitAndAttitude/orbitInformation/passDirection"
    value = xpath_get(dic, xpath)
    name = "passDirection"
    return {name: value}


def verify_hierarchy_product(folder_path):
    """
    Verify if the product includes the necessary files

    Parameters
    ----------
    folder_path : str
        Folder path containing the level 1 files

    Returns
    -------
    Nothing, or an error telling that it is missing files (and specify the concerned one)
    """
    must_have_files = [
        "product.xml",
        "lutBeta.xml",
        "lutGamma.xml",
        "lutSigma.xml",
        "imagery_*.tif",
    ]
    for file in must_have_files:
        entire_file_path = os.path.join(folder_path, file)
        list_file = glob.glob(entire_file_path)
        if len(list_file) == 0:
            raise FileNotFoundError(
                f"The file {file} is missing for the product {folder_path}"
            )


def rs2_reader(folder_path):
    """
    Principal function of the reader, that create a datatree with all the product.xml and lut xml files dataset

    Parameters
    ----------
    folder_path : str
        Folder path containing the level 1 files

    Returns
    -------
    xarray.Datatree
        datatree containing every dataset
    """
    # Verify if the product is complete
    verify_hierarchy_product(folder_path)

    # get product.xml path
    product_xml_path = os.path.join(folder_path, "product.xml")

    # get product.xml content as a dict
    with open(product_xml_path, "rb") as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()

    # Create dictionaries then dataset, and fill them in a datatree
    dic_geo = create_dic_geolocation_grid(dic)
    ds_geo = create_dataset_geolocation_grid(dic_geo, folder_path)
    ds_geo.attrs = fill_image_attribute(dic)
    dic_orbit_information = get_dic_orbit_information(dic)
    ds_orbit_info = create_dataset_orbit_information(
        dic_orbit_information["ds_attr"],
        dic_orbit_information["timestamp"],
        dic_orbit_information["xPosition"],
        dic_orbit_information["yPosition"],
        dic_orbit_information["zPosition"],
        dic_orbit_information["xVelocity"],
        dic_orbit_information["yVelocity"],
        dic_orbit_information["zVelocity"],
        folder_path,
    )
    dic_attitude_info = get_dic_attitude_info(dic)
    ds_attitude_info = create_dataset_attitude_information(
        dic_attitude_info["ds_attr"],
        dic_attitude_info["timestamp"],
        dic_attitude_info["yaw"],
        dic_attitude_info["roll"],
        dic_attitude_info["pitch"],
        folder_path,
    )
    ds_orbit_attitude_info = xr.merge([ds_attitude_info, ds_orbit_info])
    try:
        ds_orbit_attitude_info.attrs["Description"] += (
            ". " + ds_orbit_info.attrs["Description"] + "."
        )
    except Exception:
        pass
    dt = xr.DataTree()
    dt["orbitAndAttitude"] = xr.DataTree(ds_orbit_attitude_info)
    dt["geolocationGrid"] = xr.DataTree(ds_geo)
    dic_doppler_centroid = get_dict_doppler_centroid(dic)
    ds_doppler_centroid = create_dataset_doppler_centroid(
        dic_doppler_centroid["ds_attr"],
        dic_doppler_centroid["timeOfDopplerCentroidEstimate"],
        dic_doppler_centroid["dopplerAmbiguity"],
        dic_doppler_centroid["dopplerAmbiguityConfidence"],
        dic_doppler_centroid["dopplerCentroidReferenceTime"],
        dic_doppler_centroid["dopplerCentroidPolynomialPeriod"],
        dic_doppler_centroid["dopplerCentroidCoefficients"],
        dic_doppler_centroid["dopplerCentroidConfidence"],
        folder_path,
    )
    dt["imageGenerationParameters/doppler/dopplerCentroid"] = xr.DataTree(
        ds_doppler_centroid
    )
    dic_doppler_rate_values = get_dic_doppler_rate_values(dic)
    ds_doppler_rate_values = create_dataset_doppler_rate_values(
        dic_doppler_rate_values["ds_attr"],
        dic_doppler_rate_values["dopplerRateReferenceTime"],
        dic_doppler_rate_values["dopplerRateValuesCoefficients"],
        folder_path,
    )
    dt["imageGenerationParameters/doppler/dopplerRateValues"] = xr.DataTree(
        ds_doppler_rate_values
    )
    dic_chirp = get_dict_chirp(dic)
    ds_chirp = create_dataset_chirp(
        dic_chirp["pole"],
        dic_chirp["ds_attr"],
        dic_chirp["replicaQualityValid"],
        dic_chirp["crossCorrelationWidth"],
        dic_chirp["sideLobeLevel"],
        dic_chirp["integratedSideLobeRatio"],
        dic_chirp["crossCorrelationPeakLoc"],
        dic_chirp["chirpPower"],
        dic_chirp["amplitudeCoefficients"],
        dic_chirp["phaseCoefficients"],
        folder_path,
    )
    dt["imageGenerationParameters/chirp"] = xr.DataTree(ds_chirp)
    radar_parameters_dic = get_dict_radar_parameters(dic)
    ds_radar_parameters = create_dataset_radar_parameters(
        radar_parameters_dic, folder_path
    )
    dt["radarParameters"] = xr.DataTree(ds_radar_parameters)
    ds_lut = create_dataset_lut(list_lut_files(folder_path), dt, folder_path)
    dt["lut"] = xr.DataTree(ds_lut)
    dt.attrs["product_path"] = folder_path
    dt.attrs |= get_product_attributes(dic)
    dt.attrs |= get_satellite_height(dic)
    dt.attrs |= get_satellite_pass_direction(dic)
    return dt
