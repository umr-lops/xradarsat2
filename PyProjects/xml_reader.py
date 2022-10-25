# imports
import numpy as np
from lxml import etree
from lxml import objectify
import xmltodict
import xarray as xr
import time
import xsar
import datatree
import os
import glob

xpath_dict = {
    "geolocation_grid":
        {
            "xpath": "/product/imageAttributes/geographicInformation/geolocationGrid/imageTiePoint"
        },
    "orbit_information":
        {
            "xpath": "/product/sourceAttributes/orbitAndAttitude/orbitInformation",
        },
    "attitude_information":
        {
            "xpath": "/product/sourceAttributes/orbitAndAttitude/attitudeInformation"
        },
    "doppler":
        {
            "xpath": "/product/imageGenerationParameters"
        },
    "radarParameters": {
        "xpath": "/product/sourceAttributes/radarParameters"
    }
}

radar_parameters_key_dict = {
    "ds_attributes": ["acquisitionType", "beams", "polarizations", "pulses", "radarCenterFrequency", "antennaPointing",
                      "yawSteeringFlag", "geodeticFlag", "rawBitsPerSample", "pulseLength", "pulseBandwidth",
                      "adcSamplingRate"],
    "coord": {
        "pulsesReceivedPerDwell": ["beam"],
        "numberOfPulseIntervalsPerDwell": ["beam"],
        "rank": ["beam"],
        "settableGain": ["beam", "pole"],
        "pulseRepetitionFrequency": ["beam"],
        "samplesPerEchoLine": ["beam"],
        "incidenceAngleCorrection_Beta_Nought": [],
        "incidenceAngleCorrection_Sigma_Nought": [],
        "incidenceAngleCorrection_Gamma": []
    }

}


def xpath_get(mydict, xpath):
    elem = mydict
    try:
        for x in xpath.strip("/").split("/"):
            elem = elem.get(x)
    except:
        pass
    return elem


def create_dic_geolocation_grid(dictio):
    content_list = xpath_get(dictio, xpath_dict["geolocation_grid"]["xpath"])
    final_dc = {
        "latitude": {
            "values": [],
            "units": ""
        },
        "longitude": {
            "values": [],
            "units": ""
        },
        "height": {
            "values": [],
            "units": ""
        },
        "coords": {
            "lines": [],
            "pixels": []
        },
        "attr": get_line_and_pix_info(fill_image_attribute(dictio))
    }
    for element in content_list:
        final_dc["coords"]["lines"].append(parse_value(element['imageCoordinate']['line']))
        final_dc["coords"]["pixels"].append(parse_value(element['imageCoordinate']['pixel']))
        final_dc["longitude"]["values"].append(parse_value(element['geodeticCoordinate']['longitude']['#text']))
        final_dc["latitude"]["values"].append(parse_value(element['geodeticCoordinate']['latitude']['#text']))
        final_dc["height"]["values"].append(parse_value(element['geodeticCoordinate']['height']['#text']))
        final_dc["longitude"]["units"] = element['geodeticCoordinate']['longitude']['@units']
        final_dc["latitude"]["units"] = element['geodeticCoordinate']['latitude']['@units']
        final_dc["height"]["units"] = element['geodeticCoordinate']['height']['@units']
    return final_dc


def create_dataset_geolocation_grid(dictio, folder_path):
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
            xpath = os.path.join(xpath_dict['geolocation_grid']['xpath'], xpath_suffix)
            da = xr.DataArray(data=data, name=key, coords={"line": np.unique(np.array(lines)), "pixel":
                np.unique(np.array(pixs))}, dims=['line', "pixel"], attrs=({"units": unit, "xpath": xpath} |
                                                                           generate_doc_vars(xpath, folder_path)))
            ds[key] = da
    ds.attrs = generate_doc_ds(xpath_dict['geolocation_grid']['xpath'], folder_path)
    ds["line"].attrs = dictio["attr"]["line"]
    ds["pixel"].attrs = dictio["attr"]["pixel"]
    return ds


def get_line_and_pix_info(dictio):
    line = {}
    pixel = {}
    for key in dictio:
        if "PixelSpacing" in key:
            pixel[key] = dictio[key]
        if "LineSpacing" in key:
            line[key] = dictio[key]
    return {
        "line": line,
        "pixel": pixel
    }


def get_dic_orbit_information(dictio):
    content_list = xpath_get(dictio, xpath_dict["orbit_information"]["xpath"])
    ds_attr = {}
    timestamp = []
    xPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/xPosition"
        }
    }
    yPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/yPosition"
        }
    }
    zPosition = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/zPosition"
        }
    }
    xVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/xVelocity"
        }
    }
    yVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/yVelocity"
        }
    }
    zVelocity = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['geolocation_grid']['xpath']}/stateVector/zVelocity"
        }
    }
    for key in content_list:
        if isinstance(content_list[key], str):
            ds_attr[key] = content_list[key]
        elif isinstance(content_list[key], list):
            """for index in range(len(content_list[key])):
                value = content_list[key][index]
                timestamp.append(content_list[key][index]["timeStamp"])
                xPosition["values"].append(content_list[key][index]["xPosition"]["#text"])
                xPosition["unit"] = content_list[key][index]["xPosition"]["@units"]"""
            for value in content_list[key]:
                timestamp.append(np.datetime64(value["timeStamp"]))
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
        "zVelocity": zVelocity
    }


def create_dataset_orbit_information(ds_attr, timestamp, xPos, yPos, zPos, xVel, yVel, zVel, folder_path):
    ds = xr.Dataset()
    xpos_da = xr.DataArray(data=xPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(xPos["attr"] | generate_doc_vars(xPos["attr"]["xpath"], folder_path)))
    ypos_da = xr.DataArray(data=yPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(yPos["attr"] | generate_doc_vars(yPos["attr"]["xpath"], folder_path)))
    zpos_da = xr.DataArray(data=zPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(zPos["attr"] | generate_doc_vars(zPos["attr"]["xpath"], folder_path)))
    xvel_da = xr.DataArray(data=xVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(xVel["attr"] | generate_doc_vars(xVel["attr"]["xpath"], folder_path)))
    yvel_da = xr.DataArray(data=yVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(yVel["attr"] | generate_doc_vars(yVel["attr"]["xpath"], folder_path)))
    zvel_da = xr.DataArray(data=zVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(zVel["attr"] | generate_doc_vars(zVel["attr"]["xpath"], folder_path)))
    ds["xPosition"] = xpos_da
    ds["yPosition"] = ypos_da
    ds["zPosition"] = zpos_da
    ds["xVelocity"] = xvel_da
    ds["yVelocity"] = yvel_da
    ds["zVelocity"] = zvel_da
    ds.attrs = (ds_attr | generate_doc_ds(xpath_dict["orbit_information"]["xpath"], folder_path))
    ds.attrs["Description"] += ". " + generate_doc_ds('stateVector', folder_path)["Description"]
    return ds


def get_dic_attitude_info(dictio):
    content_list = xpath_get(dictio, xpath_dict["attitude_information"]["xpath"])
    ds_attr = {}
    timestamp = []
    yaw = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['attitude_information']['xpath']}/attitudeAngles/yaw"
        }
    }
    roll = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['attitude_information']['xpath']}/attitudeAngles/roll"
        }
    }
    pitch = {
        "values": [],
        "attr": {
            "units": "",
            "xpath": f"{xpath_dict['attitude_information']['xpath']}/attitudeAngles/pitch"
        }
    }
    for key in content_list:
        if isinstance(content_list[key], str):
            ds_attr[key] = content_list[key]
        elif isinstance(content_list[key], list):
            for value in content_list[key]:
                timestamp.append(np.datetime64(value["timeStamp"]))
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
        "pitch": pitch
    }


def create_dataset_attitude_information(ds_attr, timestamp, yaw, roll, pitch, folder_path):
    ds = xr.Dataset()
    interesting_files_yaw = list_xsd_files(yaw["attr"]["xpath"], folder_path)
    yaw_doc = find_doc_in_xsd_files(yaw["attr"]["xpath"], interesting_files_yaw)
    interesting_files_roll = list_xsd_files(roll["attr"]["xpath"], folder_path)
    roll_doc = find_doc_in_xsd_files(roll["attr"]["xpath"], interesting_files_roll)
    interesting_files_pitch = list_xsd_files(pitch["attr"]["xpath"], folder_path)
    pitch_doc = find_doc_in_xsd_files(pitch["attr"]["xpath"], interesting_files_pitch)
    yaw_da = xr.DataArray(data=yaw["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                          attrs=(yaw["attr"] | yaw_doc))
    roll_da = xr.DataArray(data=roll["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(roll["attr"] | roll_doc))
    pitch_da = xr.DataArray(data=pitch["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                            attrs=(pitch["attr"] | pitch_doc))
    ds["yaw"] = yaw_da
    ds["roll"] = roll_da
    ds["pitch"] = pitch_da
    ds.attrs = (ds_attr | generate_doc_ds(xpath_dict["attitude_information"]["xpath"], folder_path))
    return ds


def get_dict_doppler_centroid(dictio):
    content_list = xpath_get(dictio, xpath_dict["doppler"]["xpath"])
    ds_attr = {}
    times = []
    Ambiguity = {
        "values": [],
        "attr": {}
    }
    AmbiguityConfidence = {
        "values": [],
        "attr": {}
    }
    CentroidReferenceTime = {
        "values": [],
        "attr": {}
    }
    CentroidPolynomialPeriod = {
        "values": [],
        "attr": {}
    }
    CentroidCoefficients = {
        "values": [],
        "attr": {}
    }
    CentroidConfidence = {
        "values": [],
        "attr": {}
    }
    for key in content_list:
        xpath = xpath_dict["doppler"]["xpath"]
        if key == "dopplerCentroid":
            xpath += f"/{key}"
            for value in content_list[key]:
                times.append(np.datetime64(value["timeOfDopplerCentroidEstimate"]))
                Ambiguity["values"].append(int(value["dopplerAmbiguity"]))
                Ambiguity["attr"]["xpath"] = xpath + "/dopplerAmbiguity"
                AmbiguityConfidence["values"].append(float(value["dopplerAmbiguityConfidence"]))
                AmbiguityConfidence["attr"]["xpath"] = xpath + "/dopplerAmbiguityConfidence"
                CentroidReferenceTime["values"].append(float(value["dopplerCentroidReferenceTime"]["#text"]))
                CentroidReferenceTime["attr"]["units"] = value["dopplerCentroidReferenceTime"]["@units"]
                CentroidReferenceTime["attr"]["xpath"] = xpath + "/dopplerCentroidReferenceTime"
                CentroidPolynomialPeriod["values"].append(float(value["dopplerCentroidPolynomialPeriod"]["#text"]))
                CentroidPolynomialPeriod["attr"]["units"] = value["dopplerCentroidPolynomialPeriod"]["@units"]
                CentroidPolynomialPeriod["attr"]["xpath"] = xpath + "/dopplerCentroidPolynomialPeriod"
                CentroidCoefficients["values"].append(
                    [float(x) for x in value["dopplerCentroidCoefficients"].split(" ")])
                CentroidCoefficients["attr"]["xpath"] = xpath + "/dopplerCentroidCoefficients"
                CentroidConfidence["values"].append(float(value["dopplerCentroidConfidence"]))
                CentroidConfidence["attr"]["xpath"] = xpath + "/dopplerCentroidConfidence"
        elif len(times) != 0:
            break
    return {
        "ds_attr": ds_attr,
        "timeOfDopplerCentroidEstimate": times,
        "dopplerAmbiguity": Ambiguity,
        "dopplerAmbiguityConfidence": AmbiguityConfidence,
        "dopplerCentroidReferenceTime": CentroidReferenceTime,
        "dopplerCentroidPolynomialPeriod": CentroidPolynomialPeriod,
        "dopplerCentroidCoefficients": CentroidCoefficients,
        "dopplerCentroidConfidence": CentroidConfidence
    }


def create_dataset_doppler_centroid(ds_attr, times, Ambiguity, AmbiguityConfidence, CentroidReferenceTime,
                                    CentroidPolynomialPeriod, CentroidCoefficients, CentroidConfidence, folder_path):
    ds = xr.Dataset()
    Ambiguity_doc = find_doc_in_xsd_files(Ambiguity["attr"]["xpath"], list_xsd_files(Ambiguity["attr"]["xpath"],
                                                                                     folder_path))
    AmbiguityConfidence_doc = find_doc_in_xsd_files(AmbiguityConfidence["attr"]["xpath"],
                                                    list_xsd_files(AmbiguityConfidence["attr"]["xpath"],
                                                                   folder_path))
    CentroidReferenceTime_doc = find_doc_in_xsd_files(CentroidReferenceTime["attr"]["xpath"],
                                                      list_xsd_files(CentroidReferenceTime["attr"]["xpath"],
                                                                     folder_path))
    CentroidPolynomialPeriod_doc = find_doc_in_xsd_files(CentroidPolynomialPeriod["attr"]["xpath"],
                                                         list_xsd_files(CentroidPolynomialPeriod["attr"]["xpath"],
                                                                        folder_path))
    CentroidCoefficients_doc = find_doc_in_xsd_files(CentroidCoefficients["attr"]["xpath"],
                                                     list_xsd_files(CentroidCoefficients["attr"]["xpath"], folder_path))
    CentroidConfidence_doc = find_doc_in_xsd_files(CentroidConfidence["attr"]["xpath"],
                                                   list_xsd_files(CentroidConfidence["attr"]["xpath"], folder_path))
    ambiguity_da = xr.DataArray(data=Ambiguity["values"], coords={"timeOfDopplerCentroidEstimate": times},
                                dims=["timeOfDopplerCentroidEstimate"],
                                attrs=(Ambiguity["attr"] | Ambiguity_doc))
    ambiguityConfidence_da = xr.DataArray(data=AmbiguityConfidence["values"],
                                          coords={"timeOfDopplerCentroidEstimate": times},
                                          dims=["timeOfDopplerCentroidEstimate"],
                                          attrs=(AmbiguityConfidence["attr"] | AmbiguityConfidence_doc))
    centroidReferenceTime_da = xr.DataArray(data=CentroidReferenceTime["values"],
                                            coords={"timeOfDopplerCentroidEstimate": times},
                                            dims=["timeOfDopplerCentroidEstimate"],
                                            attrs=(CentroidReferenceTime["attr"] | CentroidReferenceTime_doc))
    centroidPolynomialPeriod_da = xr.DataArray(data=CentroidPolynomialPeriod["values"],
                                               coords={"timeOfDopplerCentroidEstimate": times},
                                               dims=["timeOfDopplerCentroidEstimate"],
                                               attrs=(CentroidPolynomialPeriod["attr"] | CentroidPolynomialPeriod_doc))
    centroidCoefficients_da = xr.DataArray(data=np.array(CentroidCoefficients["values"]),
                                           coords={"timeOfDopplerCentroidEstimate": times,
                                                   "n-Coefficients":
                                                       [i for i in
                                                        range(np.array(CentroidCoefficients["values"]).shape[1])]},
                                           dims=["timeOfDopplerCentroidEstimate", "n-Coefficients"],
                                           attrs=(CentroidCoefficients["attr"] | CentroidCoefficients_doc))
    centroidConfidence_da = xr.DataArray(data=CentroidConfidence["values"],
                                         coords={"timeOfDopplerCentroidEstimate": times},
                                         dims=["timeOfDopplerCentroidEstimate"],
                                         attrs=(CentroidConfidence["attr"] | CentroidConfidence_doc))
    ds["dopplerAmbiguity"] = ambiguity_da
    ds["dopplerAmbiguityConfidence"] = ambiguityConfidence_da
    ds["dopplerCentroidReferenceTime"] = centroidReferenceTime_da
    ds["dopplerCentroidPolynomialPeriod"] = centroidPolynomialPeriod_da
    ds["dopplerCentroidCoefficients"] = centroidCoefficients_da
    ds["dopplerCentroidConfidence"] = centroidConfidence_da
    ds.attrs = (ds_attr | generate_doc_ds(f"{xpath_dict['doppler']['xpath']}/dopplerCentroid", folder_path))
    return ds


def get_dic_doppler_rate_values(dictio):
    content_list = xpath_get(dictio, xpath_dict["doppler"]["xpath"])
    ds_attr = {}
    RateReferenceTime = {
        "values": [],
        "attr": {}
    }
    RateValuesCoefficients = {
        "values": [],
        "attr": {}
    }
    for key in content_list:
        xpath = xpath_dict["doppler"]["xpath"]
        if key == "dopplerRateValues":
            xpath += f"/{key}"
            if isinstance(content_list[key], dict):
                RateReferenceTime["values"].append(float(content_list[key]["dopplerRateReferenceTime"]["#text"]))
                RateReferenceTime["attr"]["RateReferenceTime units"] = content_list[key]["dopplerRateReferenceTime"][
                    "@units"]
                RateReferenceTime["attr"]["dopplerRateReferenceTime_xpath"] = xpath + "/dopplerRateReferenceTime"
                RateValuesCoefficients["values"].append(
                    [float(x) for x in content_list[key]["dopplerRateValuesCoefficients"].split(" ")])
                RateValuesCoefficients["attr"][
                    "dopplerRateValuesCoefficients_xpath"] = xpath + "/dopplerRateValuesCoefficients"
            elif isinstance(content_list[key], list):
                for value in content_list[key]:
                    RateReferenceTime["values"].append(float(content_list[key]["dopplerRateReferenceTime"]["#text"]))
                    RateReferenceTime["attr"]["units"] = content_list[key]["dopplerRateReferenceTime"]["@units"]
                    RateReferenceTime["attr"]["dopplerRateReferenceTime_xpath"] = xpath + "/dopplerRateReferenceTime"
                    RateValuesCoefficients["values"].append(
                        [float(x) for x in content_list[key]["dopplerRateValuesCoefficients"].split(" ")])
                    RateValuesCoefficients["attr"]["dopplerRateValuesCoefficients_xpath"] = \
                        xpath + "/dopplerRateValuesCoefficients"
        elif len(RateReferenceTime["values"]) != 0:
            break
    return {
        "ds_attr": ds_attr,
        "dopplerRateReferenceTime": RateReferenceTime,
        "dopplerRateValuesCoefficients": RateValuesCoefficients
    }


def create_dataset_doppler_rate_values(ds_attr, rateTime, rateCoefficients, folder_path):
    ds = xr.Dataset()
    rateTime_doc = find_doc_in_xsd_files(rateTime["attr"]['dopplerRateReferenceTime_xpath'],
                                         list_xsd_files(rateTime["attr"]['dopplerRateReferenceTime_xpath'],
                                                        folder_path))
    rateCoefficients_doc = find_doc_in_xsd_files(rateCoefficients["attr"]['dopplerRateValuesCoefficients_xpath'],
                                                 list_xsd_files(
                                                     rateCoefficients["attr"]['dopplerRateValuesCoefficients_xpath'],
                                                     folder_path))
    rateCoefficients_da = xr.DataArray(data=np.array(rateCoefficients["values"]),
                                       coords={"dopplerRateReferenceTime": rateTime["values"],
                                               "n-RateValuesCoefficients":
                                                   [i for i in range(np.array(rateCoefficients["values"]).shape[1])]},
                                       dims=["dopplerRateReferenceTime", "n-RateValuesCoefficients"],
                                       attrs=(rateTime["attr"] | rateCoefficients["attr"]
                                              | rateTime_doc | rateCoefficients_doc))
    ds["dopplerRateValues"] = rateCoefficients_da
    ds.attrs = (ds_attr | generate_doc_ds(f"{xpath_dict['doppler']['xpath']}/dopplerRateValues", folder_path))
    return ds


def get_dict_chirp(dictio):
    xpath = xpath_dict["doppler"]["xpath"]
    content_list = xpath_get(dictio, xpath)
    pole = {
        "values": []
    }
    ds_attr = {}
    replicaQualityValid = {
        "values": [],
        "attr": {}
    }
    crossCorrelationWidth = {
        "values": [],
        "attr": {}
    }
    sideLobeLevel = {
        "values": [],
        "attr": {}
    }
    integratedSideLobeRatio = {
        "values": [],
        "attr": {}
    }
    crossCorrelationPeakLoc = {
        "values": [],
        "attr": {}
    }
    chirpPower = {
        "values": [],
        "attr": {}
    }
    amplitudeCoefficients = {
        "values": [],
        "attr": {}
    }
    phaseCoefficients = {
        "values": [],
        "attr": {}
    }
    for key in content_list:
        if key == "chirp":
            xpath += f"/{key}"
            if isinstance(content_list[key], list):
                for value in content_list[key]:
                    pole["values"].append(value["@pole"])
                    for k in value:
                        if isinstance(value[k], str) and ("pole" not in k) and ("@" in k):
                            if value["@pole"] not in list(ds_attr.keys()):
                                ds_attr[value["@pole"]] = {}
                            ds_attr[value["@pole"]][k.replace("@", "")] = value[k]
                        elif (k == "amplitudeCoefficients") or (k == "phaseCoefficients"):
                            eval(k)["values"].append([float(x) for x in value[k].split(" ")])
                            eval(k)["attr"]["xpath"] = xpath + f"/{k}"
                        elif k == 'chirpQuality':
                            prefix_path = f"{xpath}/{k}/"
                            for var in value[k]:
                                eval(var)["attr"]["xpath"] = prefix_path + var
                                if (var == "crossCorrelationPeakLoc") or (var == "crossCorrelationWidth"):
                                    eval(var)["values"].append(float(value[k][var]))
                                elif (var == "sideLobeLevel") or (var == "integratedSideLobeRatio"):
                                    for intern_key in value[k][var]:
                                        if "@" in intern_key:
                                            eval(var)["attr"][intern_key.replace("@", "")] = value[k][var][intern_key]
                                        elif intern_key == "#text":
                                            eval(var)["values"].append(float(value[k][var][intern_key]))
                                elif var == "replicaQualityValid":
                                    eval(var)["values"].append(value[k][var])
                        elif k == "chirpPower":
                            eval(k)["attr"]["xpath"] = xpath + f"/{k}"
                            for intern_key in value[k]:
                                if "@" in intern_key:
                                    eval(k)["attr"][intern_key.replace("@", "")] = value[k][intern_key]
                                elif intern_key == "#text":
                                    eval(k)["values"].append(float(value[k][intern_key]))
            else:
                pole["values"].append(content_list[key]["@pole"])
                for value in content_list[key]:
                    if isinstance(content_list[key][value], str) and ("pole" not in value) and ("@" in value):
                        if content_list[key]["@pole"] not in list(ds_attr.keys()):
                            ds_attr[content_list[key]["@pole"]] = {}
                        ds_attr[content_list[key]["@pole"]][value.replace("@", "")] = content_list[key][value]
                    elif (value == "amplitudeCoefficients") or (value == "phaseCoefficients"):
                        eval(value)["values"].append([float(x) for x in content_list[key][value].split(" ")])
                        eval(value)["attr"]["xpath"] = xpath + f"/{value}"
                    elif value == 'chirpQuality':
                        prefix_path = f"{xpath}/{value}/"
                        for var in content_list[key][value]:
                            eval(var)["attr"]["xpath"] = prefix_path + var
                            if (var == "crossCorrelationPeakLoc") or (var == "crossCorrelationWidth") \
                                    or (var == "replicaQualityValid"):
                                eval(var)["values"].append(parse_value(content_list[key][value][var]))
                            elif (var == "sideLobeLevel") or (var == "integratedSideLobeRatio"):
                                for intern_key in content_list[key][value][var]:
                                    if "@" in intern_key:
                                        eval(var)["attr"][intern_key.replace("@", "")] = \
                                            content_list[key][value][var][intern_key]
                                    elif intern_key == "#text":
                                        eval(var)["values"].append(
                                            parse_value(content_list[key][value][var][intern_key]))
                    elif value == "chirpPower":
                        eval(value)["attr"]["xpath"] = xpath + f"/{value}"
                        for intern_key in content_list[key][value]:
                            if "@" in intern_key:
                                eval(value)["attr"][intern_key.replace("@", "")] = content_list[key][value][intern_key]
                            elif intern_key == "#text":
                                eval(value)["values"].append(parse_value(content_list[key][value][intern_key]))
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
        "phaseCoefficients": phaseCoefficients
    }


def create_dataset_chirp(pole, ds_attr, replicaQualityValid, crossCorrelationWidth, sideLobeLevel,
                         integratedSideLobeRatio, crossCorrelationPeakLoc, chirpPower,
                         amplitudeCoefficients, phaseCoefficients, folder_path):
    ds = xr.Dataset()
    ds.attrs = (ds_attr | generate_doc_ds(f"{xpath_dict['doppler']['xpath']}/chirp", folder_path))
    replicaQualityValid_da = xr.DataArray(data=replicaQualityValid["values"], coords={'pole': pole["values"]},
                                          dims=["pole"],
                                          attrs=(replicaQualityValid["attr"] |
                                                 generate_doc_vars(replicaQualityValid["attr"]["xpath"], folder_path)))
    crossCorrelationWidth_da = xr.DataArray(data=crossCorrelationWidth["values"], coords={'pole': pole["values"]},
                                            dims=["pole"],
                                            attrs=(crossCorrelationWidth["attr"] |
                                                   generate_doc_vars(crossCorrelationWidth["attr"]["xpath"],
                                                                     folder_path)))
    sideLobeLevel_da = xr.DataArray(data=sideLobeLevel["values"], coords={'pole': pole["values"]},
                                    dims=["pole"],
                                    attrs=(sideLobeLevel["attr"] | generate_doc_vars(sideLobeLevel["attr"]["xpath"],
                                                                                     folder_path)))
    integratedSideLobeRatio_da = xr.DataArray(data=integratedSideLobeRatio["values"], coords={'pole': pole["values"]},
                                              dims=["pole"],
                                              attrs=(integratedSideLobeRatio["attr"] |
                                                     generate_doc_vars(integratedSideLobeRatio["attr"]["xpath"],
                                                                       folder_path)))
    crossCorrelationPeakLoc_da = xr.DataArray(data=crossCorrelationPeakLoc["values"], coords={'pole': pole["values"]},
                                              dims=["pole"],
                                              attrs=(crossCorrelationPeakLoc["attr"] |
                                                     generate_doc_vars(crossCorrelationPeakLoc["attr"]["xpath"],
                                                                       folder_path)))
    chirpPower_da = xr.DataArray(data=chirpPower["values"], coords={'pole': pole["values"]},
                                 dims=["pole"],
                                 attrs=(chirpPower["attr"] | generate_doc_vars(chirpPower["attr"]["xpath"],
                                                                               folder_path)))
    amplitudeCoefficients_da = xr.DataArray(data=np.array(amplitudeCoefficients["values"]),
                                            coords={'pole': pole["values"],
                                                    'n-amplitudeCoefficients':
                                                        [i for i in
                                                         range(np.array(amplitudeCoefficients["values"]).shape[1])]},
                                            dims=["pole", "n-amplitudeCoefficients"],
                                            attrs=(amplitudeCoefficients["attr"] |
                                                   generate_doc_vars(amplitudeCoefficients["attr"]["xpath"],
                                                                     folder_path)))
    phaseCoefficients_da = xr.DataArray(data=np.array(phaseCoefficients["values"]),
                                        coords={'pole': pole["values"],
                                                'n-phaseCoefficients':
                                                    [i for i in range(np.array(phaseCoefficients["values"]).shape[1])]},
                                        dims=["pole", "n-phaseCoefficients"],
                                        attrs=(phaseCoefficients["attr"] | generate_doc_vars(
                                            phaseCoefficients["attr"]["xpath"], folder_path)))
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
    xpath = xpath_dict["radarParameters"]["xpath"]
    content_list = xpath_get(dictio, xpath)
    principal_dic = {
        "ds_attr": {}
    }
    ds_attr = list(radar_parameters_key_dict["ds_attributes"])
    vars = list(radar_parameters_key_dict["coord"].keys())
    for var in vars:
        if "incidenceAngleCorrection" in var:
            template_dic = {
                "noiseLevelValues": [],
                "coords": {},
                "attr": {}
            }
        else:
            template_dic = {
                "values": [],
                "coords": {},
                "attr": {}
            }
        for val in radar_parameters_key_dict["coord"][var]:
            template_dic["coords"][val] = []
        principal_dic[var] = template_dic
    for key in content_list:
        if key in ds_attr:
            if (key == "polarizations") or (key == "beams"):
                principal_dic["ds_attr"][key] = content_list[key].split(" ")
            elif isinstance(content_list[key], dict):
                for intern_key in content_list[key]:
                    if "@" in intern_key:
                        principal_dic["ds_attr"][f"{key}_{intern_key.replace('@', '')}"] \
                            = parse_value(content_list[key][intern_key])
                    else:
                        principal_dic["ds_attr"][key] = parse_value(content_list[key][intern_key])
            else:
                principal_dic["ds_attr"][key] = parse_value(content_list[key])
        elif key in vars:
            prefix_path = f"{xpath}/{key}"
            if isinstance(content_list[key], list):
                principal_dic[key]["attr"]["xpath"] = prefix_path
                for value in content_list[key]:
                    for intern_key in value:
                        if intern_key.replace("@", "") in radar_parameters_key_dict["coord"][key]:
                            principal_dic[key]["coords"][intern_key.replace("@", "")] \
                                .append(parse_value(value[intern_key]))
                        elif intern_key == "#text":
                            principal_dic[key]["values"].append(parse_value(value[intern_key]))
                        else:
                            principal_dic[key]["attr"][intern_key.replace("@", "")] = parse_value(value[intern_key])
        elif isinstance(content_list[key], list):
            prefix_path = f"{xpath}/{key}/"
            for value in content_list[key]:
                var_name = ""
                # referenceNoiseLevel case
                for k in value:
                    if "@" in k:
                        var_name = f"{k.replace('@', '')}_{value[k].replace(' ', '_')}"
                    elif isinstance(value[k], str):
                        principal_dic[var_name]["attr"][k] = parse_value(value[k])
                    elif isinstance(value[k], dict):
                        principal_dic[var_name]["attr"]["xpath"] = prefix_path + k
                        for intern_key in value[k]:
                            if "@" in intern_key:
                                principal_dic[var_name]["attr"][f"{k}_{intern_key.replace('@', '')}"] = value[k][
                                    intern_key]
                            else:
                                principal_dic[var_name]["noiseLevelValues"] \
                                    = [parse_value(x) for x in value[k][intern_key].split(" ")]
                                principal_dic[var_name]["coords"]["NbOfNoiseLevelValues"] = \
                                    np.arange(len(value[k][intern_key].split(" ")))

    return principal_dic


def create_dataset_radar_parameters(dictio, folder_path):
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
                data = dictio[key]['noiseLevelValues']
                coords = dictio[key]["coords"]
                if "Beta" in key:
                    general_ds['noiseLevelValues_BetaNought'] = xr.DataArray(data=data, dims=dims,
                                                                             attrs=(attr |
                                                                                    generate_doc_vars(attr["xpath"],
                                                                                                      folder_path)))
                elif "Sigma" in key:
                    general_ds['noiseLevelValues_SigmaNought'] = xr.DataArray(data=data, dims=dims,
                                                                              attrs=(attr |
                                                                                     generate_doc_vars(attr["xpath"],
                                                                                                       folder_path)))
                elif "Gamma" in key:
                    general_ds['noiseLevelValues_Gamma'] = xr.DataArray(data=data, dims=dims,
                                                                        attrs=(attr | generate_doc_vars(attr["xpath"],
                                                                                                        folder_path)))
            else:
                if len(dims) == 2:
                    data = create_2d_matrix(dictio[key]["coords"][dims[0]], dictio[key]["coords"][dims[1]],
                                            dictio[key]["values"])
                    coords[dims[0]] = np.unique(dictio[key]["coords"][dims[0]])
                    coords[dims[1]] = np.unique(dictio[key]["coords"][dims[1]])
                elif len(dims) == 1:
                    data = dictio[key]["values"]
                    coords[dims[0]] = dictio[key]["coords"][dims[0]]
                general_ds[key] = xr.DataArray(data=data, dims=dims, coords=coords,
                                               attrs=(attr | generate_doc_vars(attr["xpath"], folder_path)))
    general_ds.attrs |= generate_doc_ds(xpath_dict["radarParameters"]["xpath"], folder_path)
    return general_ds


def create_2d_matrix(lines, cols, vals):
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
    xpath = xpath_dict["geolocation_grid"]["xpath"].split("/geographicInformation")[0]
    content_list = xpath_get(dictio, xpath)
    attr = {}
    for key in content_list:
        if isinstance(content_list[key], str):
            attr[key] = parse_value(content_list[key])
        elif key == "rasterAttributes":
            for value in content_list[key]:
                if isinstance(content_list[key][value], str):
                    attr[f"{key}_{value}"] = parse_value(content_list[key][value])
                elif isinstance(content_list[key][value], dict):
                    dico_keys = list(content_list[key][value].keys())
                    for k in dico_keys:
                        attr[f"{key}_{value}_{k.replace('@', '').replace('#text', 'value')}"] = \
                            parse_value(content_list[key][value][k])
    return attr


def close_string_dic(string_dic):
    tmp_list = list(string_dic)
    tmp_list[-1] = '}'
    return "".join(tmp_list)


def parse_value(value):
    import ast
    try:
        return ast.literal_eval(value)
    except:
        return value


def list_xsd_files(xpath, folder_path):
    xsd_folder = folder_path + "/schemas"
    var_name = xpath.split('/')[-1]
    parent_path = os.path.dirname(xpath)
    parent_var_name = parent_path.split('/')[-1]
    grand_parent_path = os.path.dirname(parent_path)
    grand_parent_var_name = grand_parent_path.split('/')[-1]
    list_files = os.listdir(xsd_folder)
    interesting_files = [[f"{xsd_folder}/{file}" for file in list_files if (var_name in file)],
                         [f"{xsd_folder}/{file}" for file in list_files if (parent_var_name in file)],
                         [f"{xsd_folder}/{file}" for file in list_files if (grand_parent_var_name in file)]]
    return interesting_files


def find_doc_in_xsd_files(xpath, interesting_files):
    var_name = xpath.split('/')[-1]
    xsd_path = "xsd:schema/xsd:complexType/xsd:sequence/xsd:element"
    description = []
    for element in interesting_files:
        # index = interesting_files.index(element)
        if len(element) != 0:
            for pathname in element:
                with open(pathname, 'rb') as f:
                    xml_content = f.read()
                    dic = xmltodict.parse(xml_content)
                    f.close()
                content_list = xpath_get(dic, xsd_path)
                for value in content_list:
                    if value["@name"] == var_name:
                        description.append(value['xsd:annotation']['xsd:documentation']
                                           .replace("  ", "")
                                           .replace("\n", "")
                                           .replace("\t", " "))
    return {f"Description": description[0]}


def get_type_for_pole(folder_path):
    pathname = os.path.join(folder_path, "schemas", "rs2prod_chirp.xsd")
    xpath = "/xsd:schema/xsd:complexType/xsd:attribute"
    with open(pathname, 'rb') as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()
    content_list = xpath_get(dic, xpath)
    for values in content_list:
        if values["@name"] == "pole":
            return {"type": values["@type"]}


def find_doc_for_ds_in_xsd_files(xpath, interesting_files):
    if "geolocationGrid" in xpath:
        ds_name = xpath.split('/')[-2]
        interesting_files = interesting_files[1]
    else:
        ds_name = xpath.split('/')[-1]
        interesting_files = interesting_files[0]
    xsd_path = "xsd:schema/xsd:annotation/xsd:documentation"
    for element in interesting_files:
        # index = interesting_files.index(element)
        if len(element) != 0:
            with open(element, 'rb') as f:
                xml_content = f.read()
                dic = xmltodict.parse(xml_content)
                f.close()
            content_list = xpath_get(dic, xsd_path) \
                .replace("  ", "") \
                .replace("\n", "") \
                .replace("\t", " ")
    return {"Description": content_list}


def generate_doc_vars(xpath, folder_path):
    return find_doc_in_xsd_files(xpath, list_xsd_files(xpath, folder_path))


def generate_doc_ds(xpath, folder_path):
    return find_doc_for_ds_in_xsd_files(xpath, list_xsd_files(xpath, folder_path))


def list_lut_files(folder_path):
    list_files = os.listdir(folder_path)
    interesting_files = [os.path.join(folder_path, file) for file in list_files if ("lut" in file)]
    return interesting_files


def create_dataset_lut(dictio):
    final_lut_dict = {
        "attrs": {}
    }
    ds = xr.Dataset()
    for value in dictio["lut"]:
        if "@" not in value:
            if value != "gains":
                final_lut_dict["attrs"][value] = parse_value(dictio["lut"][value])
            else:
                final_lut_dict[value] = [parse_value(x) for x in dictio["lut"][value].split(" ")]
    da = xr.DataArray(data=final_lut_dict["gains"])
    ds["lut"] = da
    ds.attrs = final_lut_dict["attrs"]
    return ds


def lut_processing(files):
    lut_dict = {}
    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        with open(file, 'rb') as f:
            xml_content = f.read()
            dic = xmltodict.parse(xml_content)
            f.close()
        lut_dict[filename] = create_dataset_lut(dic)
    return lut_dict


def xml_parser(folder_path):
    # get product.xml path
    product_xml_path = os.path.join(folder_path, "product.xml")

    # get product.xml content as a dict
    with open(product_xml_path, 'rb') as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()

    # Create dictionaries then dataset, and fill them in a datatree
    dic_geo = create_dic_geolocation_grid(dic)
    ds_geo = create_dataset_geolocation_grid(dic_geo, folder_path)
    ds_geo.attrs = fill_image_attribute(dic)
    dic_orbit_information = get_dic_orbit_information(dic)
    ds_orbit_info = create_dataset_orbit_information(dic_orbit_information["ds_attr"],
                                                     dic_orbit_information["timestamp"],
                                                     dic_orbit_information["xPosition"],
                                                     dic_orbit_information["yPosition"],
                                                     dic_orbit_information["zPosition"],
                                                     dic_orbit_information["xVelocity"],
                                                     dic_orbit_information["yVelocity"],
                                                     dic_orbit_information["zVelocity"],
                                                     folder_path)
    dic_attitude_info = get_dic_attitude_info(dic)
    ds_attitude_info = create_dataset_attitude_information(dic_attitude_info["ds_attr"],
                                                           dic_attitude_info["timestamp"],
                                                           dic_attitude_info["yaw"],
                                                           dic_attitude_info["roll"],
                                                           dic_attitude_info["pitch"],
                                                           folder_path)
    ds_orbit_attitude_info = xr.merge([ds_attitude_info, ds_orbit_info])
    ds_orbit_attitude_info.attrs["Description"] += ". " + ds_orbit_info.attrs["Description"] + "."
    dt = datatree.DataTree()
    """dt["orbitAndAttitude"] = datatree.DataTree.from_dict(
        {"orbitInformation": ds_orbit_info, "attitudeInformation": ds_attitude_info})"""
    dt["orbitAndAttitude"] = datatree.DataTree(data=ds_orbit_attitude_info)
    dt["geolocationGrid"] = datatree.DataTree(data=ds_geo)
    dic_doppler_centroid = get_dict_doppler_centroid(dic)
    ds_doppler_centroid = create_dataset_doppler_centroid(dic_doppler_centroid["ds_attr"],
                                                          dic_doppler_centroid["timeOfDopplerCentroidEstimate"],
                                                          dic_doppler_centroid["dopplerAmbiguity"],
                                                          dic_doppler_centroid["dopplerAmbiguityConfidence"],
                                                          dic_doppler_centroid["dopplerCentroidReferenceTime"],
                                                          dic_doppler_centroid["dopplerCentroidPolynomialPeriod"],
                                                          dic_doppler_centroid["dopplerCentroidCoefficients"],
                                                          dic_doppler_centroid["dopplerCentroidConfidence"],
                                                          folder_path
                                                          )
    dt["imageGenerationParameters/doppler/dopplerCentroid"] = datatree.DataTree(data=ds_doppler_centroid)
    dic_doppler_rate_values = get_dic_doppler_rate_values(dic)
    ds_doppler_rate_values = create_dataset_doppler_rate_values(dic_doppler_rate_values["ds_attr"],
                                                                dic_doppler_rate_values["dopplerRateReferenceTime"],
                                                                dic_doppler_rate_values[
                                                                    "dopplerRateValuesCoefficients"], folder_path)
    dt["imageGenerationParameters/doppler/dopplerRateValues"] = datatree.DataTree(data=ds_doppler_rate_values)
    dic_chirp = get_dict_chirp(dic)
    ds_chirp = create_dataset_chirp(dic_chirp["pole"], dic_chirp["ds_attr"], dic_chirp["replicaQualityValid"],
                                    dic_chirp["crossCorrelationWidth"], dic_chirp["sideLobeLevel"],
                                    dic_chirp["integratedSideLobeRatio"], dic_chirp["crossCorrelationPeakLoc"],
                                    dic_chirp["chirpPower"], dic_chirp["amplitudeCoefficients"],
                                    dic_chirp["phaseCoefficients"], folder_path)
    dt["imageGenerationParameters/chirp"] = datatree.DataTree(data=ds_chirp)
    radar_parameters_dic = get_dict_radar_parameters(dic)
    ds_radar_parameters = create_dataset_radar_parameters(radar_parameters_dic, folder_path)
    dt["radarParameters"] = datatree.DataTree(data=ds_radar_parameters)
    lut_dict = lut_processing(list_lut_files(folder_path))
    dt["lut"] = datatree.DataTree.from_dict(lut_dict)
    return dt


if __name__ == '__main__':
    xml_parser("/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV/"
               "2010/288/RS2_OK72200_PK649463_DK111111_SCWA_20101015_210132_VV_SGF")
    # xml_parser("/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV_VH/2022/084/RS2_OK134765_PK1183680_DK1147849_SCWA_20220325_133717_VV_VH_SGF")

# TODO : create doc to fill documentation automatically ( see example on github --> Antoine messages)
# TODO : read tif images
# TODO : put xpath for coord???
# TODO : homogenize the func for each datasets
