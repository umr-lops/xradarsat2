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

path = "/home/datawork-cersat-public/cache/public/ftp/project/sarwing/SAFE_REF" \
       "/RS2_OK91548_PK811670_DK739881_SCWA_20170908_105352_VV_VH_SGF/product.xml"

xpath_dict = {
    "geolocation_grid":
        {
            "xpath": "/product/imageAttributes/geographicInformation/geolocationGrid/imageTiePoint",
            "xsdpath": "/home/datawork-cersat-public/cache/public/ftp/project/sarwing/SAFE_REF"
                       "/RS2_OK91548_PK811670_DK739881_SCWA_20170908_105352_VV_VH_SGF/schemas/rs2prod_geolocationGrid"
                       ".xsd",
            "info_xsd_path": "/xsd:schema/xsd:complexType/xsd:sequence/xsd:element/xsd:annotation/xsd:documentation"
                             "/text()",
        },
    "orbit_information":
        {
            "xpath": "/product/sourceAttributes/orbitAndAttitude/orbitInformation",
            "xsdpath": "/home/datawork-cersat-public/cache/public/ftp/project/sarwing/SAFE_REF"
                       "/RS2_OK91548_PK811670_DK739881_SCWA_20170908_105352_VV_VH_SGF/schemas/rs2prod_stateVector.xsd "
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


def get_lists_geolocation_grid(dictio):
    content_list = xpath_get(dictio, xpath_dict["geolocation_grid"]["xpath"])
    lines = []
    pixs = []
    los = []
    las = []
    hes = []
    units = ["", "", ""]
    for element in content_list:
        lines.append(element['imageCoordinate']['line'])
        pixs.append(element['imageCoordinate']['pixel'])
        los.append(element['geodeticCoordinate']['longitude']['#text'])
        las.append(element['geodeticCoordinate']['latitude']['#text'])
        hes.append(element['geodeticCoordinate']['height']['#text'])
        units[0] = element['geodeticCoordinate']['longitude']['@units']
        units[1] = element['geodeticCoordinate']['latitude']['@units']
        units[2] = element['geodeticCoordinate']['height']['@units']
    return lines, pixs, los, las, hes, units


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


def create_dataset_orbit_information(ds_attr, timestamp, xPos, yPos, zPos, xVel, yVel, zVel):
    ds = xr.Dataset()
    xpos_da = xr.DataArray(data=xPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(xPos["attr"] | generate_doc_vars(xPos["attr"]["xpath"])))
    ypos_da = xr.DataArray(data=yPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(yPos["attr"] | generate_doc_vars(yPos["attr"]["xpath"])))
    zpos_da = xr.DataArray(data=zPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(zPos["attr"] | generate_doc_vars(zPos["attr"]["xpath"])))
    xvel_da = xr.DataArray(data=xVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(xVel["attr"] | generate_doc_vars(xVel["attr"]["xpath"])))
    yvel_da = xr.DataArray(data=yVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(yVel["attr"] | generate_doc_vars(yVel["attr"]["xpath"])))
    zvel_da = xr.DataArray(data=zVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp",
                           attrs=(zVel["attr"] | generate_doc_vars(zVel["attr"]["xpath"])))
    ds["xPosition"] = xpos_da
    ds["yPosition"] = ypos_da
    ds["zPosition"] = zpos_da
    ds["xVelocity"] = xvel_da
    ds["yVelocity"] = yvel_da
    ds["zVelocity"] = zvel_da
    ds.attrs = ds_attr
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


def create_dataset_attitude_information(ds_attr, timestamp, yaw, roll, pitch):
    ds = xr.Dataset()
    interesting_files_yaw = list_xsd_files(yaw["attr"]["xpath"])
    yaw_doc = find_doc_in_xsd_files(yaw["attr"]["xpath"], interesting_files_yaw)
    interesting_files_roll = list_xsd_files(roll["attr"]["xpath"])
    roll_doc = find_doc_in_xsd_files(roll["attr"]["xpath"], interesting_files_roll)
    interesting_files_pitch = list_xsd_files(pitch["attr"]["xpath"])
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
    ds.attrs = ds_attr
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
                                    CentroidPolynomialPeriod, CentroidCoefficients, CentroidConfidence):
    ds = xr.Dataset()
    Ambiguity_doc = find_doc_in_xsd_files(Ambiguity["attr"]["xpath"], list_xsd_files(Ambiguity["attr"]["xpath"]))
    AmbiguityConfidence_doc = find_doc_in_xsd_files(AmbiguityConfidence["attr"]["xpath"],
                                                    list_xsd_files(AmbiguityConfidence["attr"]["xpath"]))
    CentroidReferenceTime_doc = find_doc_in_xsd_files(CentroidReferenceTime["attr"]["xpath"],
                                                      list_xsd_files(CentroidReferenceTime["attr"]["xpath"]))
    CentroidPolynomialPeriod_doc = find_doc_in_xsd_files(CentroidPolynomialPeriod["attr"]["xpath"],
                                                         list_xsd_files(CentroidPolynomialPeriod["attr"]["xpath"]))
    CentroidCoefficients_doc = find_doc_in_xsd_files(CentroidCoefficients["attr"]["xpath"],
                                                     list_xsd_files(CentroidCoefficients["attr"]["xpath"]))
    CentroidConfidence_doc = find_doc_in_xsd_files(CentroidConfidence["attr"]["xpath"],
                                                   list_xsd_files(CentroidConfidence["attr"]["xpath"]))
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
    ds.attrs = ds_attr
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


def create_dataset_doppler_rate_values(ds_attr, rateTime, rateCoefficients):
    ds = xr.Dataset()
    rateTime_doc = find_doc_in_xsd_files(rateTime["attr"]['dopplerRateReferenceTime_xpath'],
                                         list_xsd_files(rateTime["attr"]['dopplerRateReferenceTime_xpath']))
    rateCoefficients_doc = find_doc_in_xsd_files(rateCoefficients["attr"]['dopplerRateValuesCoefficients_xpath'],
                                                 list_xsd_files(
                                                     rateCoefficients["attr"]['dopplerRateValuesCoefficients_xpath']))
    rateCoefficients_da = xr.DataArray(data=np.array(rateCoefficients["values"]),
                                       coords={"dopplerRateReferenceTime": rateTime["values"],
                                               "n-RateValuesCoefficients":
                                                   [i for i in range(np.array(rateCoefficients["values"]).shape[1])]},
                                       dims=["dopplerRateReferenceTime", "n-RateValuesCoefficients"],
                                       attrs=(rateTime["attr"] | rateCoefficients["attr"]
                                              | rateTime_doc | rateCoefficients_doc))
    ds["dopplerRateValues"] = rateCoefficients_da
    ds.attrs = ds_attr
    return ds


def get_dict_chirp(dictio):
    xpath = xpath_dict["doppler"]["xpath"]
    content_list = xpath_get(dictio, xpath)
    pole = {
        "values": []
    }
    ds_attr = {
        "VV": {},
        "VH": {}
    }
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
            for value in content_list[key]:
                pole["values"].append(value["@pole"])
                for k in value:
                    if isinstance(value[k], str) and ("pole" not in k) and ("@" in k):
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
                         amplitudeCoefficients, phaseCoefficients):
    ds = xr.Dataset()
    ds.attrs = ds_attr
    replicaQualityValid_da = xr.DataArray(data=replicaQualityValid["values"], coords={'pole': pole["values"]},
                                          dims=["pole"],
                                          attrs=(replicaQualityValid["attr"] |
                                                 generate_doc_vars(replicaQualityValid["attr"]["xpath"])))
    crossCorrelationWidth_da = xr.DataArray(data=crossCorrelationWidth["values"], coords={'pole': pole["values"]},
                                            dims=["pole"],
                                            attrs=(crossCorrelationWidth["attr"] |
                                                   generate_doc_vars(crossCorrelationWidth["attr"]["xpath"])))
    sideLobeLevel_da = xr.DataArray(data=sideLobeLevel["values"], coords={'pole': pole["values"]},
                                    dims=["pole"],
                                    attrs=(sideLobeLevel["attr"] | generate_doc_vars(sideLobeLevel["attr"]["xpath"])))
    integratedSideLobeRatio_da = xr.DataArray(data=integratedSideLobeRatio["values"], coords={'pole': pole["values"]},
                                              dims=["pole"],
                                              attrs=(integratedSideLobeRatio["attr"] |
                                                     generate_doc_vars(integratedSideLobeRatio["attr"]["xpath"])))
    crossCorrelationPeakLoc_da = xr.DataArray(data=crossCorrelationPeakLoc["values"], coords={'pole': pole["values"]},
                                              dims=["pole"],
                                              attrs=(crossCorrelationPeakLoc["attr"] |
                                                     generate_doc_vars(crossCorrelationPeakLoc["attr"]["xpath"])))
    chirpPower_da = xr.DataArray(data=chirpPower["values"], coords={'pole': pole["values"]},
                                 dims=["pole"],
                                 attrs=(chirpPower["attr"] | generate_doc_vars(chirpPower["attr"]["xpath"])))
    amplitudeCoefficients_da = xr.DataArray(data=np.array(amplitudeCoefficients["values"]),
                                            coords={'pole': pole["values"],
                                                    'n-amplitudeCoefficients':
                                                        [i for i in
                                                         range(np.array(amplitudeCoefficients["values"]).shape[1])]},
                                            dims=["pole", "n-amplitudeCoefficients"],
                                            attrs=(amplitudeCoefficients["attr"] |
                                                   generate_doc_vars(amplitudeCoefficients["attr"]["xpath"])))
    phaseCoefficients_da = xr.DataArray(data=np.array(phaseCoefficients["values"]),
                                        coords={'pole': pole["values"],
                                                'n-phaseCoefficients':
                                                    [i for i in range(np.array(phaseCoefficients["values"]).shape[1])]},
                                        dims=["pole", "n-phaseCoefficients"],
                                        attrs=(phaseCoefficients["attr"] | generate_doc_vars(phaseCoefficients["attr"]["xpath"])))
    ds["replicaQualityValid"] = replicaQualityValid_da
    ds["crossCorrelationWidth"] = crossCorrelationWidth_da
    ds["sideLobeLevel"] = sideLobeLevel_da
    ds["integratedSideLobeRatio"] = integratedSideLobeRatio_da
    ds["crossCorrelationPeakLoc"] = crossCorrelationPeakLoc_da
    ds["chirpPower"] = chirpPower_da
    ds["amplitudeCoefficients"] = amplitudeCoefficients_da
    ds["phaseCoefficients"] = phaseCoefficients_da
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


def create_dataset_radar_parameters(dictio):
    general_ds = xr.Dataset()
    BetaNought_ds = xr.Dataset()
    SigmaNought_ds = xr.Dataset()
    Gamma_ds = xr.Dataset()
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
                    BetaNought_ds['noiseLevelValues'] = xr.DataArray(data=data, dims=dims,
                                                                     attrs=(attr | generate_doc_vars(attr["xpath"])))
                    BetaNought_ds.attrs = dictio["ds_attr"]
                elif "Sigma" in key:
                    SigmaNought_ds['noiseLevelValues'] = xr.DataArray(data=data, dims=dims,
                                                                      attrs=(attr | generate_doc_vars(attr["xpath"])))
                    SigmaNought_ds.attrs = dictio["ds_attr"]
                elif "Gamma" in key:
                    Gamma_ds['noiseLevelValues'] = xr.DataArray(data=data, dims=dims,
                                                                attrs=(attr | generate_doc_vars(attr["xpath"])))
                    Gamma_ds.attrs = dictio["ds_attr"]
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
                                               attrs=(attr | generate_doc_vars(attr["xpath"])))
    return general_ds, BetaNought_ds, SigmaNought_ds, Gamma_ds


"""def get_dic_chirp(dictio):
    xpath = xpath_dict["doppler"]["xpath"]
    content_list = xpath_get(dictio, xpath)
    principal_dic = {
        "ds_attrs": {
        "VV": {},
        "VH": {}
        },
        "pole": {
            "values": []
        }
    }
    for key in content_list:
        if key == "chirp":
            for value in content_list[key]:
                principal_dic["pole"]["values"].append(value["@pole"])
                tmp_dic_attributes = {}
                for k in value:
                    if isinstance(value[k], str) and ("pole" not in k) and ("@" in k):
                        tmp_dic_attributes[k.replace("@", "")] = value[k]
        ############################################################################
                    elif isinstance(value[k], str) and ("pole" not in k) and ("@" not in k):
                        exec("%s = %s" % ("dic_name", "{}"))
                        eval("dic_name")["values"] = []
        ############################################################################
                    elif k == "chirpPower":
                        principal_dic_keys = list(principal_dic.keys())
                        if k not in principal_dic_keys:
                            dic_string = "{"
                            for intern_key in value[k]:
                                if intern_key == "#text":
                                    dic_string += f"'{intern_key.replace('#text', 'values')}':[{value[k][intern_key]}],"
                                else:
                                    dic_string += f"'{intern_key.replace('@', '')}':'{value[k][intern_key]}',"
                            dic_string = close_string_dic(dic_string)
                            exec("%s = %s" % ("dic_name", dic_string))
                            principal_dic[k] = eval("dic_name")
                        else:
                            for intern_key in value[k]:
                                if isinstance(
                                        principal_dic[k][intern_key.replace("@", "").replace("#text", "values")],
                                        list):
                                    principal_dic[k][intern_key.replace("@", "").replace("#text", "values")].append(
                                        value[k][intern_key])
                                else:
                                    principal_dic[k][intern_key.replace("@", "").replace("#text", "values")] = \
                                    value[k][intern_key]
                    elif k == "chirpQuality":
                        for var in value[k]:
                            principal_dic_keys = list(principal_dic.keys())
                            if isinstance(value[k][var], str):
                                if var not in principal_dic_keys:
                                    exec("%s = %s" % ("dic_name", "{'values':[]}"))
                                    eval("dic_name")['values'].append(value[k][var])
                                    principal_dic[var] = eval("dic_name")
                                else:
                                    principal_dic[var]["values"].append(value[k][var])
                            elif isinstance(value[k][var], dict):
                                if var not in principal_dic_keys:
                                    dict_string = "{"
                                    for intern_key in value[k][var]:
                                        if intern_key == '#text':
                                            dict_string += f"'{intern_key.replace('#text', 'values')}':[{value[k][var][intern_key]}],"
                                        else:
                                            dict_string += f"'{intern_key.replace('@', '')}':'{value[k][var][intern_key]}',"
                                    dict_string = close_string_dic(dict_string)
                                    exec("%s = %s" % ("dic_name", dict_string))
                                    principal_dic[var] = eval("dic_name")
                                else:
                                    for intern_key in value[k][var]:
                                        if isinstance(principal_dic[var][intern_key.replace("@", "").replace("#text", "values")], list):
                                            principal_dic[var][intern_key.replace("@", "").replace("#text", "values")].append(value[k][var][intern_key])
                                        else:
                                            principal_dic[var][intern_key.replace("@", "").replace("#text", "values")] = value[k][var][intern_key]
                principal_dic["ds_attrs"][value["@pole"]] = tmp_dic_attributes
    return principal_dic"""


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


def create_data_array_geolocation_grid(data, name, coord_line, coord_pix, unit):
    if name == "longitude":
        xpath_suffix = '/geodeticCoordinate/longitude'
    elif name == "latitude":
        xpath_suffix = '/geodeticCoordinate/latitude'
    elif name == "height":
        xpath_suffix = '/geodeticCoordinate/height'
    xpath = f"{xpath_dict['geolocation_grid']['xpath']}{xpath_suffix}"
    return xr.DataArray(data=data, name=name,
                        coords={"line": np.unique(np.array(coord_line)), "pixel": np.unique(np.array(coord_pix))},
                        dims=['line', "pixel"], attrs=({"units": unit, "xpath": xpath} | generate_doc_vars(xpath)))


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


def isfloat(x):
    try:
        a = float(x)
    except (TypeError, ValueError):
        return False
    else:
        return True


def isint(x):
    try:
        a = float(x)
        b = int(a)
    except (TypeError, ValueError):
        return False
    else:
        return a == b


def parse_value(value):
    import ast
    try:
        return ast.literal_eval(value)
    except:
        return value


def list_xsd_files(xpath):
    root = "/home/datawork-cersat-public/cache/public/ftp/project/sarwing/SAFE_REF" \
           "/RS2_OK91548_PK811670_DK739881_SCWA_20170908_105352_VV_VH_SGF/schemas"
    var_name = xpath.split('/')[-1]
    parent_path = os.path.dirname(xpath)
    parent_var_name = parent_path.split('/')[-1]
    grand_parent_path = os.path.dirname(parent_path)
    grand_parent_var_name = grand_parent_path.split('/')[-1]
    list_files = os.listdir(root)
    interesting_files = [[f"{root}/{file}" for file in list_files if (var_name in file)],
                         [f"{root}/{file}" for file in list_files if (parent_var_name in file)],
                         [f"{root}/{file}" for file in list_files if (grand_parent_var_name in file)]]
    return interesting_files


def find_doc_in_xsd_files(xpath, interesting_files):
    var_name = xpath.split('/')[-1]
    xsd_path = "xsd:schema/xsd:complexType/xsd:sequence/xsd:element"
    description = []
    for element in interesting_files:
        index = interesting_files.index(element)
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
    return {f"Description_{var_name}": description[0]}


def generate_doc_vars(xpath):
    return find_doc_in_xsd_files(xpath, list_xsd_files(xpath))


def xml_parser(pathname):
    lines = []
    pixs = []
    los = []
    las = []
    hes = []
    units = ["", "", ""]
    with open(pathname, 'rb') as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()
    lines, pixs, los, las, hes, units = get_lists_geolocation_grid(dic)
    lines = [int(float(lines[k])) for k in range(len(lines))]
    pixs = [int(float(pixs[k])) for k in range(len(pixs))]
    da_los = create_data_array_geolocation_grid(create_2d_matrix(lines, pixs, los), "longitude",
                                                lines, pixs, units[0])
    da_las = create_data_array_geolocation_grid(create_2d_matrix(lines, pixs, las), "latitude",
                                                lines, pixs, units[1])
    da_hes = create_data_array_geolocation_grid(create_2d_matrix(lines, pixs, hes), "height",
                                                lines, pixs, units[2])
    with open(xpath_dict["geolocation_grid"]["xsdpath"], 'rb') as f:
        geo_xsd_content = f.read()
        geo_xsd_dic = xmltodict.parse(geo_xsd_content)
        f.close()
    ds_geo = xr.Dataset()
    ds_geo['latitude'] = da_las
    ds_geo['longitude'] = da_los
    ds_geo['height'] = da_hes
    ds_geo.attrs = {"Description": xpath_get(geo_xsd_dic, xpath_dict["geolocation_grid"]["info_xsd_path"])}
    dic_orbit_information = get_dic_orbit_information(dic)
    ds_orbit_info = create_dataset_orbit_information(dic_orbit_information["ds_attr"],
                                                     dic_orbit_information["timestamp"],
                                                     dic_orbit_information["xPosition"],
                                                     dic_orbit_information["yPosition"],
                                                     dic_orbit_information["zPosition"],
                                                     dic_orbit_information["xVelocity"],
                                                     dic_orbit_information["yVelocity"],
                                                     dic_orbit_information["zVelocity"])
    dic_attitude_info = get_dic_attitude_info(dic)
    ds_attitude_info = create_dataset_attitude_information(dic_attitude_info["ds_attr"],
                                                           dic_attitude_info["timestamp"],
                                                           dic_attitude_info["yaw"],
                                                           dic_attitude_info["roll"],
                                                           dic_attitude_info["pitch"])
    dt = datatree.DataTree()
    dt["orbitAndAttitude"] = datatree.DataTree.from_dict(
        {"orbitInformation": ds_orbit_info, "attitudeInformation": ds_attitude_info})
    dt["imageAttributes/geographicInformation/geolocationGrid"] = datatree.DataTree(data=ds_geo)
    dic_doppler_centroid = get_dict_doppler_centroid(dic)
    ds_doppler_centroid = create_dataset_doppler_centroid(dic_doppler_centroid["ds_attr"],
                                                          dic_doppler_centroid["timeOfDopplerCentroidEstimate"],
                                                          dic_doppler_centroid["dopplerAmbiguity"],
                                                          dic_doppler_centroid["dopplerAmbiguityConfidence"],
                                                          dic_doppler_centroid["dopplerCentroidReferenceTime"],
                                                          dic_doppler_centroid["dopplerCentroidPolynomialPeriod"],
                                                          dic_doppler_centroid["dopplerCentroidCoefficients"],
                                                          dic_doppler_centroid["dopplerCentroidConfidence"]
                                                          )
    dt["imageGenerationParameters/doppler/dopplerCentroid"] = datatree.DataTree(data=ds_doppler_centroid)
    dt["imageAttributes"].attrs = fill_image_attribute(dic)
    dic_doppler_rate_values = get_dic_doppler_rate_values(dic)
    ds_doppler_rate_values = create_dataset_doppler_rate_values(dic_doppler_rate_values["ds_attr"],
                                                                dic_doppler_rate_values["dopplerRateReferenceTime"],
                                                                dic_doppler_rate_values[
                                                                    "dopplerRateValuesCoefficients"])
    dt["imageGenerationParameters/doppler/dopplerRateValues"] = datatree.DataTree(data=ds_doppler_rate_values)
    dic_chirp = get_dict_chirp(dic)
    ds_chirp = create_dataset_chirp(dic_chirp["pole"], dic_chirp["ds_attr"], dic_chirp["replicaQualityValid"],
                                    dic_chirp["crossCorrelationWidth"], dic_chirp["sideLobeLevel"],
                                    dic_chirp["integratedSideLobeRatio"], dic_chirp["crossCorrelationPeakLoc"],
                                    dic_chirp["chirpPower"], dic_chirp["amplitudeCoefficients"],
                                    dic_chirp["phaseCoefficients"])
    dt["imageGenerationParameters/chirp"] = datatree.DataTree(data=ds_chirp)
    radar_parameters_dic = get_dict_radar_parameters(dic)
    ds_radar_parameters, Beta_ds, Sigma_ds, Gamma_ds = create_dataset_radar_parameters(radar_parameters_dic)
    dt["radarParameters"] = datatree.DataTree(data=ds_radar_parameters)
    dt["radarParameters/referenceNoiseLevel/incidenceAngleCorrection_Sigma_Nought"] = datatree.DataTree(data=Sigma_ds)
    dt["radarParameters/referenceNoiseLevel/incidenceAngleCorrection_Beta_Nought"] = datatree.DataTree(data=Beta_ds)
    dt["radarParameters/referenceNoiseLevel/incidenceAngleCorrection_Gamma"] = datatree.DataTree(data=Gamma_ds)
    return dt


if __name__ == '__main__':
    xml_parser(path)

"""# TODO : create doc to fill documentation automatically ( see example on github --> Antoine messages)
# TODO: fill  datasets with xsd info
# TODO : read tif images"""
# TODO : put xpath for coord???
# create a func to list xsd files which contain a var name, and if doesn't found,
# search for files which contains the parent file (+1 in xpath level)
# create func which find the xpath of the doc in the respecting file
# Then, fill dataArrays attributes with this xsd description
#TODO : automate the product.xml path and the root of xsd path (files with SAFE??)