# imports
import numpy as np
from lxml import etree
from lxml import objectify
import xmltodict
import xarray as xr
import time
import xsar
import datatree


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
    xpos_da = xr.DataArray(data=xPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=xPos["attr"])
    ypos_da = xr.DataArray(data=yPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=yPos["attr"])
    zpos_da = xr.DataArray(data=zPos["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=zPos["attr"])
    xvel_da = xr.DataArray(data=xVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=xVel["attr"])
    yvel_da = xr.DataArray(data=yVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=yVel["attr"])
    zvel_da = xr.DataArray(data=zVel["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=zVel["attr"])
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
    yaw_da = xr.DataArray(data=yaw["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=yaw["attr"])
    roll_da = xr.DataArray(data=roll["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=roll["attr"])
    pitch_da = xr.DataArray(data=pitch["values"], coords={"timeStamp": timestamp}, dims="timeStamp", attrs=pitch["attr"])
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
        "values": []
    }
    AmbiguityConfidence = {
        "values": []
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
        "values": []
    }
    CentroidConfidence = {
        "values": []
    }
    for key in content_list:
        if key == "dopplerCentroid":
            for value in content_list[key]:
                times.append(np.datetime64(value["timeOfDopplerCentroidEstimate"]))
                Ambiguity["values"].append(int(value["dopplerAmbiguity"]))
                AmbiguityConfidence["values"].append(float(value["dopplerAmbiguityConfidence"]))
                CentroidReferenceTime["values"].append(float(value["dopplerCentroidReferenceTime"]["#text"]))
                CentroidReferenceTime["attr"]["units"] = value["dopplerCentroidReferenceTime"]["@units"]
                CentroidPolynomialPeriod["values"].append(float(value["dopplerCentroidPolynomialPeriod"]["#text"]))
                CentroidPolynomialPeriod["attr"]["units"] = value["dopplerCentroidPolynomialPeriod"]["@units"]
                CentroidCoefficients["values"].append([float(x) for x in value["dopplerCentroidCoefficients"].split(" ")])
                CentroidConfidence["values"].append(float(value["dopplerCentroidConfidence"]))
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


def create_dataset_doppler_centroid(ds_attr, times, Ambiguity, AmbiguityConfidence, CentroidReferenceTime, CentroidPolynomialPeriod, CentroidCoefficients, CentroidConfidence):
    ds = xr.Dataset()
    ambiguity_da = xr.DataArray(data=Ambiguity["values"], coords={"timeOfDopplerCentroidEstimate": times}, dims=["timeOfDopplerCentroidEstimate"])
    ambiguityConfidence_da = xr.DataArray(data=AmbiguityConfidence["values"], coords={"timeOfDopplerCentroidEstimate": times}, dims=["timeOfDopplerCentroidEstimate"])
    centroidReferenceTime_da = xr.DataArray(data=CentroidReferenceTime["values"], coords={"timeOfDopplerCentroidEstimate": times}, dims=["timeOfDopplerCentroidEstimate"], attrs=CentroidReferenceTime["attr"])
    centroidPolynomialPeriod_da = xr.DataArray(data=CentroidPolynomialPeriod["values"], coords={"timeOfDopplerCentroidEstimate": times}, dims=["timeOfDopplerCentroidEstimate"], attrs=CentroidPolynomialPeriod["attr"])
    centroidCoefficients_da = xr.DataArray(data=np.array(CentroidCoefficients["values"]), coords={"timeOfDopplerCentroidEstimate": times, "n-Coefficients": [i for i in range(np.array(CentroidCoefficients["values"]).shape[1])]}, dims=["timeOfDopplerCentroidEstimate", "n-Coefficients"])
    centroidConfidence_da = xr.DataArray(data=CentroidConfidence["values"], coords={"timeOfDopplerCentroidEstimate": times}, dims=["timeOfDopplerCentroidEstimate"])
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
        "values": []
    }
    for key in content_list:
        if key == "dopplerRateValues":
            if isinstance(content_list[key], dict):
                RateReferenceTime["values"].append(float(content_list[key]["dopplerRateReferenceTime"]["#text"]))
                RateReferenceTime["attr"]["RateReferenceTime units"] = content_list[key]["dopplerRateReferenceTime"]["@units"]
                RateValuesCoefficients["values"].append([float(x) for x in content_list[key]["dopplerRateValuesCoefficients"].split(" ")])
            elif isinstance(content_list[key], list):
                for value in content_list[key]:
                    RateReferenceTime["values"].append(float(content_list[key]["dopplerRateReferenceTime"]["#text"]))
                    RateReferenceTime["attr"]["units"] = content_list[key]["dopplerRateReferenceTime"]["@units"]
                    RateValuesCoefficients["values"].append(
                        [float(x) for x in content_list[key]["dopplerRateValuesCoefficients"].split(" ")])
        elif len(RateReferenceTime["values"]) != 0:
            break
    return {
        "ds_attr": ds_attr,
        "dopplerRateReferenceTime": RateReferenceTime,
        "dopplerRateValuesCoefficients": RateValuesCoefficients
    }


def create_dataset_doppler_rate_values(ds_attr, rateTime, rateCoefficients):
    ds = xr.Dataset()
    rateCoefficients_da = xr.DataArray(data=np.array(rateCoefficients["values"]),
                                       coords={"dopplerRateReferenceTime": rateTime["values"],
                                               "n-RateValuesCoefficients":
                                                   [i for i in range(np.array(rateCoefficients["values"]).shape[1])]},
                                       dims=["dopplerRateReferenceTime", "n-RateValuesCoefficients"], attrs=rateTime["attr"])
    ds["dopplerRateValues"] = rateCoefficients_da
    ds.attrs = ds_attr
    return ds


def create_matrix_data_with_line_and_pix(lines, pixs, vals):
    height = len(np.unique(lines))
    width = len(np.unique(pixs))
    tab = np.ones((width, height)) * np.nan
    unique_lines = np.unique(lines)
    unique_pixs = np.unique(pixs)
    indexs_lines = [np.where(li == unique_lines)[0][0] for li in lines]
    indexs_pixs = [np.where(pi == unique_pixs)[0][0] for pi in pixs]
    for k in range(len(pixs)):
        tab[indexs_lines[k], indexs_pixs[k]] = vals[k]
    return np.array(tab)


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
                        dims=['line', "pixel"], attrs={"units": unit, "xpath": xpath})


def fill_image_attribute(dictio):
    #xpath ="/product/imageAttributes"
    xpath = xpath_dict["geolocation_grid"]["xpath"].split("/geographicInformation")[0]
    content_list = xpath_get(dictio, xpath)
    attr = {
        "rasterAttributes": {}
    }
    for key in content_list:
        if isinstance(content_list[key], str):
            attr[key] = content_list[key]
        elif key == "rasterAttributes":
            for value in content_list[key]:
                if isinstance(content_list[key][value], str):
                    attr[key][value] = content_list[key][value]
                elif isinstance(content_list[key][value], dict):
                    dico_keys = list(content_list[key][value].keys())
                    dico = {}
                    for k in dico_keys:
                        dico[k.replace("@", "").replace("#text", "value")] = content_list[key][value][k]
                    #dico = {"units": content_list[key][value]["@units"], "value": content_list[key][value]["#text"]}
                    attr[key][value] = dico
    return attr


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
    da_los = create_data_array_geolocation_grid(create_matrix_data_with_line_and_pix(lines, pixs, los), "longitude",
                                                lines, pixs, units[0])
    da_las = create_data_array_geolocation_grid(create_matrix_data_with_line_and_pix(lines, pixs, las), "latitude",
                                                lines, pixs, units[1])
    da_hes = create_data_array_geolocation_grid(create_matrix_data_with_line_and_pix(lines, pixs, hes), "height",
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
    """dt = datatree.DataTree.from_dict({
        "geolocationGrid": ds_geo, 
        "orbitAndAttitude": 
            {
                "orbitInformation": ds_orbit_info,
                "attitudeInformation": ds_attitude_info
            }})"""
    dt = datatree.DataTree()
    dt["orbitAndAttitude"] = datatree.DataTree.from_dict({"orbitInformation": ds_orbit_info, "attitudeInformation": ds_attitude_info})
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
    dt["doppler/dopplerCentroid"] = datatree.DataTree(data=ds_doppler_centroid)
    dt["imageAttributes"].attrs = fill_image_attribute(dic)
    dic_doppler_rate_values = get_dic_doppler_rate_values(dic)
    ds_doppler_rate_values = create_dataset_doppler_rate_values(dic_doppler_rate_values["ds_attr"],
                                                                dic_doppler_rate_values["dopplerRateReferenceTime"],
                                                                dic_doppler_rate_values["dopplerRateValuesCoefficients"])
    dt["doppler/dopplerRateValues"] = datatree.DataTree(data=ds_doppler_rate_values)
    print(dt)
    get_dic_doppler_rate_values(dic)
    return dt


if __name__ == '__main__':
    xml_parser(path)


"""# TODO : create doc to fill documentation automatically ( see example on github --> Antoine messages)"""
"""# TODO: fill  datasets with xsd info"""