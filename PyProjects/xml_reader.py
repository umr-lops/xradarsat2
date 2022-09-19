# imports
import numpy as np
from lxml import etree
from lxml import objectify
import xmltodict
import xarray as xr
import time

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
            "xpath": "/product/sourceAttributes/orbitAndAttitude/orbitInformation"
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


def get_lists_orbit_information(dictio):
    content_list = xpath_get(dictio, xpath_dict["orbit_information"]["xpath"])
    ds_attr = {}
    timestamp = []
    xPosition = {
        "values": [],
        "unit": ""
    }
    yPosition = {
        "values": [],
        "unit": ""
    }
    zPosition = {
        "values": [],
        "unit": ""
    }
    xVelocity = {
        "values": [],
        "unit": ""
    }
    yVelocity = {
        "values": [],
        "unit": ""
    }
    zVelocity = {
        "values": [],
        "unit": ""
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
                timestamp.append(value["timeStamp"])
                xPosition["values"].append(value["xPosition"]["#text"])
                xPosition["unit"] = value["xPosition"]["@units"]
                yPosition["values"].append(value["yPosition"]["#text"])
                yPosition["unit"] = value["yPosition"]["@units"]
                zPosition["values"].append(value["zPosition"]["#text"])
                zPosition["unit"] = value["zPosition"]["@units"]
                xVelocity["values"].append(value["xVelocity"]["#text"])
                xVelocity["unit"] = value["xVelocity"]["@units"]
                yVelocity["values"].append(value["yVelocity"]["#text"])
                yVelocity["unit"] = value["yVelocity"]["@units"]
                zVelocity["values"].append(value["zVelocity"]["#text"])
                zVelocity["unit"] = value["zVelocity"]["@units"]
    return ds_attr, timestamp, xPosition, yPosition, zPosition, xVelocity, yVelocity, zVelocity


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


"""
def recursive_dataset_filling_geolocation_grid(dictio, lines, pixs, los, las, hes, lo_unit, la_unit, he_unit, units):
    for keys in dictio:
        if len(lines) != 0:
            break
        elif isinstance(dictio[keys], dict) and (len(lines) == 0):
            recursive_dataset_filling_geolocation_grid(dictio[keys], lines, pixs, los, las, hes, lo_unit, la_unit, he_unit, units)

        while keys == 'imageTiePoint':
            for coord_dict in dictio[keys]:
                line = coord_dict['imageCoordinate']['line']
                pix = coord_dict['imageCoordinate']['pixel']
                lo_dict = coord_dict['geodeticCoordinate']['longitude']
                if lo_dict['@units'] != "":
                    lo_unit = lo_dict['@units']
                    units[0] = lo_dict['@units']
                lo = lo_dict['#text']
                la_dict = coord_dict['geodeticCoordinate']['latitude']
                if la_dict['@units'] != "":
                    la_unit = la_dict['@units']
                    units[1] = la_dict['@units']
                la = la_dict['#text']
                he_dict = coord_dict['geodeticCoordinate']['height']
                if he_dict['@units'] != "":
                    he_unit = he_dict['@units']
                    units[2] = he_dict['@units']
                he = he_dict['#text']
                lines.append(line)
                pixs.append(pix)
                los.append(lo)
                las.append(la)
                hes.append(he)
            break
    if len(lines) != 0:
        return lines, pixs, los, las, hes, lo_unit, la_unit, he_unit, units
# TODO : debug loop to keep units --> it returns on the recursive when line is full and put an empty string in units...why??
# TODO : insert xsd path and collect documentation for DS attributes
# TODO : history of hierarchy to create an Xpath and put it in DA attributes
# TODO : finally, when all of it is finished, transfer to notebook !!
"""


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
    return ds_geo


if __name__ == '__main__':
    # xml_parser(path)
    with open(path, 'rb') as f:
        xml_content = f.read()
        dic = xmltodict.parse(xml_content)
        f.close()
    get_lists_orbit_information(dic)
