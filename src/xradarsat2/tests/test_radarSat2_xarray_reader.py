import unittest
from unittest.mock import Mock, patch, mock_open, MagicMock
import numpy as np
import xarray as xr
from datetime import datetime
import os
import tempfile

# Import functions from the module
from xradarsat2.radarSat2_xarray_reader import (
    xpath_get,
    parse_value,
    create_2d_matrix,
    get_line_and_pix_info,
    fill_image_attribute,
    sort_list_files_and_get_pols,
    get_glob,
    get_product_attributes,
    get_satellite_height,
    get_satellite_pass_direction,
    create_dic_geolocation_grid,
    get_dic_orbit_information,
    get_dic_attitude_info,
    get_dict_doppler_centroid,
    get_dic_doppler_rate_values,
    get_dict_chirp,
    get_dict_radar_parameters,
)


class TestXpathGet(unittest.TestCase):
    """Test cases for xpath_get function"""
    
    def test_xpath_get_valid_path(self):
        test_dict = {"product": {"sourceAttributes": {"satellite": "RS2"}}}
        result = xpath_get(test_dict, "/product/sourceAttributes/satellite")
        self.assertEqual(result, "RS2")
    
    def test_xpath_get_nested_dict(self):
        test_dict = {"a": {"b": {"c": "value"}}}
        result = xpath_get(test_dict, "/a/b/c")
        self.assertEqual(result, "value")
    
    def test_xpath_get_invalid_path(self):
        test_dict = {"product": {"sourceAttributes": {"satellite": "RS2"}}}
        result = xpath_get(test_dict, "/product/invalid/path")
        self.assertIsNone(result)
    
    def test_xpath_get_empty_dict(self):
        result = xpath_get({}, "/some/path")
        self.assertIsNone(result)


class TestParseValue(unittest.TestCase):
    """Test cases for parse_value function"""
    
    def test_parse_integer(self):
        self.assertEqual(parse_value("123"), 123)
    
    def test_parse_float(self):
        self.assertEqual(parse_value("123.456"), 123.456)
    
    def test_parse_boolean(self):
        self.assertTrue(parse_value("True"))
        self.assertFalse(parse_value("False"))
    
    def test_parse_list(self):
        self.assertEqual(parse_value("[1, 2, 3]"), [1, 2, 3])
    
    def test_parse_string(self):
        self.assertEqual(parse_value("hello world"), "hello world")
    
    def test_parse_dict(self):
        result = parse_value("{'key': 'value'}")
        self.assertEqual(result, {'key': 'value'})


class TestCreate2DMatrix(unittest.TestCase):
    """Test cases for create_2d_matrix function"""
    
    def test_create_simple_matrix(self):
        lines = [0, 0, 1, 1]
        cols = [0, 1, 0, 1]
        vals = [1.0, 2.0, 3.0, 4.0]
        result = create_2d_matrix(lines, cols, vals)
        expected = np.array([[1.0, 2.0], [3.0, 4.0]])
        np.testing.assert_array_equal(result, expected)
    
    def test_create_matrix_with_gaps(self):
        lines = [0, 0, 2, 2]
        cols = [0, 1, 0, 1]
        vals = [1.0, 2.0, 3.0, 4.0]
        result = create_2d_matrix(lines, cols, vals)
        self.assertEqual(result.shape, (2, 2))
    
    def test_create_matrix_unsorted(self):
        lines = [1, 0, 1, 0]
        cols = [1, 0, 0, 1]
        vals = [4.0, 1.0, 3.0, 2.0]
        result = create_2d_matrix(lines, cols, vals)
        expected = np.array([[1.0, 2.0], [3.0, 4.0]])
        np.testing.assert_array_equal(result, expected)


class TestGetLineAndPixInfo(unittest.TestCase):
    """Test cases for get_line_and_pix_info function"""
    
    def test_get_line_and_pix_info(self):
        dictio = {
            "PixelSpacing": 12.5,
            "LineSpacing": 10.0,
            "someOtherKey": "value",
            "anotherPixelSpacing": 15.0
        }
        result = get_line_and_pix_info(dictio)
        self.assertIn("line", result)
        self.assertIn("pixel", result)
        self.assertIn("LineSpacing", result["line"])
        self.assertIn("PixelSpacing", result["pixel"])
        self.assertIn("anotherPixelSpacing", result["pixel"])


class TestFillImageAttribute(unittest.TestCase):
    """Test cases for fill_image_attribute function"""
    
    def test_fill_image_attribute_basic(self):
        dictio = {
            "product": {
                "imageAttributes": {
                    "someKey": "someValue",
                    "rasterAttributes": {
                        "numberOfLines": "1000",
                        "numberOfSamplesPerLine": "2000"
                    }
                }
            }
        }
        result = fill_image_attribute(dictio)
        self.assertIn("rasterAttributes_numberOfLines", result)
        self.assertIn("rasterAttributes_numberOfSamplesPerLine", result)


class TestSortListFilesAndGetPols(unittest.TestCase):
    """Test cases for sort_list_files_and_get_pols function"""
    
    def test_sort_cross_pol_first(self):
        list_tiff = ["/path/to/imagery_HV.tif", "/path/to/imagery_HH.tif"]
        sorted_files, pols = sort_list_files_and_get_pols(list_tiff)
        self.assertEqual(pols, ["HH", "HV"])
        self.assertTrue(sorted_files[0].endswith("HH.tif"))
    
    def test_single_pol(self):
        list_tiff = ["/path/to/imagery_VV.tif"]
        sorted_files, pols = sort_list_files_and_get_pols(list_tiff)
        self.assertEqual(pols, ["VV"])
        self.assertEqual(len(sorted_files), 1)
    
    def test_dual_pol_already_sorted(self):
        list_tiff = ["/path/to/imagery_HH.tif", "/path/to/imagery_HV.tif"]
        sorted_files, pols = sort_list_files_and_get_pols(list_tiff)
        self.assertEqual(pols, ["HH", "HV"])


class TestGetGlob(unittest.TestCase):
    """Test cases for get_glob function"""
    
    def test_get_glob_identical(self):
        strlist = ["file.txt", "file.txt"]
        result = get_glob(strlist)
        self.assertEqual(result, "file.txt")
    
    def test_get_glob_different_chars(self):
        strlist = ["file1.txt", "file2.txt", "file3.txt"]
        result = get_glob(strlist)
        self.assertIn("*", result)
    
    def test_get_glob_prefix_suffix(self):
        strlist = ["prefix_001_suffix", "prefix_002_suffix", "prefix_003_suffix"]
        result = get_glob(strlist)
        self.assertTrue(result.startswith("prefix"))
        self.assertTrue(result.endswith("suffix"))


class TestGetProductAttributes(unittest.TestCase):
    """Test cases for get_product_attributes function"""
    
    def test_get_product_attributes(self):
        dic = {
            "product": {
                "sourceAttributes": {
                    "satellite": "RADARSAT-2",
                    "inputDatasetId": "12345",
                    "rawDataStartTime": "2020-01-01T12:00:00.000000Z",
                    "radarParameters": {}
                }
            }
        }
        result = get_product_attributes(dic)
        self.assertEqual(result["satellite"], "RADARSAT-2")
        self.assertEqual(result["inputDatasetId"], "12345")
        self.assertIn("rawDataStartTime", result)
        # Check if it's a numpy datetime64 (can be array or scalar)
        self.assertTrue(isinstance(result["rawDataStartTime"], (np.datetime64, np.ndarray)))
        if isinstance(result["rawDataStartTime"], np.ndarray):
            self.assertEqual(result["rawDataStartTime"].dtype, np.dtype('datetime64[ns]'))


class TestGetSatelliteHeight(unittest.TestCase):
    """Test cases for get_satellite_height function"""
    
    def test_get_satellite_height(self):
        dic = {
            "product": {
                "imageGenerationParameters": {
                    "sarProcessingInformation": {
                        "satelliteHeight": {
                            "#text": "798000.0",
                            "@units": "m"
                        }
                    }
                }
            }
        }
        result = get_satellite_height(dic)
        self.assertIn("satelliteHeight", result)
        self.assertEqual(result["satelliteHeight"], 798000.0)
        self.assertIn("satelliteHeight_units", result)
        self.assertEqual(result["satelliteHeight_units"], "m")


class TestGetSatellitePassDirection(unittest.TestCase):
    """Test cases for get_satellite_pass_direction function"""
    
    def test_get_pass_direction_ascending(self):
        dic = {
            "product": {
                "sourceAttributes": {
                    "orbitAndAttitude": {
                        "orbitInformation": {
                            "passDirection": "Ascending"
                        }
                    }
                }
            }
        }
        result = get_satellite_pass_direction(dic)
        self.assertEqual(result["passDirection"], "Ascending")
    
    def test_get_pass_direction_descending(self):
        dic = {
            "product": {
                "sourceAttributes": {
                    "orbitAndAttitude": {
                        "orbitInformation": {
                            "passDirection": "Descending"
                        }
                    }
                }
            }
        }
        result = get_satellite_pass_direction(dic)
        self.assertEqual(result["passDirection"], "Descending")


class TestCreateDicGeolocationGrid(unittest.TestCase):
    """Test cases for create_dic_geolocation_grid function"""
    
    def test_create_dic_geolocation_grid(self):
        dictio = {
            "product": {
                "imageAttributes": {
                    "geographicInformation": {
                        "geolocationGrid": {
                            "imageTiePoint": [
                                {
                                    "imageCoordinate": {
                                        "line": "0",
                                        "pixel": "0"
                                    },
                                    "geodeticCoordinate": {
                                        "latitude": {"#text": "45.0", "@units": "deg"},
                                        "longitude": {"#text": "-75.0", "@units": "deg"},
                                        "height": {"#text": "0.0", "@units": "m"}
                                    }
                                }
                            ]
                        }
                    },
                    "rasterAttributes": {
                        "numberOfLines": "100",
                        "numberOfSamplesPerLine": "100"
                    }
                }
            }
        }
        result = create_dic_geolocation_grid(dictio)
        self.assertIn("latitude", result)
        self.assertIn("longitude", result)
        self.assertIn("height", result)
        self.assertIn("coords", result)
        self.assertEqual(len(result["latitude"]["values"]), 1)
        self.assertEqual(result["latitude"]["values"][0], 45.0)


class TestGetDicOrbitInformation(unittest.TestCase):
    """Test cases for get_dic_orbit_information function"""
    
    def test_get_dic_orbit_information(self):
        dictio = {
            "product": {
                "sourceAttributes": {
                    "orbitAndAttitude": {
                        "orbitInformation": {
                            "passDirection": "Ascending",
                            "stateVector": [
                                {
                                    "timeStamp": "2020-01-01T12:00:00.000000Z",
                                    "xPosition": {"#text": "1000000.0", "@units": "m"},
                                    "yPosition": {"#text": "2000000.0", "@units": "m"},
                                    "zPosition": {"#text": "3000000.0", "@units": "m"},
                                    "xVelocity": {"#text": "100.0", "@units": "m/s"},
                                    "yVelocity": {"#text": "200.0", "@units": "m/s"},
                                    "zVelocity": {"#text": "300.0", "@units": "m/s"}
                                }
                            ]
                        }
                    }
                }
            }
        }
        result = get_dic_orbit_information(dictio)
        self.assertIn("ds_attr", result)
        self.assertIn("timestamp", result)
        self.assertIn("xPosition", result)
        self.assertEqual(len(result["xPosition"]["values"]), 1)
        self.assertEqual(result["xPosition"]["values"][0], 1000000.0)


class TestGetDicAttitudeInfo(unittest.TestCase):
    """Test cases for get_dic_attitude_info function"""
    
    def test_get_dic_attitude_info(self):
        dictio = {
            "product": {
                "sourceAttributes": {
                    "orbitAndAttitude": {
                        "attitudeInformation": {
                            "attitudeAngles": [
                                {
                                    "timeStamp": "2020-01-01T12:00:00.000000Z",
                                    "yaw": {"#text": "0.5", "@units": "deg"},
                                    "roll": {"#text": "0.3", "@units": "deg"},
                                    "pitch": {"#text": "0.2", "@units": "deg"}
                                }
                            ]
                        }
                    }
                }
            }
        }
        result = get_dic_attitude_info(dictio)
        self.assertIn("ds_attr", result)
        self.assertIn("timestamp", result)
        self.assertIn("yaw", result)
        self.assertIn("roll", result)
        self.assertIn("pitch", result)
        self.assertEqual(len(result["yaw"]["values"]), 1)


class TestGetDictDopplerCentroid(unittest.TestCase):
    """Test cases for get_dict_doppler_centroid function"""
    
    def test_get_dict_doppler_centroid(self):
        dictio = {
            "product": {
                "imageGenerationParameters": {
                    "dopplerCentroid": [
                        {
                            "timeOfDopplerCentroidEstimate": "2020-01-01T12:00:00.000000Z",
                            "dopplerAmbiguity": "0",
                            "dopplerAmbiguityConfidence": "1.0",
                            "dopplerCentroidReferenceTime": {"#text": "0.5", "@units": "s"},
                            "dopplerCentroidPolynomialPeriod": {"#text": "1.0", "@units": "s"},
                            "dopplerCentroidCoefficients": "1.0 2.0 3.0",
                            "dopplerCentroidConfidence": "0.95"
                        }
                    ]
                }
            }
        }
        result = get_dict_doppler_centroid(dictio)
        self.assertIn("timeOfDopplerCentroidEstimate", result)
        self.assertIn("dopplerAmbiguity", result)
        self.assertEqual(len(result["dopplerCentroidCoefficients"]["values"]), 1)
        self.assertEqual(len(result["dopplerCentroidCoefficients"]["values"][0]), 3)


class TestGetDicDopplerRateValues(unittest.TestCase):
    """Test cases for get_dic_doppler_rate_values function"""
    
    def test_get_dic_doppler_rate_values_dict(self):
        dictio = {
            "product": {
                "imageGenerationParameters": {
                    "dopplerRateValues": {
                        "dopplerRateReferenceTime": {"#text": "0.5", "@units": "s"},
                        "dopplerRateValuesCoefficients": "1.0 2.0 3.0"
                    }
                }
            }
        }
        result = get_dic_doppler_rate_values(dictio)
        self.assertIn("dopplerRateReferenceTime", result)
        self.assertIn("dopplerRateValuesCoefficients", result)
        self.assertEqual(len(result["dopplerRateValuesCoefficients"]["values"]), 1)


class TestGetDictRadarParameters(unittest.TestCase):
    """Test cases for get_dict_radar_parameters function"""
    
    def test_get_dict_radar_parameters(self):
        dictio = {
            "product": {
                "sourceAttributes": {
                    "radarParameters": {
                        "acquisitionType": "Wide",
                        "beams": "W1 W2",
                        "polarizations": "HH HV",
                        "radarCenterFrequency": "5405000000.0"
                    }
                }
            }
        }
        result = get_dict_radar_parameters(dictio)
        self.assertIn("ds_attr", result)
        self.assertEqual(result["ds_attr"]["acquisitionType"], "Wide")
        self.assertEqual(result["ds_attr"]["beams"], ["W1", "W2"])
        self.assertEqual(result["ds_attr"]["polarizations"], ["HH", "HV"])


if __name__ == "__main__":
    unittest.main()
