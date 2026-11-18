"""Resampled read (`resolution` given), on small synthetic GeoTIFFs.

Each test is a concrete case met with RadarSAT-2 products, with the expected values
computed by hand in the docstring:

* a ScanSAR Wide image (50 m pixels) of 10 277 lines read at 1000 m: the last cell
  must hold its own lines;
* the same pixels read at 125 m: cells of 2.5 pixels;
* a cell straddling the swath border: the zeros outside the swath must not darken it;
* a low wind cross-pol cell: the mean must not be rounded to an integer.
"""

import warnings

import numpy as np
import pytest
import rasterio
from rasterio.enums import Resampling

from xradarsat2.radarSat2_xarray_reader import read_cells, read_valid_fraction

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)


def write_tiff(path, array):
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=array.shape[0],
        width=array.shape[1],
        count=1,
        dtype=array.dtype,
    ) as dst:
        dst.write(array, 1)
    return str(path)


def line_numbers(height, width):
    """Image whose pixels hold their line number + 1 (0 would be read as outside the swath)."""
    return np.repeat(np.arange(1, height + 1, dtype=np.uint16)[:, None], width, axis=1)


def test_last_cell_of_a_scansar_image_holds_its_own_lines(tmp_path):
    """ScanSAR Wide image of 10 277 lines of 50 m, read at 1000 m: 513 cells of 20 lines.

    Cell k holds lines 20 k to 20 k + 19, so with pixels = line number + 1:
    first cell 10.5, last cell (lines 10 240 to 10 259) 10 250.5.
    The 17 last lines (less than one cell) are not used."""
    path = write_tiff(tmp_path / "scw.tif", line_numbers(10277, 20))
    resolution = {"line": 20.0, "sample": 20.0}

    cells = read_cells(path, resolution, (513, 1), Resampling.average)[0]

    assert cells[0, 0] == 10.5
    assert cells[-1, 0] == 10250.5


def test_scansar_image_read_at_125m_has_cells_of_two_and_a_half_pixels(tmp_path):
    """50 m pixels read at 125 m: 2.5 pixels per cell.

    The pixels hold their sample number + 1. The pixel cut by a cell edge counts
    for the part of it inside the cell:
    cell 0 = samples 0, 1 and half of 2 = (1 + 2 + 0.5 x 3) / 2.5 = 1.8
    cell 1 = half of 2, 3 and 4 = (0.5 x 3 + 4 + 5) / 2.5 = 4.2
    cell 2 = 5, 6 and half of 7 = (6 + 7 + 0.5 x 8) / 2.5 = 6.8

    1003 samples give 401 cells, the last one (cell 400) on samples 1000, 1001 and
    half of 1002 = (1001 + 1002 + 0.5 x 1003) / 2.5 = 1001.8."""
    path = write_tiff(tmp_path / "scw.tif", line_numbers(1003, 10).T.copy())
    resolution = {"line": 2.5, "sample": 2.5}

    cells = read_cells(path, resolution, (4, 401), Resampling.average)[0]

    np.testing.assert_allclose(cells[0, :3], [1.8, 4.2, 6.8], rtol=1e-6)
    np.testing.assert_allclose(cells[0, -1], 1001.8, rtol=1e-6)


@pytest.mark.parametrize("border", ["range", "azimuth"])
def test_cell_across_the_swath_border_is_not_darkened(tmp_path, border):
    """Sea at DN = 50, and the swath ends at pixel 135 (zeros after it), either in
    range or in azimuth. Cells of 10 pixels.

    Cell 13 (pixels 130 to 139) is half in the swath: its value is the sea value, 50,
    and its valid_pixel_fraction 0.5. Cells 14 to 19 are outside the swath: 0."""
    array = np.full((200, 200), 50, dtype=np.uint16)
    array[:, 135:] = 0
    if border == "azimuth":
        array = array.T.copy()
    path = write_tiff(tmp_path / "border.tif", array)
    resolution = {"line": 10.0, "sample": 10.0}

    cells = read_cells(path, resolution, (20, 20), Resampling.rms)[0]
    fraction = read_valid_fraction(path, resolution, (20, 20))
    if border == "azimuth":
        cells = cells.T
        fraction = fraction.T

    np.testing.assert_allclose(cells[:, 13], 50)
    np.testing.assert_allclose(fraction[:, 13], 0.5)
    np.testing.assert_allclose(cells[:, :13], 50)
    np.testing.assert_allclose(fraction[:, :13], 1)
    np.testing.assert_array_equal(cells[:, 14:], 0)
    np.testing.assert_array_equal(fraction[:, 14:], 0)


def test_valid_pixel_fraction_of_a_cell_of_two_and_a_half_pixels_across_the_border(
    tmp_path,
):
    """2.5 pixels per cell, the swath ends at sample 136.

    Cell 54 covers samples 135 to 137.5: pixel 135 is valid, pixel 136 and half of
    pixel 137 are outside the swath. 1 valid pixel out of 2.5: valid_pixel_fraction = 0.4.
    Its value is still the sea value, 50."""
    array = np.full((200, 200), 50, dtype=np.uint16)
    array[:, 136:] = 0
    path = write_tiff(tmp_path / "border.tif", array)
    resolution = {"line": 2.5, "sample": 2.5}

    cells = read_cells(path, resolution, (80, 80), Resampling.rms)[0]
    fraction = read_valid_fraction(path, resolution, (80, 80))

    np.testing.assert_allclose(cells[:, 54], 50)
    np.testing.assert_allclose(fraction[:, 54], 0.4, rtol=1e-6)
    np.testing.assert_allclose(fraction[:, 53], 1)
    np.testing.assert_array_equal(fraction[:, 55], 0)


def test_low_wind_cross_pol_cell_is_not_rounded(tmp_path):
    """Low wind cross-pol over the ocean: DN around 20. Cell of 10 x 10 pixels alternating 20 and 21.

    rms = sqrt((20^2 + 21^2) / 2) = 20.506, kept in float32 (rounding it to 20 would
    lower sigma0 by 0.22 dB)."""
    array = (20 + np.indices((100, 100)).sum(axis=0) % 2).astype(np.uint16)
    path = write_tiff(tmp_path / "vh.tif", array)

    cells = read_cells(path, {"line": 10.0, "sample": 10.0}, (10, 10), Resampling.rms)

    assert cells.dtype == np.float32
    np.testing.assert_allclose(cells, 20.506, atol=1e-3)
