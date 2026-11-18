"""End-to-end tests of the resampled read, on the test product used by the CI.

RS2 ScanSAR Wide 2022-04-07, 10 277 lines x 10 618 samples of 50 m, VV and VH.
At 1000 m, a cell is 20 x 20 pixels: 513 x 530 cells.
The swath is tilted in the image: the corners are outside the swath (zeros),
and many cells straddle its border.

The cells read by `load_digital_number` are compared with the pixels of the
imagery files, recomputed with numpy."""

import numpy as np
import pytest
import rasterio
from rasterio.windows import Window

import xradarsat2
from xradarsat2.utils import get_test_file, load_config

CHUNKS = {"line": 5000, "sample": 5000}
CELL = 20  # pixels per cell at 1000 m


@pytest.fixture(scope="module")
def product():
    return get_test_file(load_config()["folder_path"])


def load(product, resolution):
    dt = xradarsat2.rs2_reader(product)
    dt = xradarsat2.load_digital_number(dt, resolution=resolution, chunks=dict(CHUNKS))
    return dt["digital_numbers"].to_dataset()


@pytest.fixture(scope="module")
def at_1000m(product):
    return load(product, "1000m")


def pixels_of_cell(product, pol, i, j):
    """The 20 x 20 full resolution pixels of cell (i, j)."""
    with rasterio.open(f"{product}/imagery_{pol}.tif") as src:
        return src.read(1, window=Window(CELL * j, CELL * i, CELL, CELL))


def rms_of_valid_pixels(pixels):
    valid = pixels[pixels > 0].astype(np.float64)
    return np.sqrt(np.mean(valid**2)) if valid.size else 0


def test_full_resolution_is_the_file_itself(product):
    """Without resolution, the pixels are the ones of the file, in uint16, no valid_pixel_fraction."""
    ds = load(product, None)

    assert "valid_pixel_fraction" not in ds
    assert ds.sizes["line"] == 10277
    assert ds.sizes["sample"] == 10618
    assert ds.digital_number.dtype == np.uint16
    with rasterio.open(f"{product}/imagery_VV.tif") as src:
        pixels = src.read(1, window=Window(5000, 5000, 50, 50))
    np.testing.assert_array_equal(
        ds.digital_number.sel(
            pol="VV", line=slice(5000, 5049), sample=slice(5000, 5049)
        ),
        pixels,
    )


def test_grid_at_1000m(at_1000m):
    """10 277 x 10 618 pixels of 50 m at 1000 m: 513 x 530 whole cells of 20 x 20 pixels.

    Cell coordinates are the cell centers: 9.5 for the first cell, and
    512 x 20 + 9.5 = 10 249.5 (line), 529 x 20 + 9.5 = 10 589.5 (sample) for the last.
    The 17 last lines and 18 last samples (less than one cell) are not used."""
    ds = at_1000m

    assert ds.sizes["line"] == 513
    assert ds.sizes["sample"] == 530
    assert ds.line[0] == 9.5
    assert ds.line[-1] == 10249.5
    assert ds.sample[0] == 9.5
    assert ds.sample[-1] == 10589.5
    assert ds.digital_number.dtype == np.float32
    assert ds.valid_pixel_fraction.dtype == np.float32
    assert ds.valid_pixel_fraction.dims == ("pol", "line", "sample")
    assert ds.valid_pixel_fraction.attrs["units"] == "1"


def test_four_cells_checked_by_hand(at_1000m):
    """Four cells of the scene, with their values (rms of the valid pixels):

    * corner cell (0, 0), outside the swath: 0, valid_pixel_fraction 0;
    * last cell (512, 529), outside the swath too: 0, valid_pixel_fraction 0;
    * cell in the middle of the scene (256, 265), full: VV 1796.11, VH 1014.68;
    * cell (200, 6) on the swath border, about 10 columns of zeros then 10 of sea.
      The border is not exactly the same in both polarisations:
      VV 200 valid pixels out of 400 (0.5), VV 680.41;
      VH 215 valid pixels out of 400 (0.5375), VH 192.28."""
    ds = at_1000m
    vv = ds.digital_number.sel(pol="VV").values
    vh = ds.digital_number.sel(pol="VH").values
    fraction_vv = ds.valid_pixel_fraction.sel(pol="VV").values
    fraction_vh = ds.valid_pixel_fraction.sel(pol="VH").values

    assert vv[0, 0] == 0
    assert fraction_vv[0, 0] == 0
    assert vv[512, 529] == 0
    assert fraction_vv[512, 529] == 0
    np.testing.assert_allclose(vv[256, 265], 1796.11, atol=0.01)
    np.testing.assert_allclose(vh[256, 265], 1014.68, atol=0.01)
    np.testing.assert_allclose(fraction_vv[256, 265], 1)
    np.testing.assert_allclose(vv[200, 6], 680.41, atol=0.01)
    np.testing.assert_allclose(fraction_vv[200, 6], 0.5)
    np.testing.assert_allclose(vh[200, 6], 192.28, atol=0.01)
    np.testing.assert_allclose(fraction_vh[200, 6], 0.5375)


def test_any_cell_is_the_rms_of_its_own_pixels(at_1000m, product):
    """Same check on 20 cells straddling the swath border and 20 cells taken at random:
    each cell is the rms of the non-zero pixels it covers, and valid_pixel_fraction
    the fraction of these pixels, in each polarisation."""
    ds = at_1000m
    fraction = ds.valid_pixel_fraction.sel(pol="VV").values
    n, m = fraction.shape
    rng = np.random.default_rng(0)
    border = np.argwhere((fraction > 0) & (fraction < 1))
    cells = [tuple(c) for c in border[rng.choice(len(border), 20, replace=False)]]
    cells += [(rng.integers(n), rng.integers(m)) for _ in range(20)]

    for pol in ["VV", "VH"]:
        dn = ds.digital_number.sel(pol=pol).values
        fraction = ds.valid_pixel_fraction.sel(pol=pol).values
        for i, j in cells:
            pixels = pixels_of_cell(product, pol, i, j)
            np.testing.assert_allclose(
                dn[i, j],
                rms_of_valid_pixels(pixels),
                rtol=1e-5,
                err_msg=f"{pol} cell {(i, j)}",
            )
            np.testing.assert_allclose(
                fraction[i, j],
                np.mean(pixels > 0),
                atol=1e-6,
                err_msg=f"{pol} cell {(i, j)}",
            )


@pytest.mark.parametrize(
    "resolution",
    [1000, {"line": 20, "sample": 20}],
    ids=["1000_is_in_meters", "dict_is_in_pixels"],
)
def test_1000m_can_be_written_three_ways(product, at_1000m, resolution):
    """ "1000m", 1000 (meters) and {"line": 20, "sample": 20} (pixels) give the same cells."""
    ds = load(product, resolution)

    np.testing.assert_array_equal(ds.line, at_1000m.line)
    np.testing.assert_array_equal(ds.sample, at_1000m.sample)
    np.testing.assert_array_equal(ds.digital_number, at_1000m.digital_number)
    np.testing.assert_array_equal(
        ds.valid_pixel_fraction, at_1000m.valid_pixel_fraction
    )
