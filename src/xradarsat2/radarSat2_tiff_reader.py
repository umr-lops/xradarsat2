import glob
import os

import dask
import numpy as np
import rasterio
import rioxarray
import xarray as xr
import yaml
from affine import Affine

from xradarsat2.utils import get_glob, load_config

# folder_path = "/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV/2010/288/" \
#              "RS2_OK72200_PK649463_DK111111_SCWA_20101015_210132_VV_SGF"
folder_path = (
    "/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV_VH/2021/137"
    "/RS2_OK129673_PK1136693_DK1093025_SCWA_20210517_010235_VV_VH_SGF"
)
conf = load_config()
folder_path = conf["folder_path"]


def list_tiff_files(root_path):
    return glob.glob(os.path.join(root_path, "*imagery*tif"))


def _load_digital_number(
    root_path,
    dt,
    resolution=None,
    chunks=None,
    resampling=rasterio.enums.Resampling.rms,
):
    """

    :param root_path:
    :param dt:
    :param resolution:
    :param chunks:
    :param resampling:
    :return:
    """
    tiff_files = list_tiff_files(root_path)
    map_dims = {"pol": "band", "line": "y", "sample": "x"}
    if resolution is not None:
        comment = f'resampled at "{resolution}" with {resampling.__module__}.{resampling.__class__.__name__}.{resampling.name}'
    else:
        comment = "read at full resolution"

    # arbitrary rio object, to get shape, etc ... (will not be used to read data)
    rio = rasterio.open(tiff_files[0])

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
        ).assign_coords(
            band=np.arange(len(dt["imageGenerationParameters"]["chirp"]["pole"].values))
            + 1
        )

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
                for f, pol in zip(
                    tiff_files, dt["imageGenerationParameters"]["chirp"]["pole"].values
                )
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
    dn = dn.assign_coords(pol=dt["imageGenerationParameters"]["chirp"]["pole"].values)

    descr = "not denoised"
    var_name = "digital_number"

    dn.attrs = {
        "comment": f"{descr} digital number, {comment}",
        "history": yaml.safe_dump(
            {var_name: get_glob([p.replace(root_path + "/", "") for p in tiff_files])}
        ),
    }
    ds = dn.to_dataset(name=var_name)
    dt["digital_numbers"] = xr.DataTree(data=ds)
    return dt
