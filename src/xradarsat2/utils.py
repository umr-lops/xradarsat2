import logging
import os
import re
import warnings
import zipfile

import aiohttp
import fsspec
import yaml

import xradarsat2


def get_glob(strlist):
    # from list of str, replace diff by '?'
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


def load_config():
    """

    Returns:
        conf: dict
    """
    local_config_path = os.path.join(
        os.path.dirname(xradarsat2.__file__), "localxradarsat2-config.yaml"
    )

    if os.path.exists(local_config_path):
        config_path = local_config_path
    else:
        config_path = os.path.join(
            os.path.dirname(xradarsat2.__file__), "xradarsat2-config.yaml"
        )

    logging.info("config path: %s", config_path)
    stream = open(config_path)
    conf = yaml.load(stream, Loader=yaml.CLoader)
    return conf


def get_test_file(fname):
    """
    get test file from  https://cyclobs.ifremer.fr/static/sarwing_datarmor/xsardata/
    file is unzipped and extracted to `config['data_dir']`

    Parameters
    ----------
    fname: str
        file name to get (without '.zip' extension)

    Returns
    -------
    str
        path to file, relative to `config['data_dir']`

    """
    config = {"data_dir": "/tmp"}

    def url_get(url, cache_dir=os.path.join(config["data_dir"], "fsspec_cache")):
        """
        Get fil from url, using caching.

        Parameters
        ----------
        url: str
        cache_dir: str
            Cache dir to use. default to `os.path.join(config['data_dir'], 'fsspec_cache')`

        Raises
        ------
        FileNotFoundError

        Returns
        -------
        filename: str
            The local file name

        Notes
        -----
        Due to fsspec, the returned filename won't match the remote one.
        """

        if "://" in url:
            with fsspec.open(
                f"filecache::{url}",
                https={"client_kwargs": {"timeout": aiohttp.ClientTimeout(total=3600)}},
                filecache={
                    "cache_storage": os.path.join(
                        os.path.join(config["data_dir"], "fsspec_cache")
                    )
                },
            ) as f:
                fname = f.name
        else:
            fname = url

        return fname

    res_path = config["data_dir"]
    base_url = "https://cyclobs.ifremer.fr/static/sarwing_datarmor/xsardata"
    file_url = f"{base_url}/{fname}.zip"
    if not os.path.exists(os.path.join(res_path, fname)):
        warnings.warn(f"Downloading {file_url}")
        local_file = url_get(file_url)
        warnings.warn(f"Unzipping {os.path.join(res_path, fname)}")
        with zipfile.ZipFile(local_file, "r") as zip_ref:
            zip_ref.extractall(res_path)
    return os.path.join(res_path, fname)
