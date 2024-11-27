from importlib.metadata import version

from xradarsat2.radarSat2_xarray_reader import (
    load_digital_number,  # noqa: F401
    rs2_reader,  # noqa: F401
)

try:
    __version__ = version("xradarsat2")
except Exception:
    __version__ = "999"
