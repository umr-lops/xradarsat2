from xradarsat2.radarSat2_xarray_reader import load_digital_number  # noqa: F401
from xradarsat2.radarSat2_xarray_reader import rs2_reader  # noqa: F401
from importlib.metadata import version

try:
    __version__ = version("xradarsat2")
except Exception:
    __version__ = "999"
