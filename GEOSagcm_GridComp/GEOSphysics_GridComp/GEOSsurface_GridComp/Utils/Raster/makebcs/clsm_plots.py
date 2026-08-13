#!/usr/bin/env python3
"""
Python replacement/extension for GEOS makebcs clsm_plots.pro

This script is intended to be run from the same place the IDL driver was run
(usually clsm/plots). With no arguments it reads $gfile, $workdir, $NC, and 
$NR and writes the standard plot products into the current directory.

Typical first run from clsm/plots:

    python clsm_plots.py \
        --gfile "$gfile" --workdir "$workdir" --nc "$NC" --nr "$NR" \
        --plots default --outdir .

Notes:
  * F77 unformatted files are read as sequential records with 4-byte record
    markers by default.  Use --endian/--record-marker if your files differ.
  * Cartopy is optional.  If available, this script can draw coastlines with
    --coastlines.  Otherwise plots are still generated with lon/lat axes.
   * Movie generation is optional and potentially slow. Use --plots movies
    for movies only, or --plots legacy for fixed JPGs plus movies.  
"""
from __future__ import annotations

import argparse
import dataclasses
import datetime as _dt
import glob
import math
import os
import sys
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.ticker import FixedLocator, FixedFormatter

try:
    from scipy import sparse as sp_sparse  # type: ignore
except Exception:  # pragma: no cover - optional on Discover modules
    sp_sparse = None

try:
    from scipy.stats import mode as scipy_mode  # type: ignore
except Exception:  # pragma: no cover
    scipy_mode = None

class FFMpegPipeWriter:
    """Small MP4 writer that pipes RGB frames to the system ffmpeg executable.

    Discover's GEOSpyD imageio package may not include the imageio-ffmpeg plugin,
    even when the shell ffmpeg module is available.  Calling ffmpeg through
    subprocess avoids that Python plugin dependency.  Load the module with

        module load ffmpeg/5.0

    or set CLSM_FFMPEG=/path/to/ffmpeg before running.
    """

    def __init__(self, outpath: Path, fps: int = 10):
        self.outpath = Path(outpath)
        self.fps = int(fps)
        self.proc: Optional[subprocess.Popen] = None
        self.width: Optional[int] = None
        self.height: Optional[int] = None

    def __enter__(self):
        self.outpath.parent.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            if self.proc is not None and self.proc.poll() is None:
                try:
                    self.proc.kill()
                except Exception:
                    pass
            return False
        self.close()
        return False

    @staticmethod
    def _ffmpeg_exe() -> str:
        explicit = os.environ.get("CLSM_FFMPEG")
        if explicit:
            p = Path(explicit).expanduser()
            if p.exists():
                return str(p)
        exe = shutil.which("ffmpeg")
        if exe:
            return exe
        # Useful Discover fallback if the module path is present but PATH was not
        # updated for some reason.  The preferred route is still module load.
        fallback = Path("/usr/local/other/ffmpeg/5.0/bin/ffmpeg")
        if fallback.exists():
            return str(fallback)
        raise ClsmPlotError(
            "ffmpeg executable not found. Load it with `module load ffmpeg/5.0` "
            "or set CLSM_FFMPEG=/path/to/ffmpeg before requesting movies."
        )

    @staticmethod
    def _as_rgb_uint8(frame: np.ndarray) -> np.ndarray:
        arr = np.asarray(frame)
        if arr.ndim == 2:
            arr = np.repeat(arr[:, :, None], 3, axis=2)
        if arr.ndim != 3 or arr.shape[2] not in (3, 4):
            raise ClsmPlotError(f"Movie frame must be HxWx3 or HxWx4, got shape {arr.shape}")
        arr = arr[:, :, :3]
        if arr.dtype != np.uint8:
            arr = arr.astype(np.float32, copy=False)
            if np.nanmax(arr) <= 1.0:
                arr = arr * 255.0
            arr = np.nan_to_num(arr, nan=255.0, posinf=255.0, neginf=0.0)
            arr = np.clip(arr, 0.0, 255.0).astype(np.uint8)

        # H.264 with yuv420p requires even frame dimensions.  Adding the movie
        # colorbar changed the canvas to 780x585 on Discover, which made
        # ffmpeg reject the stream ("height not divisible by 2").  Pad, rather
        # than crop, so no tick labels/colorbar pixels are lost.  Use white
        # padding to match the figure background.
        h, w = arr.shape[:2]
        new_h = h + (h % 2)
        new_w = w + (w % 2)
        if new_h != h or new_w != w:
            padded = np.full((new_h, new_w, 3), 255, dtype=np.uint8)
            padded[:h, :w, :] = arr
            arr = padded
        return np.ascontiguousarray(arr)

    def _start(self, frame: np.ndarray) -> None:
        h, w = frame.shape[:2]
        self.height, self.width = int(h), int(w)
        ffmpeg = self._ffmpeg_exe()
        cmd = [
            ffmpeg,
            "-y",
            "-loglevel", "error",
            "-f", "rawvideo",
            "-vcodec", "rawvideo",
            "-pix_fmt", "rgb24",
            "-s", f"{self.width}x{self.height}",
            "-r", str(self.fps),
            "-i", "-",
            "-an",
            "-vcodec", "libx264",
            "-preset", "medium",
            "-crf", "18",
            "-pix_fmt", "yuv420p",
            str(self.outpath),
        ]
        self.proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

    def append_data(self, frame: np.ndarray) -> None:
        rgb = self._as_rgb_uint8(frame)
        if self.proc is None:
            self._start(rgb)
        if rgb.shape[0] != self.height or rgb.shape[1] != self.width:
            raise ClsmPlotError(
                f"Movie frame size changed from {self.width}x{self.height} "
                f"to {rgb.shape[1]}x{rgb.shape[0]}"
            )
        assert self.proc is not None and self.proc.stdin is not None
        try:
            self.proc.stdin.write(rgb.tobytes())
        except BrokenPipeError as exc:
            err = b""
            if self.proc.stderr is not None:
                err = self.proc.stderr.read()
            raise ClsmPlotError(f"ffmpeg pipe closed while writing {self.outpath}: {err.decode(errors='replace')}") from exc

    def close(self) -> None:
        if self.proc is None:
            return
        assert self.proc.stdin is not None
        self.proc.stdin.close()
        ret = self.proc.wait()
        err = b""
        if self.proc.stderr is not None:
            err = self.proc.stderr.read()
        if ret != 0:
            raise ClsmPlotError(
                f"ffmpeg failed while writing {self.outpath} with exit code {ret}: "
                f"{err.decode(errors='replace')}"
            )


def open_mp4_writer(outpath: Path, fps: int = 10):
    """Open an MP4 writer. Uses system ffmpeg, not imageio-ffmpeg."""
    return FFMpegPipeWriter(Path(outpath), fps=fps)

try:
    import xarray as xr  # type: ignore
except Exception:  # pragma: no cover
    xr = None

try:
    import cartopy.crs as ccrs  # type: ignore
    import cartopy.feature as cfeature  # type: ignore
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter  # type: ignore
except Exception:  # pragma: no cover
    ccrs = None
    cfeature = None
    LongitudeFormatter = None
    LatitudeFormatter = None


# -----------------------------------------------------------------------------
# IDL-compatible color tables
# -----------------------------------------------------------------------------


def _as_rgb(rows: Sequence[Sequence[int]]) -> np.ndarray:
    return np.asarray(rows, dtype=np.float32) / 255.0


def idl_palette() -> np.ndarray:
    """Return a 256 x 3 RGB palette approximating load_colors.pro."""
    rgb = np.ones((256, 3), dtype=np.float32)

    def put(start: int, r: Sequence[int], g: Sequence[int], b: Sequence[int]) -> None:
        n = len(r)
        rgb[start:start + n, 0] = np.asarray(r) / 255.0
        rgb[start:start + n, 1] = np.asarray(g) / 255.0
        rgb[start:start + n, 2] = np.asarray(b) / 255.0

    r_drought = [0, 0, 0, 0, 47, 200, 255, 255, 255, 255, 249, 197]
    g_drought = [0, 115, 159, 210, 255, 255, 255, 255, 219, 157, 0, 0]
    b_drought = [0, 0, 0, 0, 67, 130, 255, 0, 0, 0, 0, 0]
    put(0, r_drought, g_drought, b_drought)

    r_green = [200, 150, 47, 60, 0, 0, 0, 0]
    g_green = [255, 255, 255, 230, 219, 187, 159, 131]
    b_green = [200, 150, 67, 15, 0, 0, 0, 0]
    put(20, r_green, g_green, b_green)

    r_blue = [55, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    g_blue = [255, 255, 227, 195, 167, 115, 83, 0, 0, 0]
    b_blue = [199, 255, 255, 255, 255, 255, 255, 255, 200, 130]
    put(30, r_blue, g_blue, b_blue)

    r_red = [255, 240, 255, 255, 255, 255, 255, 233, 197]
    g_red = [255, 255, 219, 187, 159, 131, 51, 23, 0]
    b_red = [153, 15, 0, 0, 0, 0, 0, 0, 0]
    put(40, r_red, g_red, b_red)

    r_grey = [245, 225, 205, 185, 165, 145, 125, 105, 85]
    g_grey = [245, 225, 205, 185, 165, 145, 125, 105, 85]
    b_grey = [245, 225, 205, 185, 165, 145, 125, 105, 85]
    put(50, r_grey, g_grey, b_grey)

    r_type = [255, 106, 202, 251, 0, 29, 77, 109, 142, 233, 255, 255, 255, 127, 164, 164, 217, 217, 204, 104, 0]
    g_type = [245, 91, 178, 154, 85, 115, 145, 165, 185, 23, 131, 131, 191, 39, 53, 53, 72, 72, 204, 104, 70]
    b_type = [215, 154, 214, 153, 0, 0, 0, 0, 13, 0, 0, 200, 0, 4, 3, 200, 1, 200, 204, 200, 200]
    put(60, r_type, g_type, b_type)

    r_lct2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 55, 120, 190, 240, 255, 255, 255, 255, 255, 233, 197, 158]
    g_lct2 = [0, 0, 0, 83, 115, 167, 195, 227, 255, 255, 255, 255, 255, 219, 187, 159, 131, 51, 23, 0, 0]
    b_lct2 = [130, 200, 255, 255, 255, 255, 255, 255, 255, 199, 135, 67, 15, 0, 0, 0, 0, 0, 0, 0, 0]
    put(140, r_lct2, g_lct2, b_lct2)

    r_veg = [233, 255, 255, 255, 210, 0, 0, 0, 204, 170, 255, 220, 205, 0, 0, 170, 0, 40, 120, 140, 190, 150, 255, 255, 0, 0, 0, 195, 255, 0]
    g_veg = [23, 131, 191, 255, 255, 255, 155, 0, 204, 240, 255, 240, 205, 100, 160, 200, 60, 100, 130, 160, 150, 100, 180, 235, 120, 150, 220, 20, 245, 70]
    b_veg = [0, 0, 0, 178, 255, 255, 255, 200, 204, 240, 100, 100, 102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 175, 90, 120, 130, 0, 215, 200]
    put(90, r_veg, g_veg, b_veg)

    r_grads_rb = [160, 110, 30, 0, 0, 0, 0, 160, 230, 230, 240, 250, 240]
    g_grads_rb = [0, 0, 60, 150, 200, 210, 220, 230, 220, 175, 130, 60, 0]
    b_grads_rb = [200, 220, 255, 255, 200, 140, 0, 50, 50, 45, 40, 60, 130]
    put(120, r_grads_rb, g_grads_rb, b_grads_rb)

    rgb[255] = (1.0, 1.0, 1.0)
    return rgb


PALETTE = idl_palette()
CONTINUOUS_COLOR_IDS = [27, 26, 25, 24, 23, 22, 21, 20, 40, 41, 42, 43, 44, 45, 46, 47, 48]
LAI_RGB = _as_rgb([
    [253, 253, 253], [224, 238, 224], [255, 255, 0], [238, 238, 0],
    [205, 205, 0], [193, 255, 193], [152, 251, 152], [0, 255, 127],
    [124, 252, 0], [0, 255, 0], [0, 238, 0], [0, 205, 0],
    [0, 139, 0], [0, 128, 0], [0, 100, 0], [48, 128, 20],
    [110, 139, 61], [85, 107, 47],
])
# White is reserved for missing / invalid / no-data only.
# It is not used as a valid data color in LAI, GREEN, NDVI, VISDF, or NIRDF.
NO_DATA_COLOR = (1.0, 1.0, 1.0)

# LAI is plotted from 0.0 to 7.5 using 0.5 increments.
# There are 16 boundaries and therefore 15 valid color intervals.
# Drop the first nearly-white color from LAI_RGB so the first valid data
# interval, 0.0-0.5, is light gray instead of white.
LAI_LEVELS = np.asarray([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5,
    3.0, 3.5, 4.0, 4.5, 5.0, 5.5,
    6.0, 6.5, 7.0, 7.5,
], dtype=np.float32)

LAI_PLOT_RGB = LAI_RGB[1:len(LAI_LEVELS)].copy()

LAI_TICKS = LAI_LEVELS.astype(float)

LAI_TICK_LABELS = [
    "0", "0.5", "1", "1.5", "2", "2.5",
    "3", "3.5", "4", "4.5", "5", "5.5",
    "6", "6.5", "7", "7.5",
]

# Fraction fields use FRACTION_LEVELS as true 0..1 bin boundaries.
# FRACTION_LEVELS has 18 boundaries, so it needs 17 colors.
# White is reserved for missing/no-data; valid zero values use light blue.
FRACTION_LEVELS = np.asarray([
    0.00, 0.025, 0.050, 0.075, 0.10,
    0.15, 0.20, 0.25, 0.30, 0.35,
    0.40, 0.45, 0.50,
    0.60, 0.70, 0.80, 0.90, 1.00,
], dtype=np.float32)

FRACTION_RGB = LAI_RGB[1:len(FRACTION_LEVELS)].copy()

FRACTION_TICKS = FRACTION_LEVELS.astype(float)

FRACTION_TICK_LABELS = [
    "0", "0.025", "0.05", "0.075", "0.1",
    "0.15", "0.2", "0.25", "0.3", "0.35",
    "0.4", "0.45", "0.5",
    "0.6", "0.7", "0.8", "0.9", "1",
]

FRACTION_RGB = LAI_RGB[1:].copy()
FRACTION_RGB[0] = _as_rgb([[210, 230, 255]])[0]

# IDL Z0 levels used by compute_zo for ascat/icarus/merged.  Keep
# labels as strings so matplotlib cannot round the sub-1 bins to repeated
# ``0`` labels on the horizontal colorbar.
Z0_LEVELS = np.asarray([
    0.02, 0.05, 0.07, 0.10, 0.30, 0.50,
    1.0, 2.0, 4.0, 6.0, 8.0, 10.0,
    50.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0,
], dtype=np.float32)
Z0_TICK_LABELS = [
    "≤0.02", "0.05", "0.07", "0.10", "0.30", "0.50",
    "1", "2", "4", "6", "8", "10",
    "50", "100", "500", "1000", "2000", "3000", "4000", "5000",
]
Z0_COLOR_IDS = [74, 77, 35, 34, 33, 32, 25, 24, 23, 22, 21, 20, 41, 42, 43, 44, 45, 46, 47, 48]

# User-tunable output quality for static JPG products. These are set from
# --dpi and --jpeg-quality in main().
PLOT_DPI = int(os.environ.get("CLSM_PLOT_DPI", "180"))
JPEG_QUALITY = int(os.environ.get("CLSM_JPEG_QUALITY", "95"))


# -----------------------------------------------------------------------------
# File readers
# -----------------------------------------------------------------------------


class ClsmPlotError(RuntimeError):
    pass


@dataclasses.dataclass(frozen=True)
class F77Layout:
    endian: str = "<"
    marker_bytes: int = 4

    @property
    def marker_dtype(self) -> np.dtype:
        if self.marker_bytes == 4:
            return np.dtype(self.endian + "i4")
        if self.marker_bytes == 8:
            return np.dtype(self.endian + "i8")
        raise ValueError("record markers must be 4 or 8 bytes")




@dataclasses.dataclass(frozen=True)
class TimeSeriesLayout:
    """Layout for LAI/GREEN/NDVI/AlbMap-style seasonal time-series files.
    
    Most make_bcs files are F77 sequential records with a header record
    containing at least the 9 IDL date fields, followed by an ncat-value data
    record. Finished BCS trees may expose renamed/symlinked files with extended
    headers or double-precision header fields. A few test copies may be raw
    header+data streams without F77 record markers.
    """    
    mode: str = "f77"
    endian: str = "<"
    marker_bytes: int = 4
    header_dtype: str = "f4"
    value_dtype: str = "f4"

    @property
    def f77_layout(self) -> F77Layout:
        return F77Layout(self.endian, self.marker_bytes)


class TimeSeriesReader:
    def __init__(self, path: Path, layout: TimeSeriesLayout, ncat: int):
        self.path = Path(path)
        self.layout = layout
        self.ncat = int(ncat)
        self.rdr = None
        self.fh = None

    def __enter__(self) -> "TimeSeriesReader":
        if self.layout.mode == "f77":
            self.rdr = FortranSequentialReader(self.path, self.layout.f77_layout)
        else:
            self.fh = self.path.open("rb")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[override]
        if self.rdr is not None:
            self.rdr.close()
        if self.fh is not None:
            self.fh.close()

    def _dt(self, kind: str) -> np.dtype:
        return np.dtype(self.layout.endian + kind)

    def read_record(self) -> Tuple[np.ndarray, np.ndarray]:
        hdt = self._dt(self.layout.header_dtype)
        vdt = self._dt(self.layout.value_dtype)
        if self.layout.mode == "f77":
            if self.rdr is None:
                raise ClsmPlotError("TimeSeriesReader not opened")
            hpayload = self.rdr.read_record_bytes()
            header = np.frombuffer(hpayload, dtype=hdt)
            if header.size < 9:
                raise ClsmPlotError(f"header record in {self.path} has {header.size} values, expected at least 9")
            vpayload = self.rdr.read_record_bytes()
            values = np.frombuffer(vpayload, dtype=vdt)
            if values.size != self.ncat:
                raise ClsmPlotError(f"data record in {self.path} has {values.size} values, expected {self.ncat}")
            # IDL only reads the first 9 values from the header record.  Extra
            # fields may be present in finalized lai_clim/green/ndvi files.
            return header[:9].astype(np.float64).copy(), values.astype(np.float32).copy()
        if self.fh is None:
            raise ClsmPlotError("TimeSeriesReader not opened")
        hbytes = self.fh.read(9 * hdt.itemsize)
        if not hbytes:
            raise EOFError(f"end of file in {self.path}")
        if len(hbytes) != 9 * hdt.itemsize:
            raise ClsmPlotError(f"short raw header in {self.path}")
        vbytes = self.fh.read(self.ncat * vdt.itemsize)
        if len(vbytes) != self.ncat * vdt.itemsize:
            raise ClsmPlotError(f"short raw data record in {self.path}")
        return np.frombuffer(hbytes, dtype=hdt).astype(np.float64).copy(), np.frombuffer(vbytes, dtype=vdt).astype(np.float32).copy()


def detect_timeseries_layout(path: Path, ncat: int) -> TimeSeriesLayout:
    path = Path(path)
    head = path.read_bytes()[:128]
    # F77 sequential: IDL reads the first 9 floats from each header record:
    #   readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    # In finished BCS trees the header record can contain more than those 9
    # values.  For example lai_clim_* often starts with a 56-byte record, i.e.
    # 14 float32 values: the first 9 are the dates and the remaining values are
    # metadata.  Therefore accept any float32/float64 header record with at least
    # 9 values, then verify that the following data record has ncat values.
    for marker_bytes in (4, 8):
        for endian in ("<", ">"):
            if len(head) < marker_bytes:
                continue
            mdtype = np.dtype(endian + ("i4" if marker_bytes == 4 else "i8"))
            nhead = int(np.frombuffer(head[:marker_bytes], dtype=mdtype)[0])
            if nhead <= 0:
                continue
            for hkind, hsize in (("f4", 4), ("f8", 8)):
                if nhead % hsize != 0 or nhead < 9 * hsize:
                    continue
                try:
                    with FortranSequentialReader(path, F77Layout(endian, marker_bytes)) as rdr:
                        hpayload = rdr.read_record_bytes()
                        header = np.frombuffer(hpayload, dtype=np.dtype(endian + hkind))
                        if header.size < 9:
                            continue
                        # Simple sanity check for the date fields read by IDL.
                        # year offsets are usually 0/1/2, months 1..12, days 1..31.
                        if not (0 <= float(header[1]) <= 12 and 0 <= float(header[7]) <= 12):
                            continue
                        vpayload = rdr.read_record_bytes()
                    if len(vpayload) == int(ncat) * 4:
                        print(f"Detected time-series layout for {path.name}: mode=f77 endian={endian} marker={marker_bytes} header_values={header.size} header_dtype={hkind} value_dtype=f4")
                        return TimeSeriesLayout("f77", endian, marker_bytes, hkind, "f4")
                    if len(vpayload) == int(ncat) * 8:
                        print(f"Detected time-series layout for {path.name}: mode=f77 endian={endian} marker={marker_bytes} header_values={header.size} header_dtype={hkind} value_dtype=f8")
                        return TimeSeriesLayout("f77", endian, marker_bytes, hkind, "f8")
                except Exception:
                    pass
    # Raw fallback: header/data/header/data with no F77 markers.  Accept both
    # the exact IDL 9-value header and the 14-value extended header.
    size = path.stat().st_size
    for nheader in (9, 14):
        rec4 = (nheader + int(ncat)) * 4
        rec8 = (nheader + int(ncat)) * 8
        if rec4 > 0 and size % rec4 == 0:
            print(f"Detected time-series layout for {path.name}: mode=raw header_values={nheader} dtype=f4")
            return TimeSeriesLayout("raw", "<", 0, "f4", "f4")
        if rec8 > 0 and size % rec8 == 0:
            print(f"Detected time-series layout for {path.name}: mode=raw header_values={nheader} dtype=f8")
            return TimeSeriesLayout("raw", "<", 0, "f8", "f8")
    raise ClsmPlotError(
        f"Could not detect LAI/GREEN/NDVI time-series layout for {path}. "
        "Expected an F77 header record with at least 9 float values followed by an ncat-value record, "
        "or a raw stream of header+ncat float values."
    )


def choose_timeseries_layout(path: Path, ncat: int, endian: str = "auto", marker_bytes: int = 0) -> TimeSeriesLayout:
    # Auto is safest because finished-layout symlinks can expose f4 or f8 headers.
    if endian == "auto" or marker_bytes == 0:
        return detect_timeseries_layout(path, ncat)
    endian_char = "<" if endian in ("little", "<") else ">"
    # When endian and record-marker are explicitly specified, assume a
    # single-precision F77 layout. Use auto detection for files that may
    # have extended or double-precision headers.    
    return TimeSeriesLayout("f77", endian_char, marker_bytes, "f4", "f4")

class FortranSequentialReader:
    """Read F77 sequential unformatted records."""

    def __init__(self, path: Path, layout: F77Layout):
        self.path = Path(path)
        self.layout = layout
        self.fh = self.path.open("rb")

    def close(self) -> None:
        self.fh.close()

    def __enter__(self) -> "FortranSequentialReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[override]
        self.close()

    def read_record_bytes(self) -> bytes:
        m = self.layout.marker_bytes
        start = self.fh.read(m)
        if not start:
            raise EOFError(f"end of file in {self.path}")
        if len(start) != m:
            raise ClsmPlotError(f"short F77 record marker in {self.path}")
        nbytes = int(np.frombuffer(start, dtype=self.layout.marker_dtype)[0])
        if nbytes < 0:
            raise ClsmPlotError(f"negative F77 record length {nbytes} in {self.path}")
        payload = self.fh.read(nbytes)
        if len(payload) != nbytes:
            raise ClsmPlotError(f"short F77 record payload in {self.path}: wanted {nbytes}, got {len(payload)}")
        end = self.fh.read(m)
        if len(end) != m:
            raise ClsmPlotError(f"short trailing F77 record marker in {self.path}")
        end_nbytes = int(np.frombuffer(end, dtype=self.layout.marker_dtype)[0])
        if end_nbytes != nbytes:
            raise ClsmPlotError(
                f"F77 marker mismatch in {self.path}: leading={nbytes}, trailing={end_nbytes}"
            )
        return payload

    def read_array(self, dtype: np.dtype | str, count: Optional[int] = None) -> np.ndarray:
        payload = self.read_record_bytes()
        dt = np.dtype(dtype).newbyteorder(self.layout.endian)
        arr = np.frombuffer(payload, dtype=dt)
        if count is not None and arr.size != count:
            raise ClsmPlotError(
                f"record in {self.path} has {arr.size} values of {dt}, expected {count}"
            )
        return arr.copy()


def detect_f77_layout(path: Path, expected_payload_bytes: int) -> F77Layout:
    head = Path(path).read_bytes()[:32]
    for marker_bytes in (4, 8):
        for endian in ("<", ">"):
            if len(head) < marker_bytes:
                continue
            dtype = np.dtype(endian + ("i4" if marker_bytes == 4 else "i8"))
            n = int(np.frombuffer(head[:marker_bytes], dtype=dtype)[0])
            if n == expected_payload_bytes:
                return F77Layout(endian=endian, marker_bytes=marker_bytes)
    raise ClsmPlotError(
        f"Could not detect F77 layout for {path}. Expected first record payload "
        f"{expected_payload_bytes} bytes. Try --endian little|big and/or --record-marker."
    )


def choose_layout(path: Path, expected_payload_bytes: int, endian: str, marker_bytes: int) -> F77Layout:
    if endian == "auto" or marker_bytes == 0:
        return detect_f77_layout(path, expected_payload_bytes)
    endian_char = "<" if endian in ("little", "<") else ">"
    return F77Layout(endian=endian_char, marker_bytes=marker_bytes)


def is_netcdf(path: Path) -> bool:
    try:
        magic = Path(path).read_bytes()[:8]
    except FileNotFoundError:
        return False
    return magic.startswith(b"CDF") or magic.startswith(b"\x89HDF\r\n\x1a\n")


def load_ascii_table(path: Path, skiprows: int = 0, min_cols: Optional[int] = None) -> np.ndarray:
    if not path.exists():
        raise ClsmPlotError(f"Missing required file: {path}")
    arr = np.loadtxt(path, comments="#", skiprows=skiprows)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if min_cols is not None and arr.shape[1] < min_cols:
        raise ClsmPlotError(f"{path} has {arr.shape[1]} columns; expected at least {min_cols}")
    return arr


def read_ncat(base_dir: Path) -> int:
    path = base_dir / "catchment.def"
    if not path.exists():
        raise ClsmPlotError(f"Missing catchment definition: {path}")
    with path.open("r") as fh:
        first = fh.readline().split()
    if not first:
        raise ClsmPlotError(f"Empty catchment definition: {path}")
    return int(float(first[0]))


def read_limits(base_dir: Path, gfile: str) -> Tuple[float, float, float, float]:
    """Return IDL-style default map limits as (min_lat, min_lon, max_lat, max_lon)."""
    default = (-60.0, -180.0, 90.0, 180.0)
    if "Pfafstetter" in gfile or "SMAP" in gfile:
        return default
    path = base_dir / "catchment.def"
    try:
        rows = load_ascii_table(path, skiprows=1, min_cols=6)
    except Exception:
        return default
    min_lon = np.nanmin(rows[:, 2])
    max_lon = np.nanmax(rows[:, 3])
    min_lat = np.nanmin(rows[:, 4])
    max_lat = np.nanmax(rows[:, 5])
    if math.ceil(max_lon) - math.floor(min_lon) < 180:
        return (math.floor(min_lat), math.floor(min_lon), math.ceil(max_lat), math.ceil(max_lon))
    return default



def list_rst_files(workdir: Path, gfile: str, explicit_rst_file: Optional[str] = None) -> List[Path]:
    """Return candidate rst files.

    IDL used a wildcard, ``<workdir>/rst/<gfile>*.rst``.  For EASE grids this
    can include several companion rasters, for example the pure EASE raster and
    the EASE/Pfafstetter combined raster.  The combined raster is often the one
    that maps pixels to the CLSM/catchment tile IDs used by ``catchment.def``.

    Therefore do *not* discard filenames containing ``Pfafstetter`` here.  Keep
    every ``<gfile>*.rst`` candidate and let ``select_rst_file`` score which one
    actually behaves like a geographic CLSM tile-id raster.  The user can still
    force any file with ``--rst-file``.
    """
    if explicit_rst_file:
        p = Path(explicit_rst_file).expanduser()
        if not p.is_absolute():
            p = Path.cwd() / p
        if not p.exists():
            raise ClsmPlotError(f"Explicit --rst-file does not exist: {p}")
        return [p.resolve()]

    rst_dir = workdir / "rst"
    pattern = str(rst_dir / f"{gfile}*.rst")
    matches = [Path(m).resolve() for m in sorted(glob.glob(pattern))]
    if not matches:
        raise ClsmPlotError(f"No restart/raster file matched {pattern}")

    # Stable unique order.
    seen = set()
    out: List[Path] = []
    for m in matches:
        if m not in seen:
            seen.add(m)
            out.append(m)
    return out


@dataclasses.dataclass(frozen=True)
class RstCandidateScore:
    path: Path
    layout: F77Layout
    valid_fraction: float
    valid_count: int
    sample_count: int
    min_value: int
    max_value: int
    first_valid_row: int
    last_valid_row: int
    valid_row_count: int
    sampled_row_count: int

    @property
    def row_span_fraction(self) -> float:
        if self.first_valid_row < 0 or self.last_valid_row < 0 or self.sampled_row_count <= 1:
            return 0.0
        return float(self.last_valid_row - self.first_valid_row) / float(max(1, self.sampled_row_count - 1))


def sample_rst_candidate(
    path: Path,
    nc: int,
    nr: int,
    ncat: int,
    endian: str,
    marker_bytes: int,
    sample_rows: int = 144,
) -> RstCandidateScore:
    """Sample an rst candidate and summarize whether it looks like CLSM tile IDs.

    A misleading raster can contain many values in 1..ncat but only over a narrow
    projected band when interpreted as lon/lat.  The original IDL maps are
    geographic/global, so for global grids the correct raster should have valid
    land rows spread across much of the sampled row range.  We therefore record
    both the number of valid tile IDs and the sampled row span.
    """
    layout = choose_layout(path, nc * 4, endian, marker_bytes)
    record_bytes = layout.marker_bytes + nc * 4 + layout.marker_bytes
    row_ids = np.unique(np.linspace(0, nr - 1, min(sample_rows, nr), dtype=np.int64))
    valid_count = 0
    sample_count = 0
    min_val: Optional[int] = None
    max_val: Optional[int] = None
    first_valid_idx = -1
    last_valid_idx = -1
    valid_row_count = 0
    marker_dtype = layout.marker_dtype
    data_dtype = np.dtype(layout.endian + "i4")

    with path.open("rb") as fh:
        for sample_idx, row in enumerate(row_ids):
            fh.seek(int(row) * record_bytes)
            marker = fh.read(layout.marker_bytes)
            if len(marker) != layout.marker_bytes:
                raise ClsmPlotError(f"short marker while sampling {path} row {row}")
            nbytes = int(np.frombuffer(marker, dtype=marker_dtype)[0])
            if nbytes != nc * 4:
                raise ClsmPlotError(
                    f"record length mismatch while sampling {path} row {row}: "
                    f"got {nbytes}, expected {nc * 4}"
                )
            payload = fh.read(nc * 4)
            if len(payload) != nc * 4:
                raise ClsmPlotError(f"short payload while sampling {path} row {row}")
            arr = np.frombuffer(payload, dtype=data_dtype)
            sample_count += int(arr.size)
            valid = (arr >= 1) & (arr <= ncat)
            row_valid = int(np.count_nonzero(valid))
            valid_count += row_valid
            if row_valid > 0:
                if first_valid_idx < 0:
                    first_valid_idx = sample_idx
                last_valid_idx = sample_idx
                valid_row_count += 1
            row_min = int(arr.min())
            row_max = int(arr.max())
            min_val = row_min if min_val is None else min(min_val, row_min)
            max_val = row_max if max_val is None else max(max_val, row_max)

    frac = float(valid_count) / float(sample_count) if sample_count else 0.0
    return RstCandidateScore(
        path=path,
        layout=layout,
        valid_fraction=frac,
        valid_count=valid_count,
        sample_count=sample_count,
        min_value=int(min_val or 0),
        max_value=int(max_val or 0),
        first_valid_row=first_valid_idx,
        last_valid_row=last_valid_idx,
        valid_row_count=valid_row_count,
        sampled_row_count=int(len(row_ids)),
    )


def select_rst_file(
    workdir: Path,
    gfile: str,
    ncat: int,
    nc: int,
    nr: int,
    endian: str,
    marker_bytes: int,
    explicit_rst_file: Optional[str] = None,
) -> Tuple[Path, F77Layout]:
    """Choose the raster that actually looks like the CLSM tile-id raster."""
    candidates = list_rst_files(workdir, gfile, explicit_rst_file)
    if len(candidates) > 1:
        print("RST candidates after filtering:")
    scored: List[Tuple[float, float, float, int, RstCandidateScore]] = []
    errors: List[str] = []
    for idx, path in enumerate(candidates):
        try:
            score = sample_rst_candidate(path, nc, nr, ncat, endian, marker_bytes)
            # Primary: row span.  Secondary: number/fraction of valid tile IDs.
            # Keep original order as the final tie breaker.
            scored.append((score.row_span_fraction, score.valid_fraction, float(score.valid_count), -idx, score))
            if len(candidates) > 1 or explicit_rst_file:
                print(
                    f"  {path.name}: valid_sample={score.valid_count}/{score.sample_count} "
                    f"({score.valid_fraction:.6f}), valid_rows={score.valid_row_count}/{score.sampled_row_count}, "
                    f"row_span={score.row_span_fraction:.3f}, min={score.min_value}, max={score.max_value}, "
                    f"endian={score.layout.endian}, marker={score.layout.marker_bytes}"
                )
        except Exception as exc:
            errors.append(f"{path.name}: {exc}")
            if len(candidates) > 1 or explicit_rst_file:
                print(f"  {path.name}: skipped ({exc})")

    if not scored:
        detail = "\n".join(errors)
        raise ClsmPlotError(f"No usable rst files found under {workdir / 'rst'}.\n{detail}")

    # Prefer candidates with the broadest sampled latitude/row coverage.  This
    # avoids selecting projection/Pfafstetter companion rasters that have valid
    # numeric ranges but do not make global geographic plots.
    scored.sort(key=lambda item: (-item[0], -item[1], -item[2], -item[3]))
    best = scored[0][4]
    if best.valid_count == 0:
        raise ClsmPlotError(
            f"Selected rst candidate has zero values in 1..ncat: {best.path}. "
            "Check --gfile/--workdir/--nc/--nr or pass --rst-file explicitly."
        )
    if len(candidates) > 1:
        print(
            f"Selected raster: {best.path} "
            f"(row_span={best.row_span_fraction:.3f}, valid_sample_fraction={best.valid_fraction:.6f})"
        )
    elif explicit_rst_file:
        print(
            f"Using explicit raster: {best.path} "
            f"(row_span={best.row_span_fraction:.3f}, valid_sample_fraction={best.valid_fraction:.6f})"
        )
    else:
        print(f"Using raster file: {best.path}")
    return best.path, best.layout

def read_nc_var(path: Path, varname: str) -> np.ndarray:
    if xr is None:
        raise ClsmPlotError("xarray/netCDF support is not available in this Python environment")
    if not path.exists():
        raise ClsmPlotError(f"Missing NetCDF file: {path}")
    with xr.open_dataset(path, decode_times=False) as ds:
        if varname not in ds:
            raise ClsmPlotError(f"Variable {varname!r} not found in {path}. Available: {list(ds.data_vars)}")
        return np.asarray(ds[varname].values)


# -----------------------------------------------------------------------------
# Grid and tile tools
# -----------------------------------------------------------------------------


def lon_lat_centers(nc: int, nr: int) -> Tuple[np.ndarray, np.ndarray]:
    lon = np.arange(nc, dtype=np.float64) * (360.0 / nc) - 180.0 + 0.5 * (360.0 / nc)
    lat = np.arange(nr, dtype=np.float64) * (180.0 / nr) - 90.0 + 0.5 * (180.0 / nr)
    return lon, lat


def mode_rows_int(blocks: np.ndarray) -> np.ndarray:
    """Return the row-wise mode for integer blocks."""
    if blocks.ndim != 2:
        raise ValueError("blocks must be 2D")
    if blocks.shape[1] == 1:
        return blocks[:, 0].astype(np.int32, copy=False)
    if scipy_mode is not None:
        result = scipy_mode(blocks, axis=1, keepdims=False)
        return np.asarray(result.mode, dtype=np.int32)
    out = np.empty(blocks.shape[0], dtype=np.int32)
    for i, row in enumerate(blocks):
        vals, counts = np.unique(row, return_counts=True)
        out[i] = vals[np.argmax(counts)]
    return out


def dominant_land_tile(blocks: np.ndarray, ncat: int) -> np.ndarray:
    """IDL-compatible dominant tile selection for a downsampled raster block.

    The IDL code ignores ocean/invalid values when a plotting cell contains at
    least one land tile.  A naive mode over values with invalid pixels replaced
    by zero would let ocean dominate coastal cells, so we first compute a fast
    mode and then repair only cells where zero won but valid land is present.
    """
    valid = (blocks >= 1) & (blocks <= ncat)
    masked = np.where(valid, blocks, 0).astype(np.int32, copy=False)
    out = mode_rows_int(masked)

    repair = (out == 0) & valid.any(axis=1)
    if np.any(repair):
        repair_idx = np.flatnonzero(repair)
        for i in repair_idx:
            vals, counts = np.unique(blocks[i, valid[i]], return_counts=True)
            out[i] = vals[np.argmax(counts)]
    return out


def build_tile_id_from_rst(
    rst_file: Path,
    nc: int,
    nr: int,
    ncat: int,
    nc_plot: int,
    nr_plot: int,
    layout: F77Layout,
    cache: Optional[Path] = None,
) -> np.ndarray:
    """Build the dominant catchment tile map.  Returns array shape (nr_plot, nc_plot)."""
    if cache and cache.exists():
        data = np.load(cache, allow_pickle=True)
        meta = dict(data["meta"].item()) if "meta" in data else {}
        rst_stat = rst_file.stat()
        cache_ok = (
            meta.get("nc") == nc
            and meta.get("nr") == nr
            and meta.get("ncat") == ncat
            and meta.get("nc_plot") == nc_plot
            and meta.get("nr_plot") == nr_plot
            and meta.get("rst_file") == str(rst_file.resolve())
            and meta.get("rst_size") == int(rst_stat.st_size)
        )
        if cache_ok:
            print(f"Reading cached tile map: {cache}")
            return np.asarray(data["tile_id"], dtype=np.int32)
        print(f"Ignoring stale cache: {cache}")

    if nc % nc_plot != 0 or nr % nr_plot != 0:
        raise ClsmPlotError(
            f"NC/NR must be integer multiples of plot grid. Got NC={nc}, NR={nr}, "
            f"plot={nc_plot}x{nr_plot}. Try --plot-nc/--plot-nr."
        )
    dx = nc // nc_plot
    dy = nr // nr_plot
    print(f"Building tile map from {rst_file}: source {nc}x{nr}, plot {nc_plot}x{nr_plot}, block {dx}x{dy}")
    tile_id = np.zeros((nr_plot, nc_plot), dtype=np.int32)

    with FortranSequentialReader(rst_file, layout) as rdr:
        for j in range(nr_plot):
            rows = np.empty((dy, nc), dtype=np.int32)
            for jj in range(dy):
                rows[jj, :] = rdr.read_array(np.int32, nc)
            if dx == 1 and dy == 1:
                vals = rows[0].copy()
                vals[(vals < 1) | (vals > ncat)] = 0
                tile_id[j, :] = vals
            else:
                block = rows.reshape(dy, nc_plot, dx).transpose(1, 0, 2).reshape(nc_plot, dx * dy)
                tile_id[j, :] = dominant_land_tile(block, ncat)
            if (j + 1) % max(1, nr_plot // 20) == 0 or j == nr_plot - 1:
                print(f"  tile map rows {j + 1}/{nr_plot}")

    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            cache,
            tile_id=tile_id,
            meta={
                "nc": nc,
                "nr": nr,
                "ncat": ncat,
                "nc_plot": nc_plot,
                "nr_plot": nr_plot,
                "rst_file": str(rst_file.resolve()),
                "rst_size": int(rst_file.stat().st_size),
            },
        )
        print(f"Wrote tile-map cache: {cache}")
    return tile_id


def vector_to_grid(tile_id: np.ndarray, vec: np.ndarray, fill_zero: bool = False) -> np.ndarray:
    vals = np.asarray(vec)
    out = np.full(tile_id.shape, np.nan, dtype=np.float32)
    mask = (tile_id >= 1) & (tile_id <= vals.shape[0])
    out[mask] = vals[tile_id[mask] - 1].astype(np.float32)
    if fill_zero:
        out[~mask] = 0.0
    return out

def build_tile_id_from_catchment_def(
    base_dir: Path,
    ncat: int,
    nc_plot: int,
    nr_plot: int,
    cache: Optional[Path] = None,
) -> np.ndarray:
    """Build a plotting tile map directly from catchment.def lon/lat boxes.

    This is a fallback for grids such as EASE where ``<gfile>*.rst``
    companion rasters can contain projection-cell or Pfafstetter IDs rather
    than the CLSM tile IDs that index the tile-parameter files.  IDL primarily
    used the rst map, but all of the global fixed-parameter maps ultimately
    need only a lon/lat plotting grid whose values are row numbers in
    ``catchment.def``/``cti_stats.dat``/``soil_param.dat``.  The catchment
    bounds in columns 3:6 of ``catchment.def`` provide that mapping.
    """
    cdef = base_dir / "catchment.def"
    if cache and cache.exists():
        data = np.load(cache, allow_pickle=True)
        meta = dict(data["meta"].item()) if "meta" in data else {}
        st = cdef.stat()
        if (
            meta.get("source") == "catchment.def"
            and meta.get("ncat") == ncat
            and meta.get("nc_plot") == nc_plot
            and meta.get("nr_plot") == nr_plot
            and meta.get("catchment_def") == str(cdef.resolve())
            and meta.get("catchment_size") == int(st.st_size)
        ):
            print(f"Reading cached catchment.def tile map: {cache}")
            return np.asarray(data["tile_id"], dtype=np.int32)
        print(f"Ignoring stale catchment.def cache: {cache}")

    rows = load_ascii_table(cdef, skiprows=1, min_cols=6)
    if rows.shape[0] < ncat:
        raise ClsmPlotError(f"{cdef} has {rows.shape[0]} rows but ncat={ncat}")
    rows = rows[:ncat]
    minlon = rows[:, 2].astype(float)
    maxlon = rows[:, 3].astype(float)
    minlat = rows[:, 4].astype(float)
    maxlat = rows[:, 5].astype(float)

    dx = 360.0 / float(nc_plot)
    dy = 180.0 / float(nr_plot)
    tile_id = np.zeros((nr_plot, nc_plot), dtype=np.int32)

    def lat_slice(lo: float, hi: float) -> Optional[slice]:
        if not np.isfinite(lo) or not np.isfinite(hi):
            return None
        lo = max(-90.0, min(90.0, lo))
        hi = max(-90.0, min(90.0, hi))
        if hi < lo:
            lo, hi = hi, lo
        # Fill cells whose area intersects the catchment box.  This guarantees
        # very small boxes still get at least one plot pixel.
        j0 = int(math.floor((lo + 90.0) / dy))
        j1 = int(math.ceil((hi + 90.0) / dy)) - 1
        j0 = max(0, min(nr_plot - 1, j0))
        j1 = max(0, min(nr_plot - 1, j1))
        if j1 < j0:
            j = max(0, min(nr_plot - 1, int(round(((lo + hi) * 0.5 + 90.0) / dy - 0.5))))
            j0 = j1 = j
        return slice(j0, j1 + 1)

    def lon_slices(lo: float, hi: float) -> List[slice]:
        if not np.isfinite(lo) or not np.isfinite(hi):
            return []
        # Normalize to [-180, 180).  If the original box spans the dateline,
        # split into two slices.
        lo0, hi0 = lo, hi
        lo = ((lo + 180.0) % 360.0) - 180.0
        hi = ((hi + 180.0) % 360.0) - 180.0
        wraps = (lo0 > hi0) or (lo > hi and abs(lo - hi) < 359.999)

        def one_slice(a: float, b: float) -> Optional[slice]:
            a = max(-180.0, min(180.0, a))
            b = max(-180.0, min(180.0, b))
            i0 = int(math.floor((a + 180.0) / dx))
            i1 = int(math.ceil((b + 180.0) / dx)) - 1
            i0 = max(0, min(nc_plot - 1, i0))
            i1 = max(0, min(nc_plot - 1, i1))
            if i1 < i0:
                i = max(0, min(nc_plot - 1, int(round(((a + b) * 0.5 + 180.0) / dx - 0.5))))
                i0 = i1 = i
            return slice(i0, i1 + 1)

        if not wraps:
            sl = one_slice(lo, hi)
            return [sl] if sl is not None else []
        out: List[slice] = []
        sl1 = one_slice(lo, 180.0)
        sl2 = one_slice(-180.0, hi)
        if sl1 is not None:
            out.append(sl1)
        if sl2 is not None:
            out.append(sl2)
        return out

    print(f"Building tile map from catchment.def boxes: plot {nc_plot}x{nr_plot}")
    for k in range(ncat):
        js = lat_slice(float(minlat[k]), float(maxlat[k]))
        if js is None:
            continue
        for is_ in lon_slices(float(minlon[k]), float(maxlon[k])):
            tile_id[js, is_] = k + 1
        if (k + 1) % max(1, ncat // 10) == 0 or k == ncat - 1:
            print(f"  catchment boxes {k + 1}/{ncat}")

    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            cache,
            tile_id=tile_id,
            meta={
                "source": "catchment.def",
                "ncat": ncat,
                "nc_plot": nc_plot,
                "nr_plot": nr_plot,
                "catchment_def": str(cdef.resolve()),
                "catchment_size": int(cdef.stat().st_size),
            },
        )
        print(f"Wrote catchment.def tile-map cache: {cache}")
    return tile_id


def catchment_spatial_match_fraction(tile_id: np.ndarray, base_dir: Path, lon: np.ndarray, lat: np.ndarray, ncat: int, max_samples: int = 200000) -> float:
    """Return fraction of sampled tile_id cells whose lon/lat falls in its catchment.def box."""
    valid = (tile_id >= 1) & (tile_id <= ncat)
    rows, cols = np.where(valid)
    if rows.size == 0:
        return 0.0
    if rows.size > max_samples:
        step = int(math.ceil(rows.size / float(max_samples)))
        rows = rows[::step]
        cols = cols[::step]
    ids = tile_id[rows, cols].astype(np.int64) - 1
    cdef = load_ascii_table(base_dir / "catchment.def", skiprows=1, min_cols=6)[:ncat]
    minlon = cdef[:, 2].astype(float)
    maxlon = cdef[:, 3].astype(float)
    minlat = cdef[:, 4].astype(float)
    maxlat = cdef[:, 5].astype(float)
    x = lon[cols]
    y = lat[rows]
    tol_lon = 360.0 / float(lon.size) + 1e-6
    tol_lat = 180.0 / float(lat.size) + 1e-6
    inlat = (y >= minlat[ids] - tol_lat) & (y <= maxlat[ids] + tol_lat)
    normal = minlon[ids] <= maxlon[ids]
    inlon_normal = (x >= minlon[ids] - tol_lon) & (x <= maxlon[ids] + tol_lon)
    inlon_wrap = (x >= minlon[ids] - tol_lon) | (x <= maxlon[ids] + tol_lon)
    inlon = np.where(normal, inlon_normal, inlon_wrap)
    return float(np.count_nonzero(inlat & inlon)) / float(rows.size)


def build_fractional_sparse_from_rst(
    rst_file: Path,
    nc: int,
    nr: int,
    ncat: int,
    nc_out: int,
    nr_out: int,
    layout: F77Layout,
    cache: Optional[Path] = None,
):
    if sp_sparse is None:
        raise ClsmPlotError("scipy.sparse is required for fractional movie aggregation")
    if cache and cache.exists():
        print(f"Reading cached fractional mapping: {cache}")
        return sp_sparse.load_npz(cache)
    if nc % nc_out != 0 or nr % nr_out != 0:
        raise ClsmPlotError(
            f"NC/NR must be integer multiples of movie grid. Got NC={nc}, NR={nr}, movie={nc_out}x{nr_out}."
        )
    dx = nc // nc_out
    dy = nr // nr_out
    row_idx: List[int] = []
    col_idx: List[int] = []
    weight: List[float] = []
    print(f"Building fractional mapping from {rst_file}: movie {nc_out}x{nr_out}, block {dx}x{dy}")
    with FortranSequentialReader(rst_file, layout) as rdr:
        for j in range(nr_out):
            rows = np.empty((dy, nc), dtype=np.int32)
            for jj in range(dy):
                rows[jj, :] = rdr.read_array(np.int32, nc)
            for i in range(nc_out):
                subset = rows[:, i * dx:(i + 1) * dx].ravel()
                valid = subset[(subset >= 1) & (subset <= ncat)]
                if valid.size == 0:
                    continue
                ids, counts = np.unique(valid, return_counts=True)
                cell = j * nc_out + i
                row_idx.extend([cell] * len(ids))
                col_idx.extend((ids - 1).astype(int).tolist())
                # IDL divides by the full block count, not only land count.
                weight.extend((counts / float(subset.size)).astype(float).tolist())
            if (j + 1) % max(1, nr_out // 20) == 0 or j == nr_out - 1:
                print(f"  fractional rows {j + 1}/{nr_out}")
    mat = sp_sparse.csr_matrix((weight, (row_idx, col_idx)), shape=(nc_out * nr_out, ncat), dtype=np.float32)
    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        sp_sparse.save_npz(cache, mat)
        print(f"Wrote fractional mapping cache: {cache}")
    return mat


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------


def _limits_to_extent(limits: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    min_lat, min_lon, max_lat, max_lon = limits
    return (min_lon, max_lon, min_lat, max_lat)


def crop_to_limits(grid: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits: Tuple[float, float, float, float]):
    min_lat, min_lon, max_lat, max_lon = limits
    lon_mask = (lon >= min_lon) & (lon <= max_lon)
    lat_mask = (lat >= min_lat) & (lat <= max_lat)
    if not lon_mask.any() or not lat_mask.any():
        return grid, lon, lat, _limits_to_extent(limits)
    sub = grid[np.ix_(lat_mask, lon_mask)]
    lon_sub = lon[lon_mask]
    lat_sub = lat[lat_mask]
    dx = 360.0 / lon.size
    dy = 180.0 / lat.size
    extent = (lon_sub[0] - dx / 2, lon_sub[-1] + dx / 2, lat_sub[0] - dy / 2, lat_sub[-1] + dy / 2)
    return sub, lon_sub, lat_sub, extent


def centers_to_edges(a: np.ndarray) -> np.ndarray:
    """Convert a 1-D regular center-coordinate array to edges."""
    a = np.asarray(a, dtype=float)
    if a.size == 0:
        return a
    if a.size == 1:
        # Fall back to a one-degree cell for pathological one-cell debug plots.
        return np.asarray([a[0] - 0.5, a[0] + 0.5], dtype=float)
    d = float(np.median(np.diff(a)))
    return np.r_[a[0] - 0.5 * d, a + 0.5 * d]


def add_land_outline(ax, sub: np.ndarray, lon_sub: np.ndarray, lat_sub: np.ndarray, coastlines: bool, linewidth: float = 0.75) -> None:
    """Draw a no-dependency land/ocean outline from the valid-data mask.

    IDL drew MAP_CONTINENTS on every panel.  When Cartopy is unavailable or
    --no-coastlines is used, the plots otherwise lose all continent outlines.
    This mask-derived outline is not a political boundary dataset, but it gives
    the same visual coast/land edge cue and works on Discover without extra
    data downloads.
    """
    if coastlines and ccrs is not None:
        return
    try:
        good = np.isfinite(sub).astype(float)
        if good.shape[0] < 2 or good.shape[1] < 2 or np.nanmax(good) <= 0:
            return
        ax.contour(lon_sub, lat_sub, good, levels=[0.5], colors="black", linewidths=linewidth)
    except Exception:
        pass


def add_horizontal_category_key(
    fig: plt.Figure,
    colors: np.ndarray,
    labels: Sequence[str],
    *,
    x0: float = 0.16,
    y0: float = 0.035,
    width: float = 0.68,
    height: float = 0.035,
    label_size: float = 7.0,
    title: str = "",
    tick_rotation: float = 45.0,
) -> None:
    """Add an IDL-style horizontal categorical color strip."""
    n = len(labels)
    if n == 0:
        return
    key_ax = fig.add_axes([x0, y0, width, height])
    cmap = ListedColormap(np.asarray(colors[:n], dtype=np.float32))
    key_ax.imshow(np.arange(n, dtype=float)[None, :], cmap=cmap, aspect="auto", extent=[0, n, 0, 1], interpolation="nearest")
    key_ax.set_yticks([])
    key_ax.set_xticks(np.arange(n) + 0.5)
    key_ax.set_xticklabels(labels, rotation=tick_rotation, fontsize=label_size, ha="right", va="top")
    key_ax.tick_params(axis="x", length=0, pad=3)
    for spine in key_ax.spines.values():
        spine.set_linewidth(0.8)
    if title:
        key_ax.set_xlabel(title, fontsize=8, labelpad=4)


def read_rst_subset(
    rst_file: Path,
    nc: int,
    nr: int,
    layout: F77Layout,
    limits: Tuple[float, float, float, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a lon/lat subset from a fixed-record F77 int32 raster.

    This is used for US-east.jpg.  The IDL code reads the full-resolution rst
    raster directly for that regional plot; using the globally downsampled map
    creates blocky catchments and apparent leakage.
    """
    min_lat, min_lon, max_lat, max_lon = limits
    dx = 360.0 / float(nc)
    dy = 180.0 / float(nr)
    i1 = max(0, int(math.floor((min_lon + 180.0) / dx)))
    i2 = min(nc - 1, int(math.ceil((max_lon + 180.0) / dx)) - 1)
    j1 = max(0, int(math.floor((min_lat + 90.0) / dy)))
    j2 = min(nr - 1, int(math.ceil((max_lat + 90.0) / dy)) - 1)
    if i2 < i1 or j2 < j1:
        raise ClsmPlotError(f"Invalid rst subset for limits={limits}")
    width = i2 - i1 + 1
    height = j2 - j1 + 1
    dtype = np.dtype(layout.endian + "i4")
    rec_bytes = layout.marker_bytes + nc * 4 + layout.marker_bytes
    out = np.empty((height, width), dtype=np.int32)
    with Path(rst_file).open("rb") as fh:
        for jj, j in enumerate(range(j1, j2 + 1)):
            fh.seek(j * rec_bytes + layout.marker_bytes + i1 * 4)
            row = np.fromfile(fh, dtype=dtype, count=width)
            if row.size != width:
                raise ClsmPlotError(f"Short read from {rst_file} row {j}: got {row.size}, wanted {width}")
            out[jj, :] = row.astype(np.int32, copy=False)
    lon_sub = -180.0 + (np.arange(i1, i2 + 1) + 0.5) * dx
    lat_sub = -90.0 + (np.arange(j1, j2 + 1) + 0.5) * dy
    return out, lon_sub, lat_sub


def build_boundary_segments_from_tile_ids(
    tile_ids: np.ndarray,
    lon_sub: np.ndarray,
    lat_sub: np.ndarray,
    ncat: int,
) -> Tuple[List[Tuple[Tuple[float, float], Tuple[float, float]]], List[Tuple[Tuple[float, float], Tuple[float, float]]]]:
    """Return vertical and horizontal catchment/coast boundary line segments.

    IDL's plot_tiles draws an oplot line at every pixel edge where the raster
    category changes, including category-to-ocean edges.  Baking one-pixel black
    edges into a full-resolution image can vanish when the figure is resampled,
    so draw true matplotlib line segments in lon/lat coordinates instead.
    """
    valid = (tile_ids >= 1) & (tile_ids <= int(ncat))
    if tile_ids.size == 0:
        return [], []
    lon_edges = centers_to_edges(lon_sub)
    lat_edges = centers_to_edges(lat_sub)
    v_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    h_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    # Vertical boundaries between neighboring columns.
    diff_v = tile_ids[:, 1:] != tile_ids[:, :-1]
    draw_v = diff_v & (valid[:, 1:] | valid[:, :-1])
    rows, cols = np.where(draw_v)
    for r, c in zip(rows.tolist(), cols.tolist()):
        x = float(lon_edges[c + 1])
        v_segments.append(((x, float(lat_edges[r])), (x, float(lat_edges[r + 1]))))
    # Left/right outside edges where a valid cell borders the regional/ocean edge.
    for r in range(tile_ids.shape[0]):
        if valid[r, 0]:
            x = float(lon_edges[0])
            v_segments.append(((x, float(lat_edges[r])), (x, float(lat_edges[r + 1]))))
        if valid[r, -1]:
            x = float(lon_edges[-1])
            v_segments.append(((x, float(lat_edges[r])), (x, float(lat_edges[r + 1]))))
    # Horizontal boundaries between neighboring rows.
    diff_h = tile_ids[1:, :] != tile_ids[:-1, :]
    draw_h = diff_h & (valid[1:, :] | valid[:-1, :])
    rows, cols = np.where(draw_h)
    for r, c in zip(rows.tolist(), cols.tolist()):
        y = float(lat_edges[r + 1])
        h_segments.append(((float(lon_edges[c]), y), (float(lon_edges[c + 1]), y)))
    # Bottom/top outside edges.
    for c in range(tile_ids.shape[1]):
        if valid[0, c]:
            y = float(lat_edges[0])
            h_segments.append(((float(lon_edges[c]), y), (float(lon_edges[c + 1]), y)))
        if valid[-1, c]:
            y = float(lat_edges[-1])
            h_segments.append(((float(lon_edges[c]), y), (float(lon_edges[c + 1]), y)))
    return v_segments, h_segments


def add_tile_boundary_lines(ax, tile_ids: np.ndarray, lon_sub: np.ndarray, lat_sub: np.ndarray, ncat: int, linewidth: float = 0.23) -> None:
    v_segments, h_segments = build_boundary_segments_from_tile_ids(tile_ids, lon_sub, lat_sub, ncat)
    if v_segments:
        ax.add_collection(LineCollection(v_segments, colors="black", linewidths=linewidth, antialiaseds=False, zorder=5))
    if h_segments:
        ax.add_collection(LineCollection(h_segments, colors="black", linewidths=linewidth, antialiaseds=False, zorder=5))


def make_axes(fig: plt.Figure, nrows: int, ncols: int, idx: int, coastlines: bool):
    if coastlines and ccrs is not None:
        ax = fig.add_subplot(nrows, ncols, idx, projection=ccrs.PlateCarree())
    else:
        ax = fig.add_subplot(nrows, ncols, idx)
    return ax


def _nice_geo_ticks(lo: float, hi: float, is_lon: bool) -> np.ndarray:
    """Return readable lon/lat tick locations for IDL-like map axes."""
    span = float(hi) - float(lo)
    if span >= 300.0:
        step = 60.0
    elif span >= 150.0:
        step = 30.0
    elif span >= 70.0:
        step = 15.0
    elif span >= 25.0:
        step = 5.0
    elif span >= 10.0:
        step = 2.0
    else:
        step = 1.0
    start = math.ceil(float(lo) / step) * step
    stop = math.floor(float(hi) / step) * step
    ticks = np.arange(start, stop + 0.5 * step, step, dtype=float)
    if ticks.size == 0:
        ticks = np.asarray([lo, hi], dtype=float)
    # Keep the full-domain endpoints when they are part of the requested map.
    if is_lon:
        if lo <= -179.999 and not np.isclose(ticks[0], -180.0):
            ticks = np.r_[-180.0, ticks]
        if hi >= 179.999 and not np.isclose(ticks[-1], 180.0):
            ticks = np.r_[ticks, 180.0]
    else:
        if lo <= -89.999 and not np.isclose(ticks[0], -90.0):
            ticks = np.r_[-90.0, ticks]
        if hi >= 89.999 and not np.isclose(ticks[-1], 90.0):
            ticks = np.r_[ticks, 90.0]
    # Avoid too many labels in narrow panels.
    if ticks.size > 9:
        ticks = ticks[:: int(math.ceil(ticks.size / 9.0))]
    return ticks


def _plain_lon_label(x: float) -> str:
    x = float(x)
    if abs(x) < 1e-9:
        return "0°"
    hemi = "E" if x > 0 else "W"
    return f"{abs(x):g}°{hemi}"


def _plain_lat_label(y: float) -> str:
    y = float(y)
    if abs(y) < 1e-9:
        return "0°"
    hemi = "N" if y > 0 else "S"
    return f"{abs(y):g}°{hemi}"


def decorate_geo(
    ax,
    limits: Tuple[float, float, float, float],
    coastlines: bool,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
) -> None:
    min_lat, min_lon, max_lat, max_lon = limits
    if coastlines and ccrs is not None and hasattr(ax, "set_extent"):
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
        ax.coastlines(linewidth=0.5)
        try:
            ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        except Exception:
            pass

        # Cartopy axes do not show normal Matplotlib lon/lat ticks unless we
        # explicitly set them.  Without this, movie frames with --coastlines
        # have coastlines but lose the longitude/latitude labels.
        try:
            xticks = _nice_geo_ticks(min_lon, max_lon, is_lon=True)
            yticks = _nice_geo_ticks(min_lat, max_lat, is_lon=False)
            ax.set_xticks(xticks, crs=ccrs.PlateCarree())
            ax.set_yticks(yticks, crs=ccrs.PlateCarree())
            if LongitudeFormatter is not None:
                ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=False))
            else:
                ax.set_xticklabels([_plain_lon_label(x) for x in xticks])
            if LatitudeFormatter is not None:
                ax.yaxis.set_major_formatter(LatitudeFormatter())
            else:
                ax.set_yticklabels([_plain_lat_label(y) for y in yticks])
            ax.tick_params(
                labelsize=7,
                bottom=show_xlabel, labelbottom=show_xlabel,
                left=show_ylabel, labelleft=show_ylabel,
                top=False, right=False,
                pad=2,
            )
            ax.set_xlabel("Longitude" if show_xlabel else "", fontsize=8, labelpad=6)
            ax.set_ylabel("Latitude" if show_ylabel else "", fontsize=8, labelpad=6)
            try:
                ax.gridlines(
                    xlocs=xticks, ylocs=yticks, draw_labels=False,
                    linewidth=0.25, color="0.35", alpha=0.35, linestyle="-",
                )
            except Exception:
                pass
        except Exception:
            # Coastlines are more important than labels; avoid failing plots if a
            # particular Cartopy build cannot format projected tick labels.
            pass
    else:
        ax.set_xlim(min_lon, max_lon)
        ax.set_ylim(min_lat, max_lat)
        ax.set_xlabel("Longitude" if show_xlabel else "", labelpad=6)
        ax.set_ylabel("Latitude" if show_ylabel else "", labelpad=6)
        ax.grid(True, linewidth=0.2, alpha=0.4)
        try:
            ax.set_aspect("equal", adjustable="box")
        except Exception:
            pass

def plot_continuous_on_ax(
    ax,
    grid: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    limits: Tuple[float, float, float, float],
    title: str,
    levels: Sequence[float],
    color_ids: Optional[Sequence[int]] = None,
    rgb: Optional[np.ndarray] = None,
    coastlines: bool = False,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
    bad_color = "white",
):
    sub, lon_sub, lat_sub, extent = crop_to_limits(grid, lon, lat, limits)
    finite = np.isfinite(sub)
    if finite.any():
        print(f"  plot {title}: finite={int(finite.sum())}/{sub.size}, min={float(np.nanmin(sub)):.6g}, max={float(np.nanmax(sub)):.6g}")
    else:
        print(f"  plot {title}: finite=0/{sub.size} -- output will be blank")
    if rgb is not None:
        cmap = ListedColormap(rgb)
    else:
        cmap = ListedColormap(PALETTE[np.asarray(color_ids or CONTINUOUS_COLOR_IDS)])
    levels_arr = np.asarray(levels, dtype=float)
    # BoundaryNorm expects one more boundary than colors.  IDL uses each listed
    # level as a filled-contour break; extend one upper boundary if needed.
    if len(levels_arr) == cmap.N:
        step = levels_arr[-1] - levels_arr[-2] if len(levels_arr) > 1 else 1.0
        boundaries = np.r_[levels_arr, levels_arr[-1] + step]
    else:
        boundaries = levels_arr
    norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    # Convert data to an explicit RGBA image before putting it
    # on the axes.  On Discover some Agg/pcolormesh combinations were producing
    # blank-looking panels even though the arrays contained valid data.  This
    # follows the successful direct-image debug path.
    cmap.set_bad(bad_color)
    rgba = cmap(norm(np.ma.masked_invalid(sub)))
    imshow_kwargs = dict(
        origin="lower",
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    if coastlines and ccrs is not None:
        imshow_kwargs["transform"] = ccrs.PlateCarree()
    ax.imshow(rgba, **imshow_kwargs)
    add_land_outline(ax, sub, lon_sub, lat_sub, coastlines)
    try:
        ax.set_aspect("equal", adjustable="box")
    except Exception:
        pass
    ax.set_title(title, fontsize=10, pad=10)
    decorate_geo(ax, limits, coastlines, show_xlabel=show_xlabel, show_ylabel=show_ylabel)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sm


def plot_indexed_on_ax(
    ax,
    grid: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    limits: Tuple[float, float, float, float],
    title: str,
    classes: Sequence[int],
    colors: Sequence[Sequence[float]] | np.ndarray,
    labels: Optional[Sequence[str]] = None,
    coastlines: bool = False,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
):
    sub, lon_sub, lat_sub, extent = crop_to_limits(grid, lon, lat, limits)
    finite = np.isfinite(sub)
    if finite.any():
        vals = np.unique(sub[finite])
        print(f"  plot {title}: finite={int(finite.sum())}/{sub.size}, unique={vals.size}, first={vals[:8]}")
    else:
        print(f"  plot {title}: finite=0/{sub.size} -- output will be blank")
    cls = np.asarray(classes)
    idx_grid = np.full(sub.shape, np.nan, dtype=np.float32)
    for k, value in enumerate(cls):
        idx_grid[sub == value] = k
    cmap = ListedColormap(np.asarray(colors, dtype=np.float32))
    norm = BoundaryNorm(np.arange(-0.5, len(cls) + 0.5, 1.0), cmap.N)
    cmap.set_bad("white")
    rgba = cmap(norm(np.ma.masked_invalid(idx_grid)))
    imshow_kwargs = dict(
        origin="lower",
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    if coastlines and ccrs is not None:
        imshow_kwargs["transform"] = ccrs.PlateCarree()
    ax.imshow(rgba, **imshow_kwargs)
    add_land_outline(ax, sub, lon_sub, lat_sub, coastlines)
    try:
        ax.set_aspect("equal", adjustable="box")
    except Exception:
        pass
    ax.set_title(title, fontsize=10, pad=10)
    decorate_geo(ax, limits, coastlines, show_xlabel=show_xlabel, show_ylabel=show_ylabel)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    if labels:
        cbar = plt.colorbar(sm, ax=ax, shrink=0.75, pad=0.02, ticks=np.arange(len(cls)))
        cbar.ax.set_yticklabels(labels)
    return sm


def save_fig(fig: plt.Figure, outpath: Path, dpi: Optional[int] = None) -> None:
    """Save a static plot with package-wide DPI/JPEG quality settings."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    use_dpi = int(PLOT_DPI if dpi is None else dpi)
    kwargs = dict(dpi=use_dpi, bbox_inches="tight", facecolor="white")
    if outpath.suffix.lower() in (".jpg", ".jpeg"):
        kwargs["pil_kwargs"] = {"quality": int(JPEG_QUALITY), "optimize": True}
    try:
        fig.savefig(outpath, **kwargs)
    except TypeError:
        # Older Matplotlib builds may not support pil_kwargs.
        kwargs.pop("pil_kwargs", None)
        fig.savefig(outpath, **kwargs)
    plt.close(fig)
    print(f"Wrote {outpath}")


def panel_continuous(
    grids: Sequence[np.ndarray],
    titles: Sequence[str],
    lon: np.ndarray,
    lat: np.ndarray,
    limits: Tuple[float, float, float, float],
    outpath: Path,
    ncols: int,
    levels_list: Sequence[Sequence[float]],
    color_ids_list: Optional[Sequence[Sequence[int]]] = None,
    rgb_list: Optional[Sequence[np.ndarray]] = None,
    figsize: Tuple[float, float] = (10, 8),
    coastlines: bool = False,
    bad_color = "white",
) -> None:
    n = len(grids)
    nrows = int(math.ceil(n / ncols))
    # Multi-panel products need more vertical real estate than the default
    # Matplotlib layout, especially once each panel has its own colorbar.
    # Make this global rather than only fixing cti.jpg.  This prevents titles
    # such as POROS/COND/T2m from colliding with the axis labels of the panel
    # above.
    compact_two_column = (ncols == 2 and n in (4, 6))
    if compact_two_column:
        # 4-panel (SoilAlb) and 6-panel (soil_param) global maps were
        # visually too tall because Cartopy/geographic aspect makes each map
        # panel wide and shallow.  Use explicit compact figure heights for
        # these layouts instead of the generic multi-panel height rule.
        if n == 4:
            figsize = (max(figsize[0], 12.0), 5.6)
        else:  # n == 6
            figsize = (max(figsize[0], 12.2), 7.35)
    elif nrows > 1:
        # Keep enough room for titles/colorbars, but avoid very large gaps.
        min_h_per_row = 3.35 if ncols <= 2 else 2.85
        figsize = (max(figsize[0], 9.5 if ncols == 1 else figsize[0]), max(figsize[1], min_h_per_row * nrows))
    fig = plt.figure(figsize=figsize)
    if compact_two_column:
        fig.subplots_adjust(hspace=0.12, wspace=0.18, top=0.955, bottom=0.105, left=0.070, right=0.985)
    elif nrows > 1:
        hspace = 0.46 if ncols <= 2 else 0.34
        wspace = 0.20 if ncols > 1 else 0.14
        fig.subplots_adjust(hspace=hspace, wspace=wspace, top=0.965, bottom=0.085, left=0.075, right=0.975)
    last_im = None
    for k, grid in enumerate(grids):
        ax = make_axes(fig, nrows, ncols, k + 1, coastlines)
        row = k // ncols
        col = k % ncols
        show_xlabel = row == (nrows - 1)
        show_ylabel = col == 0
        last_im = plot_continuous_on_ax(
            ax,
            grid,
            lon,
            lat,
            limits,
            titles[k],
            levels=levels_list[k],
            color_ids=(color_ids_list[k] if color_ids_list is not None else None),
            rgb=(rgb_list[k] if rgb_list is not None else None),
            coastlines=coastlines,
            show_xlabel=show_xlabel,
            show_ylabel=show_ylabel,
            bad_color=bad_color,
        )
        cbar = fig.colorbar(last_im, ax=ax, shrink=0.65, pad=0.02)
        cbar.ax.tick_params(labelsize=7)
    save_fig(fig, outpath)


# -----------------------------------------------------------------------------
# Plot products translated from clsm_plots.pro
# -----------------------------------------------------------------------------

def plot_tiles(
    tile_id: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    limits: Tuple[float, float, float, float],
    outdir: Path,
    coastlines: bool,
    rst_file: Optional[Path] = None,
    nc_full: Optional[int] = None,
    nr_full: Optional[int] = None,
    ncat: Optional[int] = None,
    layout: Optional[F77Layout] = None,
    gfile: str = "",
) -> None:
    n_levels = 30
    colors = PALETTE[np.arange(90, 90 + n_levels)]

    if nc_full is not None and nr_full is not None:
        raster_res = max(int(nc_full), int(nr_full))
    else:
        raster_res = 0

    rst_name = str(rst_file).upper() if rst_file is not None else ""
    grid_name = f"{gfile} {rst_name}".upper()

    is_ease = "EASE" in grid_name
    is_ease_m01 = is_ease and "M01" in grid_name
    is_ease_m03 = is_ease and "M03" in grid_name

    # Cubed-sphere logical resolution.
    # Examples:
    #   CF0180x6C_DE1440xPE0720
    #   CF2880x6C_CF2880x6C
    #   CF2160x6C-SG001_CF2160x6C
    m_cf = re.search(r"CF0*([0-9]+)X6C", grid_name)
    cf_res = int(m_cf.group(1)) if m_cf else 0

    # Lat-lon / data-ocean style names.
    # Examples:
    #   DC0288xPC0181_DE0360xPE0180
    #   DE1440xPE0720
    is_latlon = bool(
        re.search(r"(?:DC|DE)0*[0-9]+X(?:PC|PE)0*[0-9]+", grid_name)
    )

    if is_ease_m01:
        # EASE 1-km: tight zoom.
        us_east_limits = (38.35, -76.45, 38.75, -75.95)

    elif is_ease_m03:
        # EASE 3-km: moderate zoom.
        us_east_limits = (38.0, -77.2, 39.2, -75.4)

    elif is_ease:
        # EASE M09/M25/M36 and any other EASE not explicitly zoomed.
        us_east_limits = (35.0, -82.0, 42.0, -73.0)

    elif cf_res > 0 and cf_res <= 720:
        # C12 through C720: broad region, like original IDL.
        us_east_limits = (35.0, -82.0, 42.0, -73.0)

    elif cf_res >= 5760:
        # C5760: tight Chesapeake zoom.
        us_east_limits = (38.35, -76.45, 38.75, -75.95)

    elif cf_res >= 2160:
        # C2160/C2880/C3072 and fine stretched.
        us_east_limits = (38.0, -76.8, 39.0, -75.6)

    elif cf_res >= 768:
        # C768/C1000/C1080/C1120/C1152/C1440/C1536.
        us_east_limits = (37.6, -77.2, 39.2, -75.4)

    elif is_latlon:
        # Lat-lon b/c/d/e grids: broad region.
        us_east_limits = (35.0, -82.0, 42.0, -73.0)

    else:
        # Preserve IDL-like behavior for unrecognized grids/res.
        # Do NOT fall back to raster_res-based zoom here.
        us_east_limits = (35.0, -82.0, 42.0, -73.0)

    print(
        f"  plot US-east catchment tile zoom: "
        f"gfile={gfile}, cf_res={cf_res}, is_latlon={is_latlon}, "
        f"raster_res={raster_res}, limits={us_east_limits}"
    )
    if rst_file is not None and nc_full is not None and nr_full is not None and ncat is not None and layout is not None:
        sub_id, lon_sub, lat_sub = read_rst_subset(Path(rst_file), int(nc_full), int(nr_full), layout, us_east_limits)
        valid = (sub_id >= 1) & (sub_id <= int(ncat))
        idx = np.full(sub_id.shape, np.nan, dtype=np.float32)
        idx[valid] = (sub_id[valid] % n_levels).astype(np.float32)
        cmap = ListedColormap(colors)
        cmap.set_bad("white")
        rgba = cmap(np.ma.masked_invalid(idx.astype(float) / max(1, n_levels - 1)))
        dx = 360.0 / float(nc_full)
        dy = 180.0 / float(nr_full)
        extent = (lon_sub[0] - dx / 2, lon_sub[-1] + dx / 2, lat_sub[0] - dy / 2, lat_sub[-1] + dy / 2)
        fig = plt.figure(figsize=(7.0, 5.0))
        ax = make_axes(fig, 1, 1, 1, coastlines)
        imshow_kwargs = dict(origin="lower", extent=extent, interpolation="nearest", aspect="equal", zorder=1)
        if coastlines and ccrs is not None:
            imshow_kwargs["transform"] = ccrs.PlateCarree()
        ax.imshow(rgba, **imshow_kwargs)

        # Boundary linewidth is grid/resolution dependent.
        # Coarse grids need visible borders so same-color neighboring tiles
        # are still separable. Fine zooms need thinner borders.
        if is_ease_m01:
            tile_boundary_linewidth = 0.10
        elif is_ease_m03:
            tile_boundary_linewidth = 0.16
        elif is_ease:
            tile_boundary_linewidth = 0.23

        elif cf_res > 0 and cf_res <= 720:
            tile_boundary_linewidth = 0.23
        elif cf_res >= 5760:
            tile_boundary_linewidth = 0.08
        elif cf_res >= 2160:
            tile_boundary_linewidth = 0.12
        elif cf_res >= 768:
            tile_boundary_linewidth = 0.16

        elif is_latlon:
            tile_boundary_linewidth = 0.23

        else:
            tile_boundary_linewidth = 0.23

        if tile_boundary_linewidth > 0.0:
            add_tile_boundary_lines(
                ax, sub_id, lon_sub, lat_sub, int(ncat),
                linewidth=tile_boundary_linewidth
            )

        decorate_geo(ax, us_east_limits, coastlines)
        ax.grid(False)
        ax.set_title("Catchment tiles", fontsize=10, pad=8)
        save_fig(fig, outdir / "US-east.jpg")
        return

    # Fallback for unusual runs where the rst file is unavailable.
    grid = np.full(tile_id.shape, np.nan, dtype=np.float32)
    mask = tile_id > 0
    grid[mask] = tile_id[mask] % n_levels
    fig = plt.figure(figsize=(7.0, 5.0))
    ax = make_axes(fig, 1, 1, 1, coastlines)
    plot_indexed_on_ax(ax, grid, lon, lat, us_east_limits, "Catchment tiles", np.arange(n_levels), colors, coastlines=coastlines)
    save_fig(fig, outdir / "US-east.jpg")

def plot_country_codes(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "country_and_state_code.data"
    if not path.exists():
        print(f"Skipping country codes; missing {path}")
        return
    # File has numeric columns followed by text labels such as UNK; read only
    # the first three numeric columns instead of np.loadtxt-ing the whole row.
    rows = np.genfromtxt(path, comments="#", usecols=(0, 1, 2), dtype=np.float64, invalid_raise=False)
    if rows.ndim == 1:
        rows = rows.reshape(1, -1)
    if rows.shape[1] < 3:
        raise ClsmPlotError(f"{path} has fewer than three numeric columns")
    ncat = int(np.nanmax(tile_id))
    if rows.shape[0] < ncat:
        print(f"WARNING: {path} has {rows.shape[0]} rows but tile map references {ncat} tiles")
    cnt = rows[:ncat, 1].astype(np.float32)
    st = rows[:ncat, 2].astype(np.float32)
    us = cnt == 243
    cnt[us] = st[us]
    cnt[cnt == 257] = np.nan
    grid = vector_to_grid(tile_id, cnt + 1)
    # Use a reproducible random-looking palette for up to 256 codes.
    rng = np.random.default_rng(12345)
    colors = rng.random((256, 3))
    colors[0] = 0
    colors[-1] = 1
    fig = plt.figure(figsize=(10, 5))
    ax = make_axes(fig, 1, 1, 1, coastlines)
    sub, lon_sub, lat_sub, extent = crop_to_limits(grid, lon, lat, limits)
    finite = np.isfinite(sub)
    if finite.any():
        vals = np.unique(sub[finite])
        print(f"  plot Country / state codes: finite={int(finite.sum())}/{sub.size}, unique={vals.size}, first={vals[:8]}")
    else:
        print(f"  plot Country / state codes: finite=0/{sub.size} -- output will be blank")
    cmap = ListedColormap(colors)
    cmap.set_bad("white")
    # Wrap arbitrary numeric country/state codes into the available palette for
    # a stable categorical image.  The exact colors do not need to encode the
    # numeric magnitude.
    idx = np.full(sub.shape, np.nan, dtype=np.float32)
    good = np.isfinite(sub)
    idx[good] = (sub[good].astype(np.int64) % colors.shape[0]).astype(np.float32)
    norm = BoundaryNorm(np.arange(-0.5, colors.shape[0] + 0.5, 1.0), cmap.N)
    rgba = cmap(norm(np.ma.masked_invalid(idx)))
    ax.imshow(rgba, origin="lower", extent=extent, interpolation="nearest", aspect="equal")
    try:
        ax.set_aspect("equal", adjustable="box")
    except Exception:
        pass
    ax.set_title("Country / state codes")
    decorate_geo(ax, limits, coastlines)
    save_fig(fig, outdir / "Country_codes.jpg")


def plot_cti(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "cti_stats.dat"
    rows = load_ascii_table(path, skiprows=1, min_cols=7)
    cti_mean = rows[:, 2].astype(np.float32)
    cti_std = rows[:, 3].astype(np.float32)
    cti_skew = rows[:, 6].astype(np.float32)
    if not (base_dir / "CLM_veg_typs_fracs").exists():
        cti_mean = 0.961 * cti_mean - 1.957
    grids = [vector_to_grid(tile_id, v) for v in (cti_mean, cti_std, cti_skew)]
    levels = [np.linspace(6.0, 14.0, 17), np.linspace(0.0, 4.0, 17), np.linspace(-2.5, 2.5, 17)]
    panel_continuous(
        grids, ["CTI mean", "CTI std", "CTI skew"], lon, lat, limits, outdir / "cti.jpg",
        ncols=1, levels_list=levels, figsize=(10.4, 10.1), coastlines=coastlines,
    )


def plot_mosaic(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "mosaic_veg_typs_fracs"
    rows = load_ascii_table(path, min_cols=3)
    mos_type = rows[:, 2].astype(int)
    grid = vector_to_grid(tile_id, mos_type)
    vtypes = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 14, 20, 30, 40, 50, 60, 70, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230]
    r = [233, 255, 255, 255, 210, 0, 0, 0, 204, 170, 255, 220, 205, 0, 0, 170, 0, 40, 120, 140, 190, 150, 255, 255, 0, 0, 0, 195, 255, 0, 255, 0]
    g = [23, 131, 191, 255, 255, 255, 155, 0, 204, 240, 255, 240, 205, 100, 160, 200, 60, 100, 130, 160, 150, 100, 180, 235, 120, 150, 220, 20, 245, 70, 255, 0]
    b = [0, 0, 0, 178, 255, 255, 255, 200, 204, 240, 100, 100, 102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 175, 90, 120, 130, 0, 215, 200, 255, 0]
    colors = _as_rgb(list(zip(r, g, b)))
    fig = plt.figure(figsize=(10.8, 5.4))
    fig.subplots_adjust(bottom=0.20, top=0.93)
    ax = make_axes(fig, 1, 1, 1, coastlines)
    plot_indexed_on_ax(ax, grid, lon, lat, limits, "Mosaic primary vegetation type", vtypes, colors, labels=None, coastlines=coastlines)
    # The IDL map uses the full vtypes color table, but the legend intentionally
    # labels only the six broad mosaic classes.
    mos_labels = ["BL Evergreen", "BL Deciduous", "Needleleaf", "Grassland", "BL Shrubs", "Dwarf"]
    add_horizontal_category_key(fig, colors[:6], mos_labels, x0=0.33, y0=0.050, width=0.36, height=0.026, label_size=7, tick_rotation=45.0)
    save_fig(fig, outdir / "mosaic_prim.jpg")

def _read_clm_rows(base_dir: Path) -> Optional[np.ndarray]:
    path = base_dir / "CLM_veg_typs_fracs"
    if not path.exists():
        print(f"Skipping CLM/Catchment-CN vegetation plots; missing {path}")
        return None
    return load_ascii_table(path, min_cols=12)


def plot_clm(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    rows = _read_clm_rows(base_dir)
    if rows is None:
        return
    clm = rows[:, [10, 11]].astype(int)
    colors = _as_rgb(list(zip(
        [255,106,202,251,0,29,77,109,142,233,255,255,127,164,217,204,0],
        [245,91,178,154,85,115,145,165,185,23,131,191,39,53,72,204,70],
        [215,154,214,153,0,0,0,0,13,0,0,0,4,3,1,204,200],
    )))
    classes = list(range(1, 18))
    labels = ["BARE", "NLEt", "NLEB", "NLDB", "BLET", "BLEt", "BLDT", "BLDt", "BLDB", "BLEtS", "BLDtS", "BLDBS", "AC3G", "CC3G", "WC4G", "CROP"]
    for idx, name in enumerate(["PRIM", "SEC"]):
        grid = vector_to_grid(tile_id, clm[:, idx])
        fig = plt.figure(figsize=(10, 6))
        ax = make_axes(fig, 1, 1, 1, coastlines)
        fig.subplots_adjust(bottom=0.18, top=0.92)
        plot_indexed_on_ax(ax, grid, lon, lat, limits, f"CLM {name} vegetation type", classes, colors, labels=None, coastlines=coastlines)
        add_horizontal_category_key(fig, colors[:len(labels)], labels, x0=0.16, y0=0.045, width=0.68, height=0.026, label_size=7, tick_rotation=45.0)
        save_fig(fig, outdir / f"CLM_{name}_veg_typs.jpg")


def plot_carbon(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    rows = _read_clm_rows(base_dir)
    if rows is None:
        return
    cn = rows[:, [2, 3, 4, 5]].astype(int)
    colors = _as_rgb(list(zip(
        [106,202,251,0,29,77,109,142,233,255,255,255,127,164,164,217,217,204,104,0],
        [91,178,154,85,115,145,165,185,23,131,131,191,39,53,53,72,72,204,104,70],
        [154,214,153,0,0,0,0,13,0,0,200,0,4,3,200,1,200,204,200,200],
    )))
    classes = list(range(1, 21))
    labels = ["NLEt", "NLEB", "NLDB", "BLET", "BLEt", "BLDT", "BLDt", "BLDB", "BLEtS", "BLDtS", "BLDtSm", "BLDBS", "AC3G", "CC3G", "CC3Gm", "WC4G", "WC4Gm", "CROP", "CROPm"]
    for label, cols in [("PRIM", [0, 1]), ("SEC", [2, 3])]:
        fig = plt.figure(figsize=(12.0, 8.6))
        fig.subplots_adjust(hspace=0.16, bottom=0.17, top=0.955, left=0.06, right=0.98)
        for k, c in enumerate(cols):
            ax = make_axes(fig, 2, 1, k + 1, coastlines)
            grid = vector_to_grid(tile_id, cn[:, c])
            plot_indexed_on_ax(ax, grid, lon, lat, limits, f"Catchment-CN {label} vegetation {k + 1}", classes, colors, labels=None, coastlines=coastlines, show_xlabel=(k == 1))
        add_horizontal_category_key(fig, colors[:len(labels)], labels, x0=0.16, y0=0.055, width=0.68, height=0.026, label_size=7, tick_rotation=45.0)
        save_fig(fig, outdir / f"CatchmentCN_{label}_veg_typs.jpg")


def plot_ndep_t2m_soilalb(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "CLM_NDep_SoilAlb_T2m"
    if not path.exists():
        print(f"Skipping NDep/T2m/SoilAlb; missing {path}")
        return
    rows = load_ascii_table(path, min_cols=7)
    ndep, visdr, visdf, nirdr, nirdf, t2mm, t2mp = [rows[:, i].astype(np.float32) for i in range(7)]
    grids = [vector_to_grid(tile_id, v) for v in (ndep, t2mm, t2mp)]
    levels = [np.asarray(list(np.arange(15) * 4.0) + [65.0, 350.0]), np.linspace(250.0, 300.0, 17), np.linspace(250.0, 300.0, 17)]
    panel_continuous(grids, ["NDep", "T2m mean", "T2m plus"], lon, lat, limits, outdir / "CLM_Ndep_T2m.jpg", ncols=1, levels_list=levels, figsize=(10.4, 9.8), coastlines=coastlines)
    soilalb_grids = [vector_to_grid(tile_id, v) for v in (visdr, visdf, nirdr, nirdf)]
    levels_alb = [np.linspace(0.0, 0.65, 17), np.linspace(0.0, 0.65, 17), np.linspace(0.0, 1.0, 17), np.linspace(0.0, 1.0, 17)]
    panel_continuous(soilalb_grids, ["VISDR", "VISDF", "NIRDR", "NIRDF"], lon, lat, limits, outdir / "SoilAlb.jpg", ncols=2, levels_list=levels_alb, figsize=(10, 7), coastlines=coastlines)


def plot_soil(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "soil_param.dat"
    rows = load_ascii_table(path, min_cols=10)
    vals = {
        "BEE": rows[:, 4].astype(np.float32),
        "PSIS": rows[:, 5].astype(np.float32),
        "POROS": rows[:, 6].astype(np.float32),
        "COND": rows[:, 7].astype(np.float32),
        "WPWET": rows[:, 8].astype(np.float32),
        "SOILDEPTH": rows[:, 9].astype(np.float32),
    }
    vlims = {
        "BEE": (1.0, 8.0),
        "PSIS": (-1.85, -0.1),
        "POROS": (0.37, 0.8),
        "COND": (2.37e-6, 2.845e-4),
        "WPWET": (0.01, 0.45),
        "SOILDEPTH": (1334.0, 5000.0),
    }
    grids, titles, levs = [], [], []
    for name in vals:
        lo, hi = vlims[name]
        if name == "POROS":
            levels = np.r_[lo, lo + np.arange(15) * ((0.57 - lo) / 15.0), hi]
        else:
            levels = np.linspace(lo, hi, 17)
        grids.append(vector_to_grid(tile_id, vals[name]))
        titles.append(name)
        levs.append(levels)
    panel_continuous(grids, titles, lon, lat, limits, outdir / "soil_param.jpg", ncols=2, levels_list=levs, figsize=(11, 9), coastlines=coastlines)


def plot_elevation(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, coastlines: bool) -> None:
    rows = load_ascii_table(base_dir / "catchment.def", skiprows=1, min_cols=7)
    elevation = rows[:, 6].astype(np.float32)
    grid = vector_to_grid(tile_id, elevation)
    lo = float(np.nanmin(elevation))
    hi = float(np.nanmax(elevation))
    panel_continuous([grid], ["ELEVATION"], lon, lat, limits, outdir / "ELEVATION.jpg", ncols=1, levels_list=[np.linspace(lo, hi, 17)], figsize=(11.2, 5.6), coastlines=coastlines)


# -----------------------------------------------------------------------------
# Time series LAI/GREEN/albedo readers and Z0
# -----------------------------------------------------------------------------


def read_timeseries_record(rdr, ncat: int) -> Tuple[np.ndarray, np.ndarray]:
    if isinstance(rdr, TimeSeriesReader):
        return rdr.read_record()
    header = rdr.read_array(np.float32, 9)
    values = rdr.read_array(np.float32, ncat)
    return header, values


def _doy_from_header_part(yroff: float, month: float, day: float) -> float:
    y = 2001 + int(round(float(yroff)))
    m = max(1, min(12, int(round(float(month)))))
    d = max(1, min(31, int(round(float(day)))))
    dt = _dt.date(y, m, d)
    return float((dt - _dt.date(2000, 12, 31)).days)


def midpoint_doy(header: np.ndarray) -> float:
    start = _doy_from_header_part(header[0], header[1], header[2])
    end = _doy_from_header_part(header[6], header[7], header[8])
    return (end - start) / 2.0 + start


def _sanitize_lai(values: np.ndarray, *, max_lai: float = 12.0) -> np.ndarray:
    """Return LAI with impossible/fill values masked as NaN.

    LAI in these climatology files should be non-negative and normally below
    about 8.  Use 12 as a conservative upper bound so real dense-canopy values
    are retained while bad extrapolated/fill values do not contaminate Z0.
    """
    out = np.asarray(values, dtype=np.float32).copy()
    bad = (~np.isfinite(out)) | (out < 0.0) | (out > max_lai) | (np.abs(out) > 1.0e10)
    out[bad] = np.nan
    return out


def _sanitize_ndvi(values: np.ndarray) -> np.ndarray:
    """Return NDVI with impossible/fill values masked as NaN."""
    out = np.asarray(values, dtype=np.float32).copy()
    bad = (~np.isfinite(out)) | (out < -0.1) | (out > 1.1) | (np.abs(out) > 1.0e10)
    out[bad] = np.nan
    # Keep a very small tolerance for files with roundoff; clip to physical range.
    out = np.where(np.isfinite(out), np.clip(out, 0.0, 1.0), np.nan).astype(np.float32)
    return out


def _sanitize_fraction(values: np.ndarray) -> np.ndarray:
    """Return generic fractional fields (GREEN/albedo) in 0..1 with fill values masked."""
    out = np.asarray(values, dtype=np.float32).copy()
    bad = (~np.isfinite(out)) | (out < -0.01) | (out > 1.01) | (np.abs(out) > 1.0e10)
    out[bad] = np.nan
    out = np.where(np.isfinite(out), np.clip(out, 0.0, 1.0), np.nan).astype(np.float32)
    return out


def _sanitize_z0_mm(values: np.ndarray, *, max_mm: float = 10000.0) -> np.ndarray:
    """Mask impossible roughness length values in millimeters."""
    out = np.asarray(values, dtype=np.float32).copy()
    bad = (~np.isfinite(out)) | (out < 0.0) | (out > max_mm) | (np.abs(out) > 1.0e20)
    out[bad] = np.nan
    return out


def _initial_loop_year_offset(h1: np.ndarray, h2: np.ndarray) -> int:
    """Mimic the IDL loop's active year variable after the second header read.

    In clsm_plots.pro, the variable ``yr`` is overwritten by the *second* header
    before the month/day loop starts.  For files whose first interval is late
    Dec -> Jan, using the first header's year makes January a full year too
    early and causes large negative extrapolated LAI.  Use h2[0] here to match
    IDL's state at loop entry.
    """
    return int(round(float(h2[0])))


def _advance_loop_year_offset(header: np.ndarray, month: int) -> int:
    yoff = int(round(float(header[0])))
    # IDL has: if((month eq 12) and (yr eq 2)) then yr = yr -1
    if month == 12 and yoff == 2:
        yoff -= 1
    return yoff


def _valid_weighted_monthly_add(total: np.ndarray, count: np.ndarray, month_index: int, vals: np.ndarray) -> None:
    good = np.isfinite(vals)
    if np.any(good):
        total[month_index, good] += vals[good].astype(np.float64)
        count[month_index, good] += 1.0


def monthly_means_from_interpolated(path: Path, ncat: int, layout: TimeSeriesLayout, kind: str = "raw") -> np.ndarray:
    mdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    total = np.zeros((12, ncat), dtype=np.float64)
    count = np.zeros((12, ncat), dtype=np.float64)
    if not path.exists():
        raise ClsmPlotError(f"Missing time series file: {path}")
    with TimeSeriesReader(path, layout, ncat) as rdr:
        h1, v1 = read_timeseries_record(rdr, ncat)
        h2, v2 = read_timeseries_record(rdr, ncat)
        b4 = midpoint_doy(h1)
        nxt = midpoint_doy(h2)
        current_year_offset = _initial_loop_year_offset(h1, h2)
        for month in range(1, 13):
            for day in range(1, mdays[month - 1] + 1):
                now = float((_dt.date(2001 + current_year_offset, month, day) - _dt.date(2000, 12, 31)).days)
                denom = (nxt - b4) if abs(nxt - b4) > 1e-6 else 1.0
                fac1 = (now - b4) / denom
                fac2 = (nxt - now) / denom
                vals = (fac1 * v2 + fac2 * v1).astype(np.float32)
                if kind.lower() == "lai":
                    vals = _sanitize_lai(vals)
                elif kind.lower() == "ndvi":
                    vals = _sanitize_ndvi(vals)
                elif kind.lower() in ("fraction", "green", "albedo"):
                    vals = _sanitize_fraction(vals)
                else:
                    vals = np.where(np.isfinite(vals), vals, np.nan).astype(np.float32)
                _valid_weighted_monthly_add(total, count, month - 1, vals)
                if now + 0.5 >= nxt:
                    v1 = v2
                    b4 = nxt
                    try:
                        h2, v2 = read_timeseries_record(rdr, ncat)
                        nxt = midpoint_doy(h2)
                        current_year_offset = _advance_loop_year_offset(h2, month)
                    except EOFError:
                        nxt = now + 9999.0
    monthly = np.full((12, ncat), np.nan, dtype=np.float32)
    good = count > 0.0
    monthly[good] = (total[good] / count[good]).astype(np.float32)
    return monthly


def panel_continuous_shared_colorbar(
    grids: Sequence[np.ndarray],
    titles: Sequence[str],
    lon: np.ndarray,
    lat: np.ndarray,
    limits: Tuple[float, float, float, float],
    outpath: Path,
    ncols: int,
    levels: Sequence[float],
    color_ids: Optional[Sequence[int]] = None,
    rgb: Optional[np.ndarray] = None,
    figsize: Tuple[float, float] = (13, 9),
    coastlines: bool = False,
    cbar_label: str = "",
    cbar_ticks: Optional[Sequence[float]] = None,
    cbar_ticklabels: Optional[Sequence[str]] = None,
    cbar_tick_rotation: float = 90.0,
    bad_color = "white",
) -> None:
    """Multi-panel plot with one shared colorbar.

    This is better for LAI and Z0 where all panels use identical bins.
    It avoids shrinking every map panel to make room for separate colorbars.
    """
    n = len(grids)
    nrows = int(math.ceil(n / ncols))
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.955, bottom=0.145, hspace=0.30, wspace=0.14)
    last_im = None
    for k, grid in enumerate(grids):
        ax = make_axes(fig, nrows, ncols, k + 1, coastlines)
        row = k // ncols
        col = k % ncols
        last_im = plot_continuous_on_ax(
            ax,
            grid,
            lon,
            lat,
            limits,
            titles[k],
            levels=levels,
            color_ids=color_ids,
            rgb=rgb,
            coastlines=coastlines,
            show_xlabel=(row == nrows - 1),
            show_ylabel=(col == 0),
            bad_color=bad_color,
        )
    if last_im is not None:
        cax = fig.add_axes([0.14, 0.060, 0.72, 0.024])
        if cbar_ticks is not None:
            tick_values = np.asarray(cbar_ticks, dtype=float)
            cbar = fig.colorbar(
                last_im,
                cax=cax,
                orientation="horizontal",
                ticks=tick_values,
                spacing="uniform",
            )
            cbar.ax.xaxis.set_major_locator(FixedLocator(tick_values))
            if cbar_ticklabels is not None:
                if len(cbar_ticklabels) != len(tick_values):
                    raise ClsmPlotError(
                        f"Colorbar label count {len(cbar_ticklabels)} does not match "
                        f"tick count {len(tick_values)} for {outpath}"
                    )
                cbar.ax.xaxis.set_major_formatter(FixedFormatter(list(cbar_ticklabels)))
        else:
            cbar = fig.colorbar(last_im, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=7, rotation=cbar_tick_rotation, pad=4)
        for label in cbar.ax.get_xticklabels():
            label.set_horizontalalignment("right" if abs(cbar_tick_rotation) > 1.0 else "center")
            label.set_verticalalignment("top")
        if cbar_label:
            cbar.set_label(cbar_label, fontsize=8, labelpad=7)
    save_fig(fig, outpath)


def plot_monthly_timeseries(
    base_dir: Path,
    tile_id: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    limits,
    outdir: Path,
    ncat: int,
    filename: str,
    outname: str,
    product_label: str,
    layout: Optional[TimeSeriesLayout],
    coastlines: bool,
    bad_color = NO_DATA_COLOR,
    *,
    kind: str = "raw",
    levels: Sequence[float] = FRACTION_LEVELS,
    rgb: np.ndarray = FRACTION_RGB,
    cbar_label: str = "",
    cbar_ticks: Optional[Sequence[float]] = None,
    cbar_ticklabels: Optional[Sequence[str]] = None,
) -> None:
    """Plot a 12-panel monthly climatology from a tile time-series file.

    This is the same machinery used by LAI.  GREEN and NDVI are easy package
    additions because they share the same F77 time-series layout already used
    for GREEN.mp4 and merged_Z0 diagnostics.
    """
    path = base_dir / filename
    if not path.exists():
        print(f"Skipping {outname}; missing {path}")
        return
    if layout is None:
        print(f"Skipping {outname}; could not determine time-series layout for {path}")
        return
    monthly = monthly_means_from_interpolated(path, ncat, layout, kind=kind)
    names = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    titles = [f"{mon}" for mon in names]
    grids = [vector_to_grid(tile_id, monthly[m, :]) for m in range(12)]
    panel_continuous_shared_colorbar(
        grids,
        titles,
        lon,
        lat,
        limits,
        outdir / outname,
        ncols=3,
        levels=levels,
        rgb=rgb,
        figsize=(13.5, 9.2),
        coastlines=coastlines,
        cbar_label=(cbar_label or product_label),
        cbar_ticks=cbar_ticks,
        cbar_ticklabels=cbar_ticklabels,
        cbar_tick_rotation=0.0 if cbar_ticks is not None else 90.0,
        bad_color=bad_color,
    )    


def plot_lai(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, ncat: int, layout: Optional[TimeSeriesLayout], coastlines: bool) -> None:
    plot_monthly_timeseries(
        base_dir, tile_id, lon, lat, limits, outdir, ncat,
        "lai.dat", "lai.jpg", "LAI", layout, coastlines,
        kind="lai", levels=LAI_LEVELS, rgb=LAI_PLOT_RGB, cbar_label="LAI",
        cbar_ticks=LAI_TICKS,
        cbar_ticklabels=LAI_TICK_LABELS,
    )    


def plot_green(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, ncat: int, layout: Optional[TimeSeriesLayout], coastlines: bool) -> None:
    plot_monthly_timeseries(
        base_dir, tile_id, lon, lat, limits, outdir, ncat,
        "green.dat", "green.jpg", "GREEN", layout, coastlines,
        kind="fraction", levels=FRACTION_LEVELS, rgb=FRACTION_RGB,
        cbar_label="Green vegetation fraction",
        cbar_ticks=FRACTION_TICKS,
        cbar_ticklabels=FRACTION_TICK_LABELS,
    )


def plot_ndvi(base_dir: Path, tile_id: np.ndarray, lon: np.ndarray, lat: np.ndarray, limits, outdir: Path, ncat: int, layout: Optional[TimeSeriesLayout], coastlines: bool) -> None:
    plot_monthly_timeseries(
        base_dir, tile_id, lon, lat, limits, outdir, ncat,
        "ndvi.dat", "ndvi.jpg", "NDVI", layout, coastlines,
        kind="ndvi", levels=FRACTION_LEVELS, rgb=FRACTION_RGB,
        cbar_label="NDVI",
        cbar_ticks=FRACTION_TICKS,
        cbar_ticklabels=FRACTION_TICK_LABELS,
    )

def z0_value(z2ch: np.ndarray, lai: np.ndarray, scale4z0: float) -> np.ndarray:
    min_veg_height = 0.01
    z0_by_zveg = 0.13
    if scale4z0 == 2.0:
        return scale4z0 * z0_by_zveg * (z2ch - (z2ch - min_veg_height) * np.exp(-lai))
    return z0_by_zveg * (z2ch - scale4z0 * (z2ch - min_veg_height) * np.exp(-lai))


def read_vegdyn(base_dir: Path, ncat: int, layout: F77Layout) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = base_dir / "vegdyn.data"
    if not path.exists():
        raise ClsmPlotError(f"Missing vegdyn data: {path}")
    if is_netcdf(path):
        ity = read_nc_var(path, "ITY").reshape(-1).astype(np.float32)
        z2 = read_nc_var(path, "Z2CH").reshape(-1).astype(np.float32)
        asz0 = read_nc_var(path, "ASCATZ0").reshape(-1).astype(np.float32)
        return ity[:ncat], z2[:ncat], asz0[:ncat]
    with FortranSequentialReader(path, layout) as rdr:
        ity = rdr.read_array(np.float32, ncat)
        z2 = rdr.read_array(np.float32, ncat)
        asz0 = rdr.read_array(np.float32, ncat)
    return ity, z2, asz0


def plot_canoph_from_vegdyn(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, ncat: int, layout: F77Layout, coastlines: bool) -> Tuple[np.ndarray, np.ndarray]:
    _, z2, asz0 = read_vegdyn(base_dir, ncat, layout)
    grid = vector_to_grid(tile_id, z2)
    levels = np.linspace(float(np.nanmin(z2)), float(np.nanmax(z2)), 17)
    panel_continuous([grid], ["Canopy height Z2CH"], lon, lat, limits, outdir / "Canopy_Height_onTiles.jpg", ncols=1, levels_list=[levels], color_ids_list=[list(reversed(CONTINUOUS_COLOR_IDS))], figsize=(11.2, 5.6), coastlines=coastlines)
    return z2, asz0 * 1000.0


def seasonal_z0_and_ndvi(base_dir: Path, ncat: int, z2ch: np.ndarray, scale4z0: float, lai_layout: TimeSeriesLayout, ndvi_layout: Optional[TimeSeriesLayout] = None) -> Tuple[np.ndarray, np.ndarray]:
    lai_path = base_dir / "lai.dat"
    ndvi_path = base_dir / "ndvi.dat"
    if not lai_path.exists() or not ndvi_path.exists():
        raise ClsmPlotError(f"Missing {lai_path} or {ndvi_path}; required for icarus/merged Z0")
    mdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Accumulate valid daily means only.  Invalid LAI/NDVI should not become
    # zero roughness; otherwise it appears as artificial brown/underflow bins.
    zo_sum = np.zeros((ncat, 4), dtype=np.float64)
    zo_count = np.zeros((ncat, 4), dtype=np.float64)
    ndvi_sum = np.zeros((ncat, 4), dtype=np.float64)
    ndvi_count = np.zeros((ncat, 4), dtype=np.float64)

    if ndvi_layout is None:
        ndvi_layout = lai_layout
    with TimeSeriesReader(lai_path, lai_layout, ncat) as lai_rdr, TimeSeriesReader(ndvi_path, ndvi_layout, ncat) as ndvi_rdr:
        lh1, lv1 = read_timeseries_record(lai_rdr, ncat)
        lh2, lv2 = read_timeseries_record(lai_rdr, ncat)
        nh1, nv1 = read_timeseries_record(ndvi_rdr, ncat)
        nh2, nv2 = read_timeseries_record(ndvi_rdr, ncat)
        lb4, lnxt = midpoint_doy(lh1), midpoint_doy(lh2)
        nb4, nnxt = midpoint_doy(nh1), midpoint_doy(nh2)
        current_year_offset = _initial_loop_year_offset(lh1, lh2)
        for month in range(1, 13):
            season = 0 if month in (12, 1, 2) else 1 if month in (3, 4, 5) else 2 if month in (6, 7, 8) else 3
            for day in range(1, mdays[month - 1] + 1):
                now = float((_dt.date(2001 + current_year_offset, month, day) - _dt.date(2000, 12, 31)).days)
                lden = (lnxt - lb4) if abs(lnxt - lb4) > 1e-6 else 1.0
                nden = (nnxt - nb4) if abs(nnxt - nb4) > 1e-6 else 1.0
                lai = ((now - lb4) / lden) * lv2 + ((lnxt - now) / lden) * lv1
                ndvi = ((now - nb4) / nden) * nv2 + ((nnxt - now) / nden) * nv1
                lai = _sanitize_lai(lai)
                ndvi = _sanitize_ndvi(ndvi)
                zot_m = 1000.0 * z0_value(z2ch, lai, scale4z0)
                zot_m = _sanitize_z0_mm(zot_m)

                good_z = np.isfinite(zot_m)
                if np.any(good_z):
                    zo_sum[good_z, season] += zot_m[good_z]
                    zo_count[good_z, season] += 1.0
                good_n = np.isfinite(ndvi)
                if np.any(good_n):
                    ndvi_sum[good_n, season] += ndvi[good_n]
                    ndvi_count[good_n, season] += 1.0

                if now + 0.5 >= lnxt:
                    lv1 = lv2
                    lb4 = lnxt
                    try:
                        lh2, lv2 = read_timeseries_record(lai_rdr, ncat)
                        lnxt = midpoint_doy(lh2)
                        current_year_offset = _advance_loop_year_offset(lh2, month)
                    except EOFError:
                        lnxt = now + 9999.0
                if now + 0.5 >= nnxt:
                    nv1 = nv2
                    nb4 = nnxt
                    try:
                        nh2, nv2 = read_timeseries_record(ndvi_rdr, ncat)
                        nnxt = midpoint_doy(nh2)
                    except EOFError:
                        nnxt = now + 9999.0
    zo_vec_mm = np.full((ncat, 4), np.nan, dtype=np.float32)
    ndvi_vec = np.full((ncat, 4), np.nan, dtype=np.float32)
    goodz = zo_count > 0.0
    goodn = ndvi_count > 0.0
    zo_vec_mm[goodz] = (zo_sum[goodz] / zo_count[goodz]).astype(np.float32)
    ndvi_vec[goodn] = (ndvi_sum[goodn] / ndvi_count[goodn]).astype(np.float32)
    return zo_vec_mm, ndvi_vec


def plot_z0(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, ncat: int, z2: np.ndarray, asz0_mm: np.ndarray, layout: Optional[TimeSeriesLayout], coastlines: bool, products: Sequence[str], ndvi_layout: Optional[TimeSeriesLayout] = None) -> None:
    scale4z0 = 0.5 if (base_dir / "CLM_veg_typs_fracs").exists() else 2.0
    need_lai = any(p in ("icarus", "merged") for p in products)
    zo_vec_mm = ndvi_vec = None
    if need_lai:
        if layout is None:
            print("Skipping icarus/merged Z0: could not determine LAI/NDVI time-series layout")
        else:
            try:
                zo_vec_mm, ndvi_vec = seasonal_z0_and_ndvi(base_dir, ncat, z2, scale4z0, layout, ndvi_layout)
            except Exception as exc:
                print(f"Skipping icarus/merged Z0: {exc}")
    colors = Z0_COLOR_IDS
    levels = Z0_LEVELS
    sea_label = ["DJF", "MAM", "JJA", "SON"]
    seasons_to_plot = [0, 2]
    asz0_mm = _sanitize_z0_mm(asz0_mm)
    for pname in products:
        grids: List[np.ndarray] = []
        titles: List[str] = []
        for season in seasons_to_plot:
            if pname == "ascat":
                data = asz0_mm.copy()
            elif pname == "icarus" and zo_vec_mm is not None:
                data = _sanitize_z0_mm(zo_vec_mm[:, season])
            elif pname == "merged" and zo_vec_mm is not None and ndvi_vec is not None:
                icarus = _sanitize_z0_mm(zo_vec_mm[:, season])
                ndvi = _sanitize_ndvi(ndvi_vec[:, season])
                data = icarus.copy()
                # IDL uses ASZ0 where seasonal NDVI <= 0.2.  Also use ASZ0
                # when NDVI or icarus is invalid so bad seasonal values do not
                # appear as artificial lowest-bin/brown areas.
                use_ascat = (~np.isfinite(data)) | (~np.isfinite(ndvi)) | (ndvi <= 0.2)
                data[use_ascat] = asz0_mm[use_ascat]
                data = _sanitize_z0_mm(data)
            else:
                continue
            grids.append(vector_to_grid(tile_id, data))
            titles.append(f"{pname}: {sea_label[season]}")
        if grids:
            panel_continuous_shared_colorbar(
                grids,
                titles,
                lon,
                lat,
                limits,
                outdir / f"{pname}_Z0.jpg",
                ncols=1,
                levels=levels,
                color_ids=colors,
                figsize=(11.0, 8.0),
                coastlines=coastlines,
                cbar_label="Z0 (mm)",
                cbar_ticks=Z0_LEVELS,
                cbar_ticklabels=Z0_TICK_LABELS,
                cbar_tick_rotation=45.0,
            )



# -----------------------------------------------------------------------------
# Irrigation products - NOT TESTED IN Python Package due to file unavailability.
# -----------------------------------------------------------------------------


def plot_irrig_method(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "irrig.dat"
    if not path.exists():
        print(f"Skipping irrigation method; missing {path}")
        return
    names = ["SPRINKLERFR", "DRIPFR", "FLOODFR"]
    titles = ["SPRINKLER FRACTION", "DRIP FRACTION", "FLOOD FRACTION"]
    grids = []
    for name in names:
        vec = read_nc_var(path, name).reshape(-1).astype(np.float32)
        vec[vec > 1.0] = np.nan
        g = vector_to_grid(tile_id, vec)
        g[g == 0.0] = np.nan
        grids.append(g)
    levels = [np.arange(21) * 0.05] * 3
    panel_continuous(grids, titles, lon, lat, limits, outdir / "IrrigMethod.png", ncols=1, levels_list=levels, color_ids_list=[list(range(140, 161))] * 3, figsize=(9, 11), coastlines=coastlines)


def plot_lai_minmax(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "irrig.dat"
    if not path.exists():
        print(f"Skipping LAI min/max; missing {path}")
        return
    grids = []
    for name in ["LAIMIN", "LAIMAX"]:
        vec = read_nc_var(path, name).reshape(-1).astype(np.float32)
        vec[vec > 100.0] = np.nan
        g = vector_to_grid(tile_id, vec)
        g[g == 0.0] = np.nan
        grids.append(g)
    panel_continuous(grids, ["LAI Minimum", "LAI Maximum"], lon, lat, limits, outdir / "LAI_minmax.png", ncols=1, levels_list=[LAI_LEVELS, LAI_LEVELS], rgb_list=[LAI_PLOT_RGB, LAI_PLOT_RGB], figsize=(9, 11), coastlines=coastlines)


def plot_irrig_fractions(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, coastlines: bool) -> None:
    path = base_dir / "irrig.dat"
    if not path.exists():
        print(f"Skipping irrigation fractions; missing {path}")
        return
    names = ["IRRIGFRAC", "PADDYFRAC", "RAINFEDFRAC"]
    titles = ["IRRIGATED CROP FRACTION", "PADDY FRACTION", "RAINFED FRACTION"]
    grids = []
    for name in names:
        vec = read_nc_var(path, name).reshape(-1).astype(np.float32)
        vec[vec > 1.0] = np.nan
        g = vector_to_grid(tile_id, vec)
        g[g == 0.0] = np.nan
        grids.append(g)
    levels = [np.arange(21) * 0.025] * 3
    panel_continuous(grids, titles, lon, lat, limits, outdir / "GIA-Hybrid_IrrigFracs.png", ncols=1, levels_list=levels, color_ids_list=[list(range(140, 161))] * 3, figsize=(9, 11), coastlines=coastlines)


def _decode_crop_names(arr: np.ndarray) -> List[str]:
    if arr.dtype.kind in {"S", "U"}:
        return [str(x).strip(" b'\x00") for x in arr.reshape(-1)]
    # xarray sometimes returns char arrays.  Try joining along the last dimension.
    if arr.ndim >= 2 and arr.dtype.kind in {"S", "U"}:
        return ["".join(map(str, row)).strip() for row in arr.reshape(arr.shape[0], -1)]
    return [f"crop_{i + 1:02d}" for i in range(26)]


def plot_crop_times(base_dir: Path, tile_id: np.ndarray, lon, lat, limits, outdir: Path, coastlines: bool) -> None:
    """A compact Python version of plot_crop_times.

    The IDL routine creates multiple 4-row pages with crop fraction, planting
    day, harvest day, and irrigation type for up to 26 crops.  This version keeps
    that structure but uses raster panels for fraction and scatter panels for DOY
    and irrigation type.
    """
    path = base_dir / "irrig.dat"
    if not path.exists():
        print(f"Skipping crop times; missing {path}")
        return
    if xr is None:
        print("Skipping crop times; xarray is unavailable")
        return
    with xr.open_dataset(path, decode_times=False) as ds:
        required = ["IRRIGPLANT", "IRRIGHARVEST", "CROPIRRIGFRAC", "IRRIGTYPE"]
        missing = [v for v in required if v not in ds]
        if missing:
            print(f"Skipping crop times; missing variables in {path}: {missing}")
            return
        plantv = np.asarray(ds["IRRIGPLANT"].values)
        harvestv = np.asarray(ds["IRRIGHARVEST"].values)
        fracv = np.asarray(ds["CROPIRRIGFRAC"].values)
        irrigtypev = np.asarray(ds["IRRIGTYPE"].values)
        crop_names = _decode_crop_names(np.asarray(ds["CROPCLASSNAME"].values)) if "CROPCLASSNAME" in ds else [f"crop_{i + 1:02d}" for i in range(26)]

    ncat = int(np.nanmax(tile_id))
    # Normalize expected dimensions so tile is first, crop is last where possible.
    plantv = np.asarray(plantv)
    harvestv = np.asarray(harvestv)
    fracv = np.asarray(fracv)
    irrigtypev = np.asarray(irrigtypev)
    # Best effort: assume first dimension is tile.  If not, try moving the ncat dimension first.
    def move_tile_first(a: np.ndarray) -> np.ndarray:
        for ax, size in enumerate(a.shape):
            if size >= ncat:
                return np.moveaxis(a, ax, 0)
        return a
    plantv = move_tile_first(plantv)
    harvestv = move_tile_first(harvestv)
    fracv = move_tile_first(fracv)
    irrigtypev = move_tile_first(irrigtypev)
    ncrops = min(26, fracv.shape[-1])
    levels_frac = np.arange(21) * 0.025
    doy_levels = np.array([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366, 370])
    page = 1
    row_in_page = 0
    fig = None

    def new_page():
        return plt.figure(figsize=(13, 11))

    for crop in range(ncrops):
        # Extract tile vectors.  Common shapes are (tile, season, crop) for plant/harvest and (tile, crop) for frac/type.
        frac = fracv[:ncat, crop].astype(np.float32) if fracv.ndim == 2 else fracv[:ncat, ..., crop].reshape(ncat, -1)[:, 0].astype(np.float32)
        ityp = irrigtypev[:ncat, crop].astype(np.float32) if irrigtypev.ndim == 2 else irrigtypev[:ncat, ..., crop].reshape(ncat, -1)[:, 0].astype(np.float32)
        plant = plantv[:ncat, 0, crop].astype(np.float32) if plantv.ndim >= 3 else plantv[:ncat, crop].astype(np.float32)
        harvest = harvestv[:ncat, 0, crop].astype(np.float32) if harvestv.ndim >= 3 else harvestv[:ncat, crop].astype(np.float32)
        if not np.isfinite(frac).any() or np.nanmax(frac) <= 0:
            continue
        if fig is None:
            fig = new_page()
        panels = [vector_to_grid(tile_id, frac), vector_to_grid(tile_id, plant), vector_to_grid(tile_id, harvest), vector_to_grid(tile_id, ityp)]
        titles = ["frac", "DOY plant", "DOY harvest", "IRRIGTYPE"]
        levels = [levels_frac, doy_levels, doy_levels, [1, 2, 3, 4]]
        color_ids = [list(range(140, 161)), [69, 145, 64, 66, 70, 71, 73, 75, 76, 78, 80, 113, 114, 116, 117], [69, 145, 64, 66, 70, 71, 73, 75, 76, 78, 80, 113, 114, 116, 117], [69, 64, 80, 255]]
        for col in range(4):
            ax = make_axes(fig, 4, 4, row_in_page * 4 + col + 1, coastlines)
            plot_continuous_on_ax(ax, panels[col], lon, lat, limits, f"{crop_names[crop] if crop < len(crop_names) else crop}: {titles[col]}", levels[col], color_ids[col], coastlines=coastlines)
        row_in_page += 1
        if row_in_page == 4:
            save_fig(fig, outdir / f"gia_irrig_params_{page:02d}.png")
            fig = None
            row_in_page = 0
            page += 1
    if fig is not None:
        save_fig(fig, outdir / f"gia_irrig_params_{page:02d}.png")


# -----------------------------------------------------------------------------
# Movies
# -----------------------------------------------------------------------------


def _format_movie_tick_label(x: float) -> str:
    """Compact numeric labels for per-frame movie colorbars."""
    x = float(x)
    if abs(x) < 1.0e-10:
        return "0"
    if abs(x - round(x)) < 1.0e-10:
        return str(int(round(x)))
    if abs(x) < 0.1:
        return f"{x:.3f}".rstrip("0").rstrip(".")
    if abs(x) < 1.0:
        return f"{x:.2f}".rstrip("0").rstrip(".")
    return f"{x:g}"


def movie_colorbar_ticks(vname: str) -> Tuple[np.ndarray, List[str], str]:
    """Return movie colorbar ticks consistent with static plots."""
    label_map = {
        "GREEN": "Green vegetation fraction",
        "VISDF": "VIS diffuse albedo",
        "NIRDF": "NIR diffuse albedo",
        "NDVI": "NDVI",
    }

    if vname == "LAI":
        return LAI_TICKS, list(LAI_TICK_LABELS), "LAI"

    return FRACTION_TICKS, list(FRACTION_TICK_LABELS), label_map.get(vname, vname)

def add_movie_colorbar(fig: plt.Figure, ax, sm: ScalarMappable, vname: str, levels: Sequence[float]) -> None:
    """Add one horizontal colorbar to every movie frame.

    The movie frame is captured from the raw canvas, so the colorbar must be
    drawn into the figure before ``buffer_rgba`` is read.  Use fixed ticks so
    every frame has a stable scale and readable labels.
    """
    ticks, labels, label = movie_colorbar_ticks(vname)
    cax = fig.add_axes([0.08, 0.075, 0.84, 0.030])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=ticks, spacing="uniform")
    cbar.ax.xaxis.set_major_locator(FixedLocator(ticks))
    cbar.ax.xaxis.set_major_formatter(FixedFormatter(labels))
    cbar.ax.tick_params(labelsize=5, rotation=0, pad=2)
    cbar.set_label(label, fontsize=8, labelpad=3)

def make_movie(
    base_dir: Path,
    rst_file: Path,
    nc: int,
    nr: int,
    ncat: int,
    layout: Optional[TimeSeriesLayout],
    rst_layout: F77Layout,
    outdir: Path,
    gfile: str,
    vname: str,
    nc_movie: int,
    nr_movie: int,
    limits: Tuple[float, float, float, float],
    cache_dir: Path,
    coastlines: bool,
) -> None:
    if layout is None:
        print(f"Skipping {vname} movie; could not determine time-series layout")
        return
    mapping_cache = cache_dir / f"fractional_{gfile}_{nc_movie}x{nr_movie}.npz"
    # The movie aggregation map is built from the integer raster (.rst), so it
    # must use the raster F77 layout.  The seasonal LAI/GREEN/AlbMap layout is
    # a different object and does not have marker_dtype.  Passing it here caused
    # movie-only runs to fail with: 'TimeSeriesLayout' object has no attribute
    # 'marker_dtype'.
    mat = build_fractional_sparse_from_rst(rst_file, nc, nr, ncat, nc_movie, nr_movie, rst_layout, mapping_cache)
    # Rows with no contributing land/catchment tiles are ocean/no-data.
    # Sparse matrix multiplication returns 0.0 for empty rows, which would
    # otherwise be plotted as a valid zero/low-value color in movies.  Static
    # climatology plots get NaN from vector_to_grid() for invalid cells; do
    # the equivalent here so cmap.set_bad(NO_DATA_COLOR) is used consistently.
    movie_cell_has_data = np.asarray(mat.getnnz(axis=1)).ravel() > 0
    filename_map = {
        "LAI": "lai.dat",
        "GREEN": "green.dat",
        "VISDF": "AlbMap.WS.8-day.tile.0.3_0.7.dat",
        "NIRDF": "AlbMap.WS.8-day.tile.0.7_5.0.dat",
        "NDVI": "ndvi.dat",
    }
    path = base_dir / filename_map[vname]
    if not path.exists():
        print(f"Skipping {vname} movie; missing {path}")
        return
    levels = LAI_LEVELS if vname == "LAI" else FRACTION_LEVELS
    rgb = LAI_PLOT_RGB if vname == "LAI" else FRACTION_RGB
    bad_color = NO_DATA_COLOR    
    lon, lat = lon_lat_centers(nc_movie, nr_movie)
    mdays = [31,28,31,30,31,30,31,31,30,31,30,31]
    outpath = outdir / f"{vname}.mp4"
    print(f"Writing movie {outpath}")
    with TimeSeriesReader(path, layout, ncat) as rdr:
        h1, v1 = read_timeseries_record(rdr, ncat)
        h2, v2 = read_timeseries_record(rdr, ncat)
        b4 = midpoint_doy(h1)
        nxt = midpoint_doy(h2)
        # Use the same year-offset logic and sanitation as the LAI/Z0 static plots.
        # Using h1[0] here can extrapolate January from the wrong year for Dec->Jan
        # climatology records, producing negative LAI in the movie even when lai.jpg is OK.
        current_year_offset = _initial_loop_year_offset(h1, h2)
        with open_mp4_writer(outpath, fps=10) as writer:
            for month in range(1, 13):
                for day in range(1, mdays[month - 1] + 1):
                    now = float((_dt.date(2001 + current_year_offset, month, day) - _dt.date(2000, 12, 31)).days)
                    denom = (nxt - b4) if abs(nxt - b4) > 1e-6 else 1.0
                    vec = ((now - b4) / denom) * v2 + ((nxt - now) / denom) * v1
                    if vname == "LAI":
                        vec = _sanitize_lai(vec)
                    elif vname == "NDVI":
                        vec = _sanitize_ndvi(vec)
                    else:
                        # GREEN, VISDF, and NIRDF are fractional fields.  Keep real
                        # values in [0,1] and mask impossible/fill values.
                        vec = np.asarray(vec, dtype=np.float32)
                        bad = (~np.isfinite(vec)) | (vec < 0.0) | (vec > 1.0) | (np.abs(vec) > 1.0e10)
                        vec = vec.copy()
                        vec[bad] = np.nan
                    flat = np.asarray(mat @ vec.astype(np.float32), dtype=np.float32).ravel()
                    flat[~movie_cell_has_data] = np.nan
                    grid = flat.reshape(nr_movie, nc_movie)
                    fig = plt.figure(figsize=(7.8, 5.85), dpi=100)
                    # Leave room for lon/lat tick labels and a fixed colorbar.
                    # Movie frames are captured from the raw canvas, not through
                    # savefig(..., bbox_inches="tight"), so margins must be explicit.
                    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.205, top=0.90)
                    ax = make_axes(fig, 1, 1, 1, coastlines)
                    date_stamp = f"{2001 + current_year_offset:04d}{month:02d}{day:02d}"
                    sm = plot_continuous_on_ax(ax, grid, lon, lat, limits, f"{vname}: {date_stamp}", levels, rgb=rgb, coastlines=coastlines, bad_color=bad_color)                    
                    add_movie_colorbar(fig, ax, sm, vname, levels)
                    fig.canvas.draw()
                    frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
                    writer.append_data(frame)
                    plt.close(fig)
                    if now + 0.5 >= nxt:
                        v1 = v2
                        b4 = nxt
                        try:
                            h2, v2 = read_timeseries_record(rdr, ncat)
                            nxt = midpoint_doy(h2)
                            current_year_offset = _advance_loop_year_offset(h2, month)
                        except EOFError:
                            nxt = now + 9999.0
    print(f"Wrote {outpath}")


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------


DEFAULT_PLOTS = [
    "tiles", "country", "cti", "mosaic", "clm", "carbon", "ndep", "soil",
    "elevation", "lai", "green", "ndvi", "canopy", "z0"
]
MOVIE_PLOTS = ["movies"]
LEGACY_PLOTS = DEFAULT_PLOTS + MOVIE_PLOTS
# Irrigation products - NOT TESTED IN Python Package due to file unavailability.
# These routines exist as legacy/experimental helpers, but the active legacy IDL
# driver does not call them by default.  Keep them explicitly requestable without
# making --plots all unexpectedly produce unvalidated products.
EXPERIMENTAL_PLOTS = ["irrig_method", "lai_minmax", "irrig_fractions", "crop_times"]
VALID_PLOTS = LEGACY_PLOTS + EXPERIMENTAL_PLOTS
# User-facing "all" is intentionally current legacy parity, not experimental extras.
ALL_PLOTS = LEGACY_PLOTS


def parse_plot_list(text: str) -> List[str]:
    text = text.strip().lower()
    if text in ("default", "main", "fixed", "images"):
        return DEFAULT_PLOTS.copy()
    if text in ("legacy", "legacy_idl", "legacy-idl", "idl", "idl_default", "idl-default"):
        return LEGACY_PLOTS.copy()
    if text in ("all", "everything"):
        return ALL_PLOTS.copy()
    if text in ("quick", "smoke"):
        return ["tiles", "cti", "elevation"]
    if text in ("experimental", "extras"):
        return EXPERIMENTAL_PLOTS.copy()
    plots = []
    aliases = {
        "veg": "mosaic",
        "ndep_t2m": "ndep",
        "soilalb": "ndep",
        "irrig": "irrig_fractions",
        "irrigation": "irrig_fractions",
    }
    for item in text.split(","):
        item = item.strip().lower().replace("-", "_")
        if not item:
            continue
        plots.append(aliases.get(item, item))
    unknown = [p for p in plots if p not in VALID_PLOTS]
    if unknown:
        raise ClsmPlotError(
            f"Unknown plot option(s): {unknown}. "
            f"Valid modes: quick, default, movies, legacy, all, experimental. "
            f"Valid explicit items: {VALID_PLOTS}"
        )
    return plots


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Drop-in Python replacement for IDL clsm_plots.pro")
    p.add_argument("--gfile", default=os.environ.get("gfile"), help="Grid/file stem used to find workdir/rst/<gfile>*.rst. Defaults to $gfile.")
    p.add_argument("--workdir", default=os.environ.get("workdir"), help="BCS work directory containing rst/. Defaults to $workdir.")
    p.add_argument("--nc", type=int, default=int(os.environ["NC"]) if os.environ.get("NC") else None, help="Full raster NC. Defaults to $NC.")
    p.add_argument("--nr", type=int, default=int(os.environ["NR"]) if os.environ.get("NR") else None, help="Full raster NR. Defaults to $NR.")
    p.add_argument("--base-dir", default="..", help="Directory containing catchment.def, cti_stats.dat, soil_param.dat, etc. Default: ..")
    p.add_argument("--outdir", default=".", help="Directory for plot outputs. Default: current directory.")
    p.add_argument("--plots", default=os.environ.get("CLSM_PLOTS", "legacy"), help=f"Comma list, or quick/default/movies/legacy/all/experimental. No-argument default is legacy (current legacy IDL-equivalent outputs: fixed JPGs + movies); override with $CLSM_PLOTS. Valid explicit items: {','.join(VALID_PLOTS)}")
    p.add_argument("--plot-nc", type=int, default=4320, help="Output longitude cells for the main tile map. IDL default: 4320")
    p.add_argument("--plot-nr", type=int, default=2160, help="Output latitude cells for the main tile map. IDL default: 2160")
    p.add_argument("--movie-nc", type=int, default=720, help="Movie longitude cells. IDL default: 720")
    p.add_argument("--movie-nr", type=int, default=360, help="Movie latitude cells. IDL default: 360")
    p.add_argument("--dpi", type=int, default=int(os.environ.get("CLSM_PLOT_DPI", "180")), help="DPI for static JPG plots. Default: 180; may also be set with $CLSM_PLOT_DPI.")
    p.add_argument("--jpeg-quality", type=int, default=int(os.environ.get("CLSM_JPEG_QUALITY", "95")), help="JPEG quality for static JPG plots. Default: 95; may also be set with $CLSM_JPEG_QUALITY.")
    p.add_argument("--endian", choices=["auto", "little", "big", "<", ">"], default="auto", help="Endian for F77 binary files. Default: auto")
    p.add_argument("--record-marker", type=int, choices=[0, 4, 8], default=0, help="F77 record marker bytes. 0 means auto; otherwise 4 or 8.")
    p.add_argument("--cache-dir", default=None, help="Directory for tile/mapping caches. Default: <workdir>/.clsm_plot_cache, so cache is not moved into final clsm/plots.")
    p.add_argument("--rst-file", default=None, help="Explicit raster file to use instead of auto-selecting workdir/rst/<gfile>*.rst. Useful when both <gfile>.rst and <gfile>-Pfafstetter.rst are present.")
    p.add_argument("--tile-source", choices=["auto", "rst", "catchment"], default=os.environ.get("CLSM_TILE_SOURCE", "auto"), help="How to build the plotting tile map. rst reproduces the IDL raster path; catchment builds lon/lat boxes directly from catchment.def; auto uses rst unless its IDs fail a catchment.def spatial-consistency check. Default: auto; may also be set with $CLSM_TILE_SOURCE")
    p.add_argument("--no-cache", action="store_true", help="Do not read/write tile-map cache")
    p.add_argument("--coastlines", dest="coastlines", action="store_true", default=(ccrs is not None), help="Draw coastlines if Cartopy is installed. Default: on when Cartopy is available.")
    p.add_argument("--no-coastlines", dest="coastlines", action="store_false", help="Disable Cartopy coastlines even if Cartopy is installed.")
    p.add_argument("--z0-products", default="ascat,icarus,merged", help="Comma list of Z0 products: ascat,icarus,merged")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    global PLOT_DPI, JPEG_QUALITY
    PLOT_DPI = int(args.dpi)
    JPEG_QUALITY = int(args.jpeg_quality)
    missing = [name for name in ("gfile", "workdir", "nc", "nr") if getattr(args, name) in (None, "")]
    if missing:
        raise ClsmPlotError(f"Missing required settings: {missing}. Provide args or set gfile/workdir/NC/NR environment variables.")

    gfile = str(args.gfile)
    workdir = Path(str(args.workdir)).expanduser().resolve()
    base_dir = Path(args.base_dir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir).expanduser().resolve() if args.cache_dir else workdir / ".clsm_plot_cache"
    plots = parse_plot_list(args.plots)

    ncat = read_ncat(base_dir)
    limits = read_limits(base_dir, gfile)
    rst_file, layout = select_rst_file(
        workdir, gfile, ncat, int(args.nc), int(args.nr), args.endian, args.record_marker, args.rst_file
    )
    print(f"ncat={ncat}; limits={limits}; F77 layout=endian {layout.endian}, marker {layout.marker_bytes} bytes")

    lon, lat = lon_lat_centers(args.plot_nc, args.plot_nr)
    tile_cache = None if args.no_cache else cache_dir / f"tile_id_{gfile}_{args.plot_nc}x{args.plot_nr}.npz"
    catch_cache = None if args.no_cache else cache_dir / f"tile_id_from_catchment_def_{gfile}_{args.plot_nc}x{args.plot_nr}.npz"

    if args.tile_source == "catchment":
        tile_id = build_tile_id_from_catchment_def(base_dir, ncat, args.plot_nc, args.plot_nr, catch_cache)
    else:
        tile_id = build_tile_id_from_rst(rst_file, int(args.nc), int(args.nr), ncat, args.plot_nc, args.plot_nr, layout, tile_cache)
        match = catchment_spatial_match_fraction(tile_id, base_dir, lon, lat, ncat)
        print(f"RST/catchment.def spatial match fraction: {match:.6f}")
        if args.tile_source == "auto" and match < 0.05:
            print(
                "RST tile IDs do not spatially match catchment.def boxes well; "
                "falling back to catchment.def-derived plotting tile map. "
                "Use --tile-source rst to force IDL-style rst mapping."
            )
            tile_id = build_tile_id_from_catchment_def(base_dir, ncat, args.plot_nc, args.plot_nr, catch_cache)

    valid_tile_cells = int(np.count_nonzero((tile_id >= 1) & (tile_id <= ncat)))
    total_tile_cells = int(tile_id.size)
    print(f"Valid plotting cells: {valid_tile_cells}/{total_tile_cells}")
    if valid_tile_cells == 0:
        raise ClsmPlotError(
            "Tile map has zero valid CLSM tile ids. This usually means the wrong rst file was used "
            "or catchment.def boxes could not be mapped to the plotting grid. Delete the cache, "
            "rerun with --no-cache, or try --tile-source catchment / --rst-file <file>."
        )

    if "tiles" in plots:
        plot_tiles(
            tile_id, lon, lat, limits, outdir, args.coastlines,
            rst_file, int(args.nc), int(args.nr), ncat, layout,
            gfile=gfile,
        )        
    if "country" in plots:
        plot_country_codes(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "cti" in plots:
        plot_cti(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "mosaic" in plots:
        plot_mosaic(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "clm" in plots:
        plot_clm(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "carbon" in plots:
        plot_carbon(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "ndep" in plots:
        plot_ndep_t2m_soilalb(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "soil" in plots:
        plot_soil(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "elevation" in plots:
        plot_elevation(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "lai" in plots:
        try:
            ts_layout = choose_timeseries_layout(base_dir / "lai.dat", ncat, args.endian, args.record_marker) if (base_dir / "lai.dat").exists() else None
            plot_lai(base_dir, tile_id, lon, lat, limits, outdir, ncat, ts_layout, args.coastlines)
        except Exception as exc:
            print(f"Skipping lai.jpg after error: {exc}")
    if "green" in plots:
        try:
            green_layout = choose_timeseries_layout(base_dir / "green.dat", ncat, args.endian, args.record_marker) if (base_dir / "green.dat").exists() else None
            plot_green(base_dir, tile_id, lon, lat, limits, outdir, ncat, green_layout, args.coastlines)
        except Exception as exc:
            print(f"Skipping green.jpg after error: {exc}")
    if "ndvi" in plots:
        try:
            ndvi_layout_static = choose_timeseries_layout(base_dir / "ndvi.dat", ncat, args.endian, args.record_marker) if (base_dir / "ndvi.dat").exists() else None
            plot_ndvi(base_dir, tile_id, lon, lat, limits, outdir, ncat, ndvi_layout_static, args.coastlines)
        except Exception as exc:
            print(f"Skipping ndvi.jpg after error: {exc}")

    z2 = asz0_mm = None
    if "canopy" in plots or "z0" in plots:
        try:
            z2, asz0_mm = plot_canoph_from_vegdyn(base_dir, tile_id, lon, lat, limits, outdir, ncat, layout, args.coastlines)
        except Exception as exc:
            print(f"Skipping Canopy_Height_onTiles.jpg / ASZ0 inputs after error: {exc}")
    if "z0" in plots and z2 is not None and asz0_mm is not None:
        z0_products = [x.strip().lower() for x in args.z0_products.split(",") if x.strip()]
        ts_layout = None
        ndvi_layout = None
        if (base_dir / "lai.dat").exists():
            try:
                ts_layout = choose_timeseries_layout(base_dir / "lai.dat", ncat, args.endian, args.record_marker)
            except Exception as exc:
                print(f"Could not determine LAI layout for icarus/merged Z0: {exc}")
        if (base_dir / "ndvi.dat").exists():
            try:
                ndvi_layout = choose_timeseries_layout(base_dir / "ndvi.dat", ncat, args.endian, args.record_marker)
            except Exception as exc:
                print(f"Could not determine NDVI layout for icarus/merged Z0: {exc}")
        try:
            plot_z0(base_dir, tile_id, lon, lat, limits, outdir, ncat, z2, asz0_mm, ts_layout, args.coastlines, z0_products, ndvi_layout=ndvi_layout)
        except Exception as exc:
            print(f"Skipping Z0 plots after error: {exc}")

    if "irrig_method" in plots:
        plot_irrig_method(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "lai_minmax" in plots:
        plot_lai_minmax(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "irrig_fractions" in plots:
        plot_irrig_fractions(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)
    if "crop_times" in plots:
        plot_crop_times(base_dir, tile_id, lon, lat, limits, outdir, args.coastlines)

    if "movies" in plots:
        for vname in ["LAI", "GREEN", "VISDF", "NIRDF", "NDVI"]:
            try:
                movie_layout = choose_timeseries_layout(base_dir / {"LAI":"lai.dat", "GREEN":"green.dat", "VISDF":"AlbMap.WS.8-day.tile.0.3_0.7.dat", "NIRDF":"AlbMap.WS.8-day.tile.0.7_5.0.dat", "NDVI":"ndvi.dat"}[vname], ncat, args.endian, args.record_marker) if (base_dir / {"LAI":"lai.dat", "GREEN":"green.dat", "VISDF":"AlbMap.WS.8-day.tile.0.3_0.7.dat", "NIRDF":"AlbMap.WS.8-day.tile.0.7_5.0.dat", "NDVI":"ndvi.dat"}[vname]).exists() else None
                make_movie(base_dir, rst_file, int(args.nc), int(args.nr), ncat, movie_layout, layout, outdir, gfile, vname, args.movie_nc, args.movie_nr, limits, cache_dir, args.coastlines)
            except Exception as exc:
                print(f"Skipping {vname}.mp4 after error: {exc}")

    print("Done.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ClsmPlotError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)

