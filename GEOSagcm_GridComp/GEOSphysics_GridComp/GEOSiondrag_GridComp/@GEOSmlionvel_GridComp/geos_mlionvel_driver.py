"""
Title: GEOS ML Ion Velocity Runtime Driver
Purpose: Run the trained UI/VI ion-velocity model through MAPL_PythonBridge
         and remap predicted ion winds to the current GEOS vertical grid.
Author: Andrew Lee
Affiliation: NASA/CUA
Email: andrew.won.lee@nasa.gov
Date Created: 2026-09-30
Last Modified: 2026-09-30
"""

import math
import os
import sys
import traceback

import numpy as np
import torch
import torch.nn as nn

from MAPL_PythonBridge import UserCode, get_MAPLPy


# ==================
# USER CONFIG
# ==================

# PyTorch runtime settings. CPU inference matches the current ML-radiation path.
TORCH_NUM_THREADS = 1
TORCH_INTEROP_THREADS = 1
COLUMN_BATCH_SIZE = 1024
DEVICE = "cpu"

# Runtime diagnostics.
PRINT_MODEL_SUMMARY = True
PRINT_OUTPUT_MINMAX = True

# Runtime data directory. The default assumes this driver is installed in
# <GEOS_PREFIX>/lib/Python and model data are installed in
# <GEOS_PREFIX>/share/GEOSmlionvel.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = os.path.abspath(
    os.path.join(SCRIPT_DIR, "..", "..", "share", "GEOSmlionvel")
)
DATA_DIR = os.path.abspath(
    os.environ.get("GEOS_MLIONVEL_DATA_DIR", DEFAULT_DATA_DIR)
)

MODEL_PATH = os.path.join(
    DATA_DIR,
    "models",
    "mlionvel_uivi_colfilm_2020-2022.pth",
)
NORM_PATH = os.path.join(
    DATA_DIR,
    "norms",
    "NormStats_IonVel_UIVI_2020-2022_train.npz",
)

# The GEOS run script stages this file into the runtime directory.
F107_AP_PATH = os.path.abspath(
    os.environ.get("GEOS_MLIONVEL_F107_AP_PATH", "F107_ap_appended.txt")
)

# Output field names exposed by GEOS_IonDragGridComp.
UI_EXPORT_NAME = "UI_IONDRAG"
VI_EXPORT_NAME = "VI_IONDRAG"

# Numerical controls.
MIN_STD = 1.0e-6
MIN_PRESSURE_PA = 1.0e-30

# Ap-to-Kp conversion used by the existing GEOS-MLT ML-radiation driver.
AP_LEVELS = np.array(
    [
        0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
        39, 48, 56, 67, 80, 94, 111, 132, 154, 179,
        207, 236, 300, 400,
    ],
    dtype=np.float32,
)

KP_LEVELS = np.array(
    [
        0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0,
        4.0 / 3.0, 5.0 / 3.0, 2.0,
        7.0 / 3.0, 8.0 / 3.0, 3.0,
        10.0 / 3.0, 11.0 / 3.0, 4.0,
        13.0 / 3.0, 14.0 / 3.0, 5.0,
        16.0 / 3.0, 17.0 / 3.0, 6.0,
        19.0 / 3.0, 20.0 / 3.0, 7.0,
        22.0 / 3.0, 23.0 / 3.0, 8.0,
        25.0 / 3.0, 26.0 / 3.0, 9.0,
    ],
    dtype=np.float32,
)


torch.set_num_threads(TORCH_NUM_THREADS)
try:
    torch.set_num_interop_threads(TORCH_INTEROP_THREADS)
except RuntimeError:
    pass


# ==================
# Cached Runtime State
# ==================

_MODEL = None
_CHECKPOINT = None
_NORM = None
_INDEX_TABLE = None

_FEATURE_NAMES = None
_CONDITIONING_NAMES = None
_MODEL_CONFIG = None
_TARGET_NAMES = None

_RANK = int(
    os.environ.get(
        "OMPI_COMM_WORLD_RANK",
        os.environ.get("PMI_RANK", os.environ.get("SLURM_PROCID", "0")),
    )
)


# ==================
# I/O & Config
# ==================

def log(message):
    """Write one runtime message to stderr."""
    sys.stderr.write(str(message) + "\n")
    sys.stderr.flush()


def rank0_log(message):
    """Write a runtime message only from MPI rank zero."""
    if _RANK == 0:
        log(message)


def array_summary(array, name):
    """Return a compact finite/min/max/mean summary."""
    values = np.asarray(array)

    finite = np.isfinite(values)
    finite_count = int(np.sum(finite))

    if finite_count == 0:
        return (
            f"{name}: shape={values.shape} "
            f"finite=0/{values.size}"
        )

    finite_values = values[finite]

    return (
        f"{name}: shape={values.shape} "
        f"finite={finite_count}/{values.size} "
        f"min={float(np.min(finite_values)):.8e} "
        f"max={float(np.max(finite_values)):.8e} "
        f"mean={float(np.mean(finite_values)):.8e}"
    )


def decode_string_array(values):
    """Convert a string array to a Python list of strings."""
    output = []

    for value in np.asarray(values).reshape(-1):
        if isinstance(value, bytes):
            output.append(value.decode("utf-8"))
        else:
            output.append(str(value))

    return output


def load_normalization_file(path):
    """Load feature and target normalization data from the runtime NPZ file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Normalization file not found: {path}")

    try:
        with np.load(path, allow_pickle=False) as dataset:
            required_vars = [
                "input_var_names",
                "input_mean",
                "input_std",
                "target_var_names",
                "lev_hpa",
                "lev_log_1d",
                "y_mean_lev",
                "y_std_lev",
                "y_valid_lev",
            ]

            for name in required_vars:
                if name not in dataset:
                    raise RuntimeError(
                        f"Normalization NPZ is missing required variable: {name}"
                    )

            input_names = decode_string_array(
                dataset["input_var_names"]
            )
            input_means = np.asarray(
                dataset["input_mean"],
                dtype=np.float64,
            )
            input_stds = np.asarray(
                dataset["input_std"],
                dtype=np.float64,
            )

            target_names = decode_string_array(
                dataset["target_var_names"]
            )

            lev_hpa = np.asarray(
                dataset["lev_hpa"],
                dtype=np.float32,
            )
            lev_log = np.asarray(
                dataset["lev_log_1d"],
                dtype=np.float32,
            )

            y_mean_all = np.asarray(
                dataset["y_mean_lev"],
                dtype=np.float32,
            )
            y_std_all = np.asarray(
                dataset["y_std_lev"],
                dtype=np.float32,
            )
            y_valid_all = np.asarray(
                dataset["y_valid_lev"],
                dtype=np.int32,
            ) > 0
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            f"Failed to read ML ion-velocity normalization NPZ: {path}"
        ) from exc

    input_stats = {}

    for name, mean, std in zip(input_names, input_means, input_stds):
        if not np.isfinite(std) or std < MIN_STD:
            raise RuntimeError(
                f"Invalid normalization standard deviation for {name}: {std}"
            )

        input_stats[name] = {
            "mean": float(mean),
            "std": float(std),
        }

    required_targets = ["UI", "VI"]
    target_indices = []

    for name in required_targets:
        if name not in target_names:
            raise RuntimeError(
                f"Required target {name} is missing from normalization file. "
                f"Found: {target_names}"
            )

        target_indices.append(target_names.index(name))

    y_mean = y_mean_all[target_indices, :]
    y_std = y_std_all[target_indices, :]
    y_valid = y_valid_all[target_indices, :]

    valid_std = np.isfinite(y_std) & (y_std >= MIN_STD)
    y_valid = y_valid & valid_std & np.isfinite(y_mean)

    if not np.all(y_valid):
        bad_count = int(np.size(y_valid) - np.sum(y_valid))
        raise RuntimeError(
            f"Normalization file contains {bad_count} invalid UI/VI levels"
        )

    if not np.all(np.isfinite(lev_hpa)) or np.any(lev_hpa <= 0.0):
        raise RuntimeError(
            "Normalization pressure levels must be finite and positive"
        )

    if not np.all(np.isfinite(lev_log)):
        raise RuntimeError("Normalization log-pressure levels must be finite")

    return {
        "input_stats": input_stats,
        "target_names": required_targets,
        "lev_hpa": lev_hpa,
        "lev_pa": lev_hpa * 100.0,
        "lev_log": lev_log,
        "y_mean": y_mean,
        "y_std": y_std,
    }


def load_checkpoint(path):
    """Load and validate the trained PyTorch checkpoint."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"ML ion-velocity checkpoint not found: {path}")

    try:
        checkpoint = torch.load(
            path,
            map_location=DEVICE,
            weights_only=False,
        )
    except TypeError:
        checkpoint = torch.load(path, map_location=DEVICE)

    required_keys = [
        "model_state_dict",
        "target_vars",
        "feature_names",
        "conditioning_names",
        "model_config",
    ]

    for key in required_keys:
        if key not in checkpoint:
            raise RuntimeError(
                f"Checkpoint is missing required key: {key}"
            )

    if list(checkpoint["target_vars"]) != ["UI", "VI"]:
        raise RuntimeError(
            "Checkpoint target order must be ['UI', 'VI']; "
            f"found {checkpoint['target_vars']}"
        )

    model_config = dict(checkpoint["model_config"])
    feature_names = list(checkpoint["feature_names"])
    conditioning_names = list(checkpoint["conditioning_names"])

    expected_config = {
        "in_channels": len(feature_names),
        "cond_dim": len(conditioning_names),
        "out_channels": 2,
    }

    for key, expected_value in expected_config.items():
        if int(model_config[key]) != int(expected_value):
            raise RuntimeError(
                f"Checkpoint model_config[{key}]={model_config[key]} "
                f"but expected {expected_value}"
            )

    if "lev_log" not in feature_names:
        raise RuntimeError(
            "Checkpoint feature list does not contain lev_log"
        )

    if "lev_log" in conditioning_names:
        raise RuntimeError(
            "Checkpoint conditioning list must not contain lev_log"
        )

    return (
        checkpoint,
        feature_names,
        conditioning_names,
        model_config,
    )


# ==================
# Processing & Analysis
# ==================

def is_leap_year(year):
    """Return True when year is a Gregorian leap year."""
    return (
        year % 4 == 0
        and (year % 100 != 0 or year % 400 == 0)
    )


def ap_to_kp(ap):
    """Convert Ap to the nearest standard Kp bin."""
    index = int(
        np.argmin(
            np.abs(AP_LEVELS - np.float32(ap))
        )
    )

    return float(KP_LEVELS[index])


def load_f107_ap_table():
    """Load hourly Ap, F10.7, and F10.7A forcing."""
    global _INDEX_TABLE

    if _INDEX_TABLE is not None:
        return _INDEX_TABLE

    if not os.path.isfile(F107_AP_PATH):
        raise FileNotFoundError(
            f"F107/AP index file not found: {F107_AP_PATH}"
        )

    table = {}

    with open(F107_AP_PATH, "r", encoding="utf-8") as file_obj:
        for line in file_obj:
            parts = line.split()

            if not parts or parts[0].startswith("#"):
                continue

            try:
                if len(parts) == 6:
                    year = int(parts[0])
                    doy = int(parts[1])
                    hour = int(parts[2])
                    ap = float(parts[3])
                    f107 = float(parts[4])
                    f107a = float(parts[5])
                elif len(parts) >= 7:
                    year = int(parts[1])
                    doy = int(parts[2])
                    hour = int(parts[3])
                    ap = float(parts[4])
                    f107 = float(parts[5])
                    f107a = float(parts[6])
                else:
                    continue
            except ValueError:
                continue

            key = (year * 1000 + doy) * 24 + hour

            table[key] = {
                "ap": ap,
                "kp": ap_to_kp(ap),
                "f107": f107,
                "f107a": f107a,
            }

    if not table:
        raise RuntimeError(
            f"No valid records found in {F107_AP_PATH}"
        )

    _INDEX_TABLE = table

    return _INDEX_TABLE


def get_space_weather_indices(year, doy, decimal_hour):
    """Get the same hourly forcing interval used by the GEOS MSIS wrapper."""
    table = load_f107_ap_table()

    hour = int(float(decimal_hour))
    hour = max(0, min(23, hour))

    key = (int(year) * 1000 + int(doy)) * 24 + hour

    if key not in table:
        raise RuntimeError(
            "No F107/AP record for "
            f"year={year}, doy={doy}, hour={hour}"
        )

    return table[key]


def normalize_scalar(name, value, input_stats):
    """Apply the scalar normalization used during training."""
    if name not in input_stats:
        raise RuntimeError(
            f"Normalization statistics are missing required input: {name}"
        )

    mean = input_stats[name]["mean"]
    std = input_stats[name]["std"]

    return (float(value) - mean) / std


def build_scalar_feature_values(
    feature_names,
    year,
    doy,
    decimal_hour,
    lat_rad,
    lon_rad,
    space_weather,
    input_stats,
):
    """Build all non-vertical feature arrays in checkpoint order."""
    num_columns = lat_rad.size
    days_in_year = 366.0 if is_leap_year(year) else 365.0

    scalar_values = {}

    for harmonic in (1, 2):
        angle = (
            2.0
            * math.pi
            * harmonic
            * float(decimal_hour)
            / 24.0
        )

        scalar_values[f"ut_cos_h{harmonic}"] = np.full(
            num_columns,
            math.cos(angle),
            dtype=np.float32,
        )
        scalar_values[f"ut_sin_h{harmonic}"] = np.full(
            num_columns,
            math.sin(angle),
            dtype=np.float32,
        )

    for harmonic in (1, 2):
        angle = (
            2.0
            * math.pi
            * harmonic
            * (float(doy) - 1.0)
            / days_in_year
        )

        scalar_values[f"doy_cos_h{harmonic}"] = np.full(
            num_columns,
            math.cos(angle),
            dtype=np.float32,
        )
        scalar_values[f"doy_sin_h{harmonic}"] = np.full(
            num_columns,
            math.sin(angle),
            dtype=np.float32,
        )

    scalar_values["lat_sin"] = np.sin(lat_rad).astype(np.float32)
    scalar_values["lat_cos"] = np.cos(lat_rad).astype(np.float32)
    scalar_values["lon_sin"] = np.sin(lon_rad).astype(np.float32)
    scalar_values["lon_cos"] = np.cos(lon_rad).astype(np.float32)

    for name in ("f107", "f107a", "kp", "ap"):
        normalized = normalize_scalar(
            name,
            space_weather[name],
            input_stats,
        )

        scalar_values[name] = np.full(
            num_columns,
            normalized,
            dtype=np.float32,
        )

    supported = set(scalar_values.keys()) | {"lev_log"}

    unsupported = [
        name
        for name in feature_names
        if name not in supported
    ]

    if unsupported:
        raise RuntimeError(
            "Unsupported checkpoint input features: "
            + ", ".join(unsupported)
        )

    return scalar_values


def build_batch_inputs(
    column_indices,
    feature_names,
    conditioning_names,
    scalar_values,
    lev_log_norm,
):
    """Build one model batch using the exact checkpoint feature order."""
    num_batch = len(column_indices)
    num_levels = len(lev_log_norm)

    feature_arrays = []

    for name in feature_names:
        if name == "lev_log":
            feature = np.repeat(
                lev_log_norm[None, :],
                num_batch,
                axis=0,
            )
        else:
            values = scalar_values[name][column_indices]
            feature = np.repeat(
                values[:, None],
                num_levels,
                axis=1,
            )

        feature_arrays.append(feature.astype(np.float32))

    cond_arrays = [
        scalar_values[name][column_indices].astype(np.float32)
        for name in conditioning_names
    ]

    x = np.stack(feature_arrays, axis=1).astype(np.float32)
    cond = np.stack(cond_arrays, axis=1).astype(np.float32)

    return x, cond


def interpolate_profile_log_pressure(
    source_pressure_pa,
    source_values,
    target_pressure_pa,
):
    """Interpolate one profile in log-pressure space with endpoint clamping."""
    source_pressure_pa = np.asarray(
        source_pressure_pa,
        dtype=np.float64,
    )
    source_values = np.asarray(
        source_values,
        dtype=np.float64,
    )
    target_pressure_pa = np.asarray(
        target_pressure_pa,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(source_pressure_pa)
        & np.isfinite(source_values)
        & (source_pressure_pa > 0.0)
    )

    source_pressure = source_pressure_pa[valid]
    values = source_values[valid]

    if source_pressure.size < 2:
        if source_pressure.size == 1:
            return np.full(
                target_pressure_pa.shape,
                values[0],
                dtype=np.float32,
            )

        return np.full(
            target_pressure_pa.shape,
            np.nan,
            dtype=np.float32,
        )

    source_logp = np.log(source_pressure)
    target_logp = np.log(
        np.maximum(target_pressure_pa, MIN_PRESSURE_PA)
    )

    order = np.argsort(source_logp)
    source_logp = source_logp[order]
    values = values[order]

    output = np.interp(
        target_logp,
        source_logp,
        values,
    )

    return output.astype(np.float32)


def remap_predictions_to_geos(
    prediction_ml,
    p_train_pa,
    p_mid_geos,
):
    """Remap UI/VI predictions from the training grid to GEOS midlevels."""
    num_columns = prediction_ml.shape[0]
    num_targets = prediction_ml.shape[1]
    num_geos_levels = p_mid_geos.shape[1]

    output = np.zeros(
        (num_columns, num_targets, num_geos_levels),
        dtype=np.float32,
    )

    for column_index in range(num_columns):
        for target_index in range(num_targets):
            output[column_index, target_index, :] = (
                interpolate_profile_log_pressure(
                    source_pressure_pa=p_train_pa,
                    source_values=prediction_ml[
                        column_index,
                        target_index,
                        :,
                    ],
                    target_pressure_pa=p_mid_geos[
                        column_index,
                        :,
                    ],
                )
            )

    return output


# ==================
# Model
# ==================

class FiLMResBlock1D(nn.Module):
    """Residual 1-D convolution block with FiLM conditioning."""

    def __init__(self, channels, kernel_size):
        super().__init__()

        padding = kernel_size // 2

        self.conv1 = nn.Conv1d(
            channels,
            channels,
            kernel_size=kernel_size,
            padding=padding,
        )
        self.conv2 = nn.Conv1d(
            channels,
            channels,
            kernel_size=kernel_size,
            padding=padding,
        )
        self.activation = nn.ReLU(inplace=True)

    def forward(self, x, gamma, beta):
        output = self.conv1(x)
        output = (1.0 + gamma) * output + beta
        output = self.activation(output)
        output = self.conv2(output)
        output = output + x
        output = self.activation(output)

        return output


class ColumnFiLMCNN1D(nn.Module):
    """Column CNN used by the UI/VI training workflow."""

    def __init__(
        self,
        in_channels,
        cond_dim,
        out_channels,
        hidden_channels,
        num_blocks,
        kernel_size,
    ):
        super().__init__()

        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd")

        padding = kernel_size // 2

        self.hidden_channels = hidden_channels
        self.num_blocks = num_blocks

        self.in_conv = nn.Conv1d(
            in_channels,
            hidden_channels,
            kernel_size=kernel_size,
            padding=padding,
        )

        self.blocks = nn.ModuleList(
            [
                FiLMResBlock1D(
                    hidden_channels,
                    kernel_size,
                )
                for _ in range(num_blocks)
            ]
        )

        self.cond_mlp = nn.Sequential(
            nn.Linear(
                cond_dim,
                hidden_channels * 2,
            ),
            nn.ReLU(inplace=True),
            nn.Linear(
                hidden_channels * 2,
                hidden_channels * 2 * num_blocks,
            ),
        )

        self.out_conv = nn.Conv1d(
            hidden_channels,
            out_channels,
            kernel_size=1,
        )

        self.activation = nn.ReLU(inplace=True)

    def forward(self, x, cond):
        hidden = self.activation(
            self.in_conv(x)
        )

        gamma_beta = self.cond_mlp(cond)
        gamma_beta = gamma_beta.reshape(
            x.shape[0],
            self.num_blocks,
            2 * self.hidden_channels,
        )

        for block_index, block in enumerate(self.blocks):
            gamma, beta = gamma_beta[
                :,
                block_index,
                :,
            ].chunk(2, dim=-1)

            hidden = block(
                hidden,
                gamma.unsqueeze(-1),
                beta.unsqueeze(-1),
            )

        return self.out_conv(hidden)


def initialize_model():
    """Load the checkpoint and normalization data once per process."""
    global _MODEL
    global _CHECKPOINT
    global _NORM
    global _FEATURE_NAMES
    global _CONDITIONING_NAMES
    global _MODEL_CONFIG
    global _TARGET_NAMES

    if _MODEL is not None:
        return

    (
        checkpoint,
        feature_names,
        conditioning_names,
        model_config,
    ) = load_checkpoint(MODEL_PATH)

    norm_data = load_normalization_file(NORM_PATH)

    if norm_data["target_names"] != list(checkpoint["target_vars"]):
        raise RuntimeError(
            "Checkpoint and normalization target order do not match"
        )

    if len(norm_data["lev_hpa"]) != norm_data["y_mean"].shape[1]:
        raise RuntimeError(
            "Normalization level count does not match target statistics"
        )

    model = ColumnFiLMCNN1D(
        in_channels=int(model_config["in_channels"]),
        cond_dim=int(model_config["cond_dim"]),
        out_channels=int(model_config["out_channels"]),
        hidden_channels=int(model_config["hidden_channels"]),
        num_blocks=int(model_config["num_blocks"]),
        kernel_size=int(model_config["kernel_size"]),
    )

    model.load_state_dict(
        checkpoint["model_state_dict"]
    )
    model.to(DEVICE)
    model.eval()

    _MODEL = model
    _CHECKPOINT = checkpoint
    _NORM = norm_data
    _FEATURE_NAMES = feature_names
    _CONDITIONING_NAMES = conditioning_names
    _MODEL_CONFIG = model_config
    _TARGET_NAMES = list(checkpoint["target_vars"])

    if PRINT_MODEL_SUMMARY:
        rank0_log(
            "[MLIONVEL] loaded checkpoint: "
            f"epoch={checkpoint.get('epoch', 'unknown')} "
            f"in={model_config['in_channels']} "
            f"cond={model_config['cond_dim']} "
            f"out={model_config['out_channels']} "
            f"levels={len(norm_data['lev_hpa'])}"
        )
        rank0_log(
            "[MLIONVEL] training pressure range [hPa]: "
            f"{float(np.min(norm_data['lev_hpa'])):.8e} to "
            f"{float(np.max(norm_data['lev_hpa'])):.8e}"
        )


def predict_ion_velocity(
    ple,
    lats_rad,
    lons_rad,
    year,
    doy,
    decimal_hour,
):
    """Predict UI and VI and remap them to the current GEOS grid."""
    initialize_model()

    ple = np.asarray(ple, dtype=np.float32)
    lats_rad = np.asarray(lats_rad, dtype=np.float32)
    lons_rad = np.asarray(lons_rad, dtype=np.float32)

    if ple.ndim != 3:
        raise RuntimeError(
            f"PLE must be three-dimensional; got shape={ple.shape}"
        )

    im, jm, num_interfaces = ple.shape
    lm = num_interfaces - 1

    if lats_rad.shape != (im, jm):
        raise RuntimeError(
            f"Latitude shape mismatch: {lats_rad.shape} versus {(im, jm)}"
        )

    if lons_rad.shape != (im, jm):
        raise RuntimeError(
            f"Longitude shape mismatch: {lons_rad.shape} versus {(im, jm)}"
        )

    if not np.all(np.isfinite(lats_rad)):
        raise RuntimeError("Non-finite GEOS latitude encountered")

    if not np.all(np.isfinite(lons_rad)):
        raise RuntimeError("Non-finite GEOS longitude encountered")

    p_mid = 0.5 * (
        ple[:, :, 0:lm]
        + ple[:, :, 1:lm + 1]
    )

    p_mid = np.maximum(
        p_mid,
        MIN_PRESSURE_PA,
    ).astype(np.float32)

    num_columns = im * jm

    lat_flat = lats_rad.reshape(num_columns)
    lon_flat = lons_rad.reshape(num_columns)
    p_mid_flat = p_mid.reshape(num_columns, lm)

    space_weather = get_space_weather_indices(
        year,
        doy,
        decimal_hour,
    )

    scalar_values = build_scalar_feature_values(
        feature_names=_FEATURE_NAMES,
        year=year,
        doy=doy,
        decimal_hour=decimal_hour,
        lat_rad=lat_flat,
        lon_rad=lon_flat,
        space_weather=space_weather,
        input_stats=_NORM["input_stats"],
    )

    lev_log_norm = np.asarray(
        _NORM["lev_log"],
        dtype=np.float32,
    )

    lev_log_norm = (
        lev_log_norm
        - _NORM["input_stats"]["lev_log"]["mean"]
    ) / _NORM["input_stats"]["lev_log"]["std"]

    pred_geos = np.zeros(
        (num_columns, 2, lm),
        dtype=np.float32,
    )

    inference_context = (
        torch.inference_mode
        if hasattr(torch, "inference_mode")
        else torch.no_grad
    )

    with inference_context():
        for start in range(0, num_columns, COLUMN_BATCH_SIZE):
            end = min(
                start + COLUMN_BATCH_SIZE,
                num_columns,
            )

            column_indices = np.arange(
                start,
                end,
                dtype=np.int64,
            )

            x_np, cond_np = build_batch_inputs(
                column_indices=column_indices,
                feature_names=_FEATURE_NAMES,
                conditioning_names=_CONDITIONING_NAMES,
                scalar_values=scalar_values,
                lev_log_norm=lev_log_norm,
            )

            x_t = torch.from_numpy(x_np).to(DEVICE)
            cond_t = torch.from_numpy(cond_np).to(DEVICE)

            pred_norm = _MODEL(
                x_t,
                cond_t,
            ).cpu().numpy()

            pred_phys = (
                pred_norm
                * _NORM["y_std"][None, :, :]
                + _NORM["y_mean"][None, :, :]
            ).astype(np.float32)

            pred_geos[start:end, :, :] = (
                remap_predictions_to_geos(
                    prediction_ml=pred_phys,
                    p_train_pa=_NORM["lev_pa"],
                    p_mid_geos=p_mid_flat[start:end, :],
                )
            )

    ui = pred_geos[:, 0, :].reshape(im, jm, lm)
    vi = pred_geos[:, 1, :].reshape(im, jm, lm)

    if not np.all(np.isfinite(ui)):
        raise RuntimeError(
            "Non-finite UI values produced by ML ion-velocity model"
        )

    if not np.all(np.isfinite(vi)):
        raise RuntimeError(
            "Non-finite VI values produced by ML ion-velocity model"
        )

    return ui, vi, space_weather


# ==================
# MAPL PythonBridge
# ==================

class MLIonVelocityDriver(UserCode):
    """MAPL PythonBridge entry point for ML ion velocity."""

    def __init__(self):
        pass

    def init(
        self,
        grid_comp,
        import_state,
        export_state,
    ):
        try:
            rank0_log(
                "[MLIONVEL] Python bridge initialization"
            )
            initialize_model()
        except Exception as exc:
            log(
                "[MLIONVEL] EXCEPTION during initialization: "
                f"{repr(exc)}"
            )
            log(traceback.format_exc())
            raise

    def run(
        self,
        grid_comp,
        import_state,
        export_state,
    ):
        pass

    def run_with_internal(
        self,
        grid_comp,
        import_state,
        export_state,
        internal_state,
    ):
        try:
            initialize_model()

            mapl_py = get_MAPLPy()
            im, jm, lm = mapl_py.grid_dims

            ple = mapl_py.get_pointer(
                name="PLE",
                state=import_state,
                dims=[im, jm, lm + 1],
            )

            lats = mapl_py.get_pointer(
                name="MLION_LATS",
                state=internal_state,
                dims=[im, jm],
            )
            lons = mapl_py.get_pointer(
                name="MLION_LONS",
                state=internal_state,
                dims=[im, jm],
            )
            year_2d = mapl_py.get_pointer(
                name="MLION_YY",
                state=internal_state,
                dims=[im, jm],
            )
            doy_2d = mapl_py.get_pointer(
                name="MLION_DOY",
                state=internal_state,
                dims=[im, jm],
            )
            hour_2d = mapl_py.get_pointer(
                name="MLION_HH",
                state=internal_state,
                dims=[im, jm],
            )

            if ple is None:
                raise RuntimeError(
                    "PLE import pointer is None"
                )

            if (
                lats is None
                or lons is None
                or year_2d is None
                or doy_2d is None
                or hour_2d is None
            ):
                raise RuntimeError(
                    "One or more ML ion-velocity internal pointers are None"
                )

            year = int(year_2d[0, 0])
            doy = int(doy_2d[0, 0])
            decimal_hour = float(hour_2d[0, 0])

            ui, vi, space_weather = predict_ion_velocity(
                ple=ple,
                lats_rad=lats,
                lons_rad=lons,
                year=year,
                doy=doy,
                decimal_hour=decimal_hour,
            )

            ui_out = mapl_py.get_pointer(
                name=UI_EXPORT_NAME,
                state=export_state,
                dims=[im, jm, lm],
            )
            vi_out = mapl_py.get_pointer(
                name=VI_EXPORT_NAME,
                state=export_state,
                dims=[im, jm, lm],
            )

            if ui_out is None:
                raise RuntimeError(
                    f"Export pointer is None: {UI_EXPORT_NAME}"
                )

            if vi_out is None:
                raise RuntimeError(
                    f"Export pointer is None: {VI_EXPORT_NAME}"
                )

            ui_out[:, :, :] = ui[:, :, :]
            vi_out[:, :, :] = vi[:, :, :]

            if PRINT_OUTPUT_MINMAX and _RANK == 0:
                rank0_log(
                    array_summary(
                        ui_out,
                        "[MLIONVEL] UI_IONDRAG [m/s]",
                    )
                )
                rank0_log(
                    array_summary(
                        vi_out,
                        "[MLIONVEL] VI_IONDRAG [m/s]",
                    )
                )
                rank0_log(
                    "[MLIONVEL] forcing "
                    f"year={year} doy={doy} "
                    f"hour={decimal_hour:.3f} "
                    f"F107={space_weather['f107']:.3f} "
                    f"F107A={space_weather['f107a']:.3f} "
                    f"Kp={space_weather['kp']:.3f} "
                    f"Ap={space_weather['ap']:.3f}"
                )

        except Exception as exc:
            log(
                "[MLIONVEL] EXCEPTION in run_with_internal: "
                f"{repr(exc)}"
            )
            log(traceback.format_exc())
            raise

    def finalize(
        self,
        grid_comp,
        import_state,
        export_state,
    ):
        pass


CODE = MLIonVelocityDriver()
