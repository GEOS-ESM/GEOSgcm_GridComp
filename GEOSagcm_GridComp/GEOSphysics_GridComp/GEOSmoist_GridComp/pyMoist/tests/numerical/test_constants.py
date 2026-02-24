import warnings
from pathlib import Path

import xarray as xr

import pyMoist.constants as const


def _dtype_is_equal(np_var, py_or_np_var):
    if isinstance(py_or_np_var, bool) or isinstance(py_or_np_var, int):
        return True
    return np_var.dtype == py_or_np_var.dtype


def _get_constant_from_module():
    # All public module var
    const_module_var = [item for item in dir(const) if not item.startswith("_")]
    # Get rid of the imports
    imports = ["np", "Float"]
    for i in imports:
        const_module_var.remove(i)
    # Remove non testable constants
    non_testable_const = ["MAPL_UNDEF", "NCNST", "N_MODES", "FLOAT_TINY"]
    for nc in non_testable_const:
        const_module_var.remove(nc)
    return const_module_var


def test_constants_match_fortran() -> None:
    this_dir_path = Path(__file__).resolve().parent
    data_dir = this_dir_path / "data"

    # Reference Fortran constants
    ds = xr.open_mfdataset(data_dir.glob("Constants.*.nc"))
    # This modules' constants
    const_module_var = _get_constant_from_module()

    unchecked_vars = set()
    for v in const_module_var:
        if v not in ds.keys():
            unchecked_vars.add(v)
            continue

        f90 = ds[v][0, 0].to_numpy()
        py = getattr(const, v)
        if _dtype_is_equal(f90, py):
            assert getattr(const, v) == f90, f"Failed to validate {v}: (f90: {f90} | py: {py})"
        else:
            warnings.warn(
                f"🔶 {v} {f90.dtype} != {py.dtype} " f"(value f90: {f90} | py: {py})",
            )

    # Fail for unchecked vars
    assert not unchecked_vars, f"Unchecked var: {unchecked_vars}"
