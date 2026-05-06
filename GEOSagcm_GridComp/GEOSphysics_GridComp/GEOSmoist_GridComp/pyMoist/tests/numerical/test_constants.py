import os
from types import ModuleType

import xarray as xr

import pyMoist.constants as shared_const


def _check_type(type_A, type_B):
    if isinstance(type_B, bool):
        return True
    return type_A == type_B


def _get_constant_from_module(my_module: ModuleType) -> list[str]:
    # All public module var
    module_var = [item for item in dir(my_module) if not item.startswith("_")]
    # Get rid of the imports
    imports = ["np", "Float", "Int"]
    for i in imports:
        module_var.remove(i)
    # Remove non testable constants
    non_testable_const = ["MAPL_UNDEF", "NCNST", "N_MODES", "FLOAT_TINY"]
    for nc in non_testable_const:
        module_var.remove(nc)
    return module_var


def _load_refrence_nc() -> xr.Dataset:
    this_dir_path = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.abspath(os.path.join(this_dir_path, "./data"))
    print(f"Looking in {data_dir}")
    ds = xr.open_mfdataset(f"{data_dir}/Constants.*.nc")

    return ds


def _run_python_vs_fortran(fortran_value_as_dataset: xr.Dataset, my_module: ModuleType) -> None:
    unchecked_vars = set()
    failures = set()
    const_module_var = _get_constant_from_module(my_module)
    for v in const_module_var:
        if v in fortran_value_as_dataset.keys():
            # Check type
            f90 = fortran_value_as_dataset[v][0, 0].to_numpy()
            py = getattr(my_module, v)
            if hasattr(py, "dtype"):
                py_dtype = py.dtype
            else:
                py_dtype = type(py)
            if not _check_type(f90.dtype, py_dtype):
                print(f"🔴 Type differs {v} f90/py {f90.dtype} != {py_dtype} (value f90: {f90} | py: {py})")
                failures.add(v)
                continue
            # Check value
            if f90 != py:
                print(f"❌ Value differs {v} f90/py {f90.dtype} != {py_dtype}")
                failures.add(v)
                continue

            print(f"✅ {v} (f90: {f90} | py: {py})")
        else:
            unchecked_vars.add(v)

    # Fail for unchecked vars
    if unchecked_vars != set():
        assert False, f"Unchecked var: {unchecked_vars}"

    assert len(failures) == 0, f"Constants {failures} are wrong"


def test_shared_constants() -> None:
    ds = _load_refrence_nc()
    _run_python_vs_fortran(ds, shared_const)


if __name__ == "__main__":
    test_shared_constants()
