import pyMoist.constants as const

import numpy as np


def _check_dtype(np_var, py_or_np_var):
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


if __name__ == "__main__":
    import xarray as xr

    test_data_dir = "/home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/test_data/11.5.2/Moist/TBC_C24_L72_Debug"

    ds = xr.open_mfdataset(f"{test_data_dir}/Constants.*.nc")
    const_module_var = _get_constant_from_module()
    unchecked_vars = set()
    for v in const_module_var:
        if v in ds.keys():
            f90 = ds[v][0, 0].to_numpy()
            py = getattr(const, v)
            if not _check_dtype(f90, py):
                print(
                    f"🔶 {v} {f90.dtype} != {py.dtype} "
                    f"(value f90: {f90} | py: {py})"
                )
            else:
                print(
                    f"{'✅' if getattr(const, v) == f90 else '❌'} {v} "
                    f"(f90: {f90} | py: {py})"
                )
        else:
            unchecked_vars.add(v)
    print(f"Unchecked var: {unchecked_vars}")
