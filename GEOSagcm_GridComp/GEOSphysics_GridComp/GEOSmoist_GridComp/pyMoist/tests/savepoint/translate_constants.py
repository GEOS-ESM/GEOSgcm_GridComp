from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
import pyMoist.constants as const


# class TranslateCompute(TranslateFortranData2Py):
#     def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
#         super().__init__(grid, stencil_factory)
#         self.stencil_factory = stencil_factory
#         self.quantity_factory = grid.quantity_factory
#         self._grid = grid

#         # FloatField Inputs
#         self.in_vars["data_vars"] = {const.}

#         # Float Inputs
#         self.in_vars["parameters"] = [

#         ]

#         # FloatField Outputs
#         self.out_vars = {const.
#             "QSAT": self.grid.compute_dict(),
#         }

#     def compute(self, inputs):

#         const.MAPL_PI_R8

#         return {const.
#             "QSAT": qsat.QSat.view[:],
#         }

if __name__ == "__main__":
    import xarray as xr

    ds = xr.open_dataset(
        "/home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/test_data/11.5.2/Moist/TBC_C24_L72_Debug/Constants-Out.nc"
    )

    all_public_const_module_var = [
        item for item in dir(const) if not item.startswith("__")
    ]
    unchecked_vars = set()
    for v in all_public_const_module_var:
        if v in ds.keys():
            ref = ds[v][0, 0].to_numpy()
            py = getattr(const, v)
            print(
                f"{'✅' if getattr(const, v) == ref else '❌'} {v} (ref: {ref} | py: {py})"
            )
        else:
            unchecked_vars.add(v)
    print(f"Unchecked var: {unchecked_vars}")
