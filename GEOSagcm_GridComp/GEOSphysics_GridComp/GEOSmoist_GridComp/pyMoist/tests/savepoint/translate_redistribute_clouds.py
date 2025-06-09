from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.redistribute_clouds import RedistributeClouds
from ndsl.stencils.testing.grid import Grid


class TranslateRedistributeClouds(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = RedistributeClouds(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "cloud_fraction": grid.compute_dict() | {"serialname": "RAD_CF"},
            "convective_cloud_fraction": grid.compute_dict() | {"serialname": "CLCN"},
            "large_scale_cloud_fraction": grid.compute_dict() | {"serialname": "CLLS"},
            "liquid": grid.compute_dict() | {"serialname": "RAD_QL"},
            "convective_liquid": grid.compute_dict() | {"serialname": "QLCN"},
            "large_scale_liquid": grid.compute_dict() | {"serialname": "QLLS"},
            "ice": grid.compute_dict() | {"serialname": "RAD_QI"},
            "convective_ice": grid.compute_dict() | {"serialname": "QICN"},
            "large_scale_ice": grid.compute_dict() | {"serialname": "QILS"},
            "vapor": grid.compute_dict() | {"serialname": "RAD_QV"},
            "temperature": grid.compute_dict() | {"serialname": "T"},
        }

        # FloatField Outputs
        self.out_vars = self.in_vars["data_vars"].copy()

    # Calculated Outputs
    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return {
            "cloud_fraction": inputs["cloud_fraction"],
            "convective_cloud_fraction": inputs["convective_cloud_fraction"],
            "large_scale_cloud_fraction": inputs["large_scale_cloud_fraction"],
            "liquid": inputs["liquid"],
            "convective_liquid": inputs["convective_liquid"],
            "large_scale_liquid": inputs["large_scale_liquid"],
            "ice": inputs["ice"],
            "convective_ice": inputs["convective_ice"],
            "large_scale_ice": inputs["large_scale_ice"],
            "vapor": inputs["vapor"],
            "temperature": inputs["temperature"],
        }
