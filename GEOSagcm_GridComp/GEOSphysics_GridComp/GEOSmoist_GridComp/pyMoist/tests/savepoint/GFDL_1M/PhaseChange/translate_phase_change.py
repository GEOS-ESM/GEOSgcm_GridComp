from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.PhaseChange.config import PhaseChangeConfiguration
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange


class Translatephase_change(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "estimated_inversion_strength": grid.compute_dict() | {"serialname": "EIS"},
            "p_mb": grid.compute_dict() | {"serialname": "PLmb"},
            "k_lcl": grid.compute_dict() | {"serialname": "KLCL"},
            "p_interface_mb": grid.compute_dict() | {"serialname": "PLEmb"},
            "area": grid.compute_dict() | {"serialname": "AREA"},
            "convection_fraction": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "surface_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "convective_liquid": grid.compute_dict() | {"serialname": "QLCN"},
            "convective_ice": grid.compute_dict() | {"serialname": "QICN"},
            "large_scale_liquid": grid.compute_dict() | {"serialname": "QLLS"},
            "large_scale_ice": grid.compute_dict() | {"serialname": "QILS"},
            "vapor": grid.compute_dict() | {"serialname": "Q"},
            "large_scale_cloud_fraction": grid.compute_dict() | {"serialname": "CLLS"},
            "convective_cloud_fraction": grid.compute_dict() | {"serialname": "CLCN"},
            "nactl": grid.compute_dict() | {"serialname": "NACTL"},
            "nacti": grid.compute_dict() | {"serialname": "NACTI"},
            "qsat": grid.compute_dict() | {"serialname": "QST"},
        }

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["estimated_inversion_strength"],
            self.out_vars["p_mb"],
            self.out_vars["k_lcl"],
            self.out_vars["p_interface_mb"],
            self.out_vars["area"],
            self.out_vars["convection_fraction"],
            self.out_vars["surface_type"],
            self.out_vars["nactl"],
            self.out_vars["nacti"],
            self.out_vars["qsat"],
        )
        self.out_vars.update(
            {
                "RHX": grid.compute_dict(),
                "EVAPC": grid.compute_dict(),
                "SUBLC": grid.compute_dict(),
                "RHCRIT3D": grid.compute_dict(),
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("phase_change-constants")

    def compute_from_storage(self, inputs):
        self.phase_change_config = PhaseChangeConfiguration(
            DT_MOIST=self.constants["DT_MOIST"],
            PDF_SHAPE=self.constants["PDFSHAPE"],
            CCW_EVAP_EFF=self.constants["CCW_EVAP_EFF"],
            CCI_EVAP_EFF=self.constants["CCI_EVAP_EFF"],
            TURNRHCRIT_PARAM=self.constants["TURNRHCRIT_PARAM"],
            DW_LAND=self.constants["dw_land"],
            DW_OCEAN=self.constants["dw_ocean"],
            DO_QA=self.constants["do_qa"],
            DO_MELT_FREEZE=self.constants["LMELTFRZ"],
            USE_BERGERON=bool(self.constants["USE_BERGERON"]),
        )

        phase_change = PhaseChange(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            phase_change_config=self.phase_change_config,
        )

        phase_change(
            phase_change_config=self.phase_change_config,
            **inputs,
        )

        inputs.update(
            {
                "RHX": phase_change.outputs.rhx,
                "EVAPC": phase_change.outputs.evapc,
                "SUBLC": phase_change.outputs.sublc,
                "RHCRIT3D": phase_change.outputs.rh_crit,
            }
        )
        return inputs
