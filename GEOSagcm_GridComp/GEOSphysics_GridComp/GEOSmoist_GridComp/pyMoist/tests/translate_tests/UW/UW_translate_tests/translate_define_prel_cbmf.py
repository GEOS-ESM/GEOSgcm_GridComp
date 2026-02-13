from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import calc_cumulus_base_mass_flux, define_prel_krel
from pyMoist.UW.config import UWConfiguration


class TranslateDefinePrelCbmf(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "cin": {},
            "cinlcl": {},
            "dp0": {},
            "exnifc0": {},
            "kinv": {},
            "klcl": {},
            "pifc0": {},
            "plcl": {},
            "thv0bot": {},
            "thv0lcl": {},
            "thv0top": {},
            "tkeavg": {},
            "rkfre": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "cbmf": self.grid.compute_dict(),
            "krel": self.grid.compute_dict(),
            "prel": self.grid.compute_dict(),
            "thv0rel": self.grid.compute_dict(),
            "ufrcinv": self.grid.compute_dict(),
            "wcrit": self.grid.compute_dict(),
            "winv": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        self._define_prel_krel = self.stencil_factory.from_dims_halo(
            func=define_prel_krel,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._calc_cumulus_base_mass_flux = self.stencil_factory.from_dims_halo(
            func=calc_cumulus_base_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "use_CINcin": config.use_CINcin,
                "rbuoy": config.rbuoy,
                "epsvarw": config.epsvarw,
                "dt": config.dt,
                "mumin1": config.mumin1,
                "rmaxfrac": config.rmaxfrac,
            },
        )

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        kinv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        klcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(klcl.view[:], inputs["klcl"] - 1)
        cin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(cin.view[:], inputs["cin"])
        cinlcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(cinlcl.view[:], inputs["cinlcl"])
        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        plcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(plcl.view[:], inputs["plcl"])
        thv0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0lcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0lcl.view[:], inputs["thv0lcl"])
        thv0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        tkeavg = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(tkeavg.view[:], inputs["tkeavg"])
        rkfre = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rkfre.view[:], inputs["rkfre"])

        # Outputs
        cbmf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        prel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thv0rel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ufrcinv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        wcrit = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        winv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        rho0inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        stop_cin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._define_prel_krel(
            condensation=condensation,
            iteration=iter_test,
            klcl=klcl,
            kinv=kinv,
            pifc0=pifc0,
            thv0bot=thv0bot,
            plcl=plcl,
            thv0lcl=thv0lcl,
            krel=krel,
            prel=prel,
            thv0rel=thv0rel,
        )

        self._calc_cumulus_base_mass_flux(
            condensation=condensation,
            iteration=iter_test,
            cin_IJ=cin,
            cinlcl_IJ=cinlcl,
            RKFRE=rkfre,
            tkeavg=tkeavg,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            kinv=kinv,
            pifc0=pifc0,
            thv0top=thv0top,
            exnifc0=exnifc0,
            dp0=dp0,
            winv=winv,
            cbmf=cbmf,
            rho0inv=rho0inv,
            ufrcinv=ufrcinv,
            wcrit=wcrit,
            cush_inout=cush_inout,
        )

        return {
            "cbmf": cbmf.view[:],
            "krel": krel.view[:],  # krel will fail by 1
            "prel": prel.view[:],
            "thv0rel": thv0rel.view[:],
            "ufrcinv": ufrcinv.view[:],
            "wcrit": wcrit.view[:],
            "winv": winv.view[:],
        }
