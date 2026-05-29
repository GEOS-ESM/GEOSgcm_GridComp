from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import find_pbl_averages, find_pbl_height
from pyMoist.convection.UW.config import UWConfiguration


class TranslateFindPbl(TranslateFortranData2Py):
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
            "cush": {},
            "kpbl_in": {},
            "pifc0": {},
            "qt0": {},
            "thvl0": {},
            "thvl0bot": {},
            "thvl0top": {},
            "tke": {},
            "u0": {},
            "v0": {},
            "zmid0": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "kinv": self.grid.compute_dict(),
            "qtavg": self.grid.compute_dict(),
            "thvlavg": self.grid.compute_dict(),
            "thvlmin": self.grid.compute_dict(),
            "tkeavg": self.grid.compute_dict(),
            "uavg": self.grid.compute_dict(),
            "vavg": self.grid.compute_dict(),
            "tscaleh": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")
        self.constants["JASON"] = True

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": config.NCNST,
            }
        )

        self._find_pbl_height = self.stencil_factory.from_dims_halo(
            func=find_pbl_height,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "k0": config.k0,
            },
        )

        self._find_pbl_averages = self.stencil_factory.from_dims_halo(
            func=find_pbl_averages,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "qtsrchgt": config.qtsrchgt,
            },
        )

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cush.view[:, :], inputs["cush"])
        kpbl_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpbl_in.view[:, :], inputs["kpbl_in"])
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:, :, :], inputs["pifc0"])
        qt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qt0.view[:, :, :], inputs["qt0"])
        thvl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thvl0.view[:, :, :], inputs["thvl0"])
        thvl0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thvl0bot.view[:, :, :], inputs["thvl0bot"])
        thvl0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thvl0top.view[:, :, :], inputs["thvl0top"])
        tke = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(tke.view[:, :, :], inputs["tke"])
        u0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(u0.view[:, :, :], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(v0.view[:, :, :], inputs["v0"])
        zmid0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(zmid0.view[:, :, :], inputs["zmid0"])

        # Outputs
        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        qtavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvlavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvlmin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tkeavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        uavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vavg = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tscaleh = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        dcm_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qvten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        uten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fer_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._find_pbl_height(
            iteration=iter_test,
            kpbl_in=kpbl_in,
            condensation=condensation,
            kinv=kinv,
            tscaleh=tscaleh,
            cush=cush,
            cush_inout=cush_inout,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qlten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
        )

        self._find_pbl_averages(
            condensation=condensation,
            thvl0bot=thvl0bot,
            thvl0top=thvl0top,
            kinv=kinv,
            pifc0=pifc0,
            tke_in=tke,
            u0=u0,
            v0=v0,
            thvl0=thvl0,
            zmid0=zmid0,
            qt0=qt0,
            thvlmin=thvlmin,
            tkeavg=tkeavg,
            uavg=uavg,
            vavg=vavg,
            thvlavg=thvlavg,
            qtavg=qtavg,
            iteration=iter_test,
        )

        return {
            "kinv": kinv.view[:],
            "qtavg": qtavg.view[:],
            "thvlavg": thvlavg.view[:],
            "thvlmin": thvlmin.view[:],
            "tkeavg": tkeavg.view[:],
            "uavg": uavg.view[:],
            "vavg": vavg.view[:],
            "tscaleh": tscaleh.view[:],  # Its okay if tscaleh fails here, its not used until iteration = 2
        }
