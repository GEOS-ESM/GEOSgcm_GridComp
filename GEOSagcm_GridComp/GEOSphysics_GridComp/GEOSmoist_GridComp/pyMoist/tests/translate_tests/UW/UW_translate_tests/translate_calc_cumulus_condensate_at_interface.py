from f90nml import Namelist
from gt4py.cartesian.gtscript import int32
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Int
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants

import pyMoist.constants as constants
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import calc_cumulus_condensate_at_interface
from pyMoist.UW.config import UWConfiguration


class TranslateCalcCumulusCondensate(TranslateFortranData2Py):
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
            "kpen": {},
            "krel": {},
            "pifc0": {},
            "ppen": {},
            "qcubelow": {},
            "qiubelow": {},
            "qlubelow": {},
            "qtu": {},
            "prel": {},
            "qtu_top": {},
            "rcwp": {},
            "riwp": {},
            "rlwp": {},
            "thlu": {},
            "thlu_top": {},
            "ufrc": {},
            "ufrclcl": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qcubelow": self.grid.compute_dict(),
            "qiubelow": self.grid.compute_dict(),
            "qlubelow": self.grid.compute_dict(),
            "rcwp": self.grid.compute_dict(),
            "riwp": self.grid.compute_dict(),
            "rlwp": self.grid.compute_dict(),
            "testvar3D_1": self.grid.compute_dict(),
            "testvar3D_2": self.grid.compute_dict(),
            "testvar3D_3": self.grid.compute_dict(),
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

        self._calc_cumulus_condensate_at_interfaces = self.stencil_factory.from_dims_halo(
            func=calc_cumulus_condensate_at_interface,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "criqc": config.criqc,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        prel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        prel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        ppen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        ppen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(ppen.view[:], inputs["ppen"])
        qtu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qtu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(qtu_top.view[:], inputs["qtu_top"])
        thlu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        thlu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(thlu_top.view[:], inputs["thlu_top"])
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(ufrc.view[:], inputs["ufrc"])
        ufrclcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ufrclcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ufrclcl.view[:], inputs["ufrclcl"])

        # Outputs
        qcubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qcubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(qcubelow.view[:], inputs["qcubelow"])
        qiubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qiubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(qiubelow.view[:], inputs["qiubelow"])
        qlubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qlubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(qlubelow.view[:], inputs["qlubelow"])
        rcwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        rcwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rcwp.view[:], inputs["rcwp"])
        riwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        riwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(riwp.view[:], inputs["riwp"])
        rlwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        rlwp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rlwp.view[:], inputs["rlwp"])

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qcu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qcu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        fer_out= self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        fdr_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qldet_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qidet_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        dcm_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qvten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qlten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qiten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        sten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        uten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        vten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qrten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        qsten_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        cufrc_out=self.quantity_factory.zeros(dims=[I_DIM, J_DIM,K_DIM], units="n/a")
        
        testvar3D_1 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_2 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_3 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_cumulus_condensate_at_interfaces(
            condensation=condensation,
            krel=krel,
            kpen=kpen,
            pifc0=pifc0,
            ppen=ppen,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            thlu=thlu,
            qtu=qtu,
            ese=self.ese,
            esx=self.esx,
            ufrc=ufrc,
            ufrclcl=ufrclcl,
            prel=prel,
            qcu=qcu,
            qlu=qlu,
            qiu=qiu,
            qcubelow=qcubelow,
            qlubelow=qlubelow,
            qiubelow=qiubelow,
            rcwp=rcwp,
            rlwp=rlwp,
            riwp=riwp,
            cufrc=cufrc,
            iteration=iter_test,
            cush=cush,
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
            cush_inout=cush_inout,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
            testvar3D_1=testvar3D_1,
            testvar3D_2=testvar3D_2,
            testvar3D_3=testvar3D_3,
        )

        return {
            "qcubelow": qcubelow.view[:],
            "qiubelow": qiubelow.view[:],
            "qlubelow": qlubelow.view[:],
            "rcwp": rcwp.view[:],
            "riwp": riwp.view[:],
            "rlwp": rlwp.view[:],
            "testvar3D_1": testvar3D_1.view[:],
            "testvar3D_2": testvar3D_2.view[:],
            "testvar3D_3": testvar3D_3.view[:],
        }
