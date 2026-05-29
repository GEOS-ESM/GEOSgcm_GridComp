from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import recalc_condensate
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateRecalcCondensate(TranslateFortranData2Py):
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
            "exnifc0": {},
            "fer": {},
            "kbup": {},
            "kpen": {},
            "pifc0": {},
            "ppen": {},
            "qt0": {},
            "qtu": {},
            "ssqt0": {},
            "ssthl0": {},
            "thl0": {},
            "thlu": {},
            "thv0bot": {},
            "thv0top": {},
            "zifc0": {},
            "krel": {},
            "dwten": {},
            "diten": {},
            "umf": {},
            "emf": {},
            "ufrc": {},
            "xco": {},
            "fdr": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "cldhgt": self.grid.compute_dict(),
            "diten": self.grid.compute_dict(),
            "dwten": self.grid.compute_dict(),
            "emf": self.grid.compute_dict(),
            "fdr": self.grid.compute_dict(),
            "qtu_top": self.grid.compute_dict(),
            "thlu_top": self.grid.compute_dict(),
            "ufrc": self.grid.compute_dict(),
            "umf": self.grid.compute_dict(),
            "xco": self.grid.compute_dict(),
            "testvar3D_1": self.grid.compute_dict(),
            "testvar3D_2": self.grid.compute_dict(),
            "testvar3D_3": self.grid.compute_dict(),
            "testvar3D_4": self.grid.compute_dict(),
            "testvar3D_5": self.grid.compute_dict(),
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

        self._recalc_condensate = self.stencil_factory.from_dims_halo(
            func=recalc_condensate,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"criqc": config.criqc, "k0": config.k0},
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        exnifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        fer = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(fer.view[:], inputs["fer"])
        kbup_IJ = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup_IJ.view[:], inputs["kbup"] - 1)
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        ppen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(ppen.view[:], inputs["ppen"])
        qt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        ssqt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        ssthl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        thl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        thv0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        zifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(zifc0.view[:], inputs["zifc0"])
        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)

        # Outputs

        cldhgt = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        diten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(diten.view[:], inputs["diten"])
        dwten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dwten.view[:], inputs["dwten"])
        emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(emf.view[:], inputs["emf"])
        fdr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(fdr.view[:], inputs["fdr"])
        qtu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        thlu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(ufrc.view[:], inputs["ufrc"])
        umf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        xco = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(xco.view[:], inputs["xco"])
        kbup = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        dwten_temp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        diten_temp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        umf_temp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        fer_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
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

        testvar3D_1 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_2 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_3 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_4 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_5 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._recalc_condensate(
            condensation=condensation,
            fer=fer,
            kpen=kpen,
            ppen=ppen,
            thlu=thlu,
            thl0=thl0,
            ssthl0=ssthl0,
            qtu=qtu,
            qt0=qt0,
            ssqt0=ssqt0,
            pifc0=pifc0,
            ese=self.ese,
            esx=self.esx,
            thv0bot=thv0bot,
            thv0top=thv0top,
            exnifc0=exnifc0,
            zifc0=zifc0,
            kbup_IJ=kbup_IJ,
            kbup=kbup,
            krel=krel,
            umf_zint=umf,
            emf=emf,
            ufrc=ufrc,
            dwten=dwten,
            diten=diten,
            dwten_temp=dwten_temp,
            diten_temp=diten_temp,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            cldhgt=cldhgt,
            umf_temp=umf_temp,
            fdr=fdr,
            xco=xco,
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
            testvar3D_4=testvar3D_4,
            testvar3D_5=testvar3D_5,
        )

        return {
            "cldhgt": cldhgt.view[:],
            "diten": diten.view[:],
            "dwten": dwten.view[:],
            "emf": emf.view[:],
            "fdr": fdr.view[:],
            "qtu_top": qtu_top.view[:],
            "thlu_top": thlu_top.view[:],
            "ufrc": ufrc.view[:],
            "umf": umf.view[:],
            "xco": xco.view[:],
            "fer": fer.view[:],
            "testvar3D_1": testvar3D_1.view[:],
            "testvar3D_2": testvar3D_2.view[:],
            "testvar3D_3": testvar3D_3.view[:],
            "testvar3D_4": testvar3D_4.view[:],
            "testvar3D_5": testvar3D_5.view[:],
        }
