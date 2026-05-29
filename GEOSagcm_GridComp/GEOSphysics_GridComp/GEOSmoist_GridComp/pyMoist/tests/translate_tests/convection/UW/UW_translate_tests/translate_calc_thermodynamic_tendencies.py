from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import calc_thermodynamic_tendencies
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class TranslateThermodynamicTendencies(TranslateFortranData2Py):
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
            "diten": {},
            "dwten": {},
            "dp0": {},
            "kpen": {},
            "qtflx": {},
            "slflx": {},
            "u0": {},
            "v0": {},
            "uf": {},
            "vf": {},
            "uflx": {},
            "vflx": {},
            "umf": {},
            "pifc0": {},
            "ppen": {},
            "prel": {},
            "qtu": {},
            "thlu": {},
            "thlu_top": {},
            "qtu_top": {},
            "emf": {},
            "ql0": {},
            "qi0": {},
            "pmid0": {},
            "kbup": {},
            "krel": {},
            "thlu_emf": {},
            "qtu_emf": {},
            "qlten_sink": {},
            "qiten_sink": {},
            "fdr": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "slten": self.grid.compute_dict(),
            "qc": self.grid.compute_dict(),
            "sten": self.grid.compute_dict(),
            "qiten": self.grid.compute_dict(),
            "qlten": self.grid.compute_dict(),
            "qvten": self.grid.compute_dict(),
            "qrten": self.grid.compute_dict(),
            "qsten": self.grid.compute_dict(),
            "testvar3D_1": self.grid.compute_dict(),
            "testvar3D_2": self.grid.compute_dict(),
            "testvar3D_3": self.grid.compute_dict(),
            "testvar3D_4": self.grid.compute_dict(),
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

        self._calc_thermodynamic_tendencies = self.stencil_factory.from_dims_halo(
            func=calc_thermodynamic_tendencies,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "frc_rasn": config.frc_rasn,
                "dt": config.dt,
            },
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        diten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(diten.view[:], inputs["diten"])
        dwten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dwten.view[:], inputs["dwten"])
        dp0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx.view[:], inputs["qtflx"])
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx.view[:], inputs["slflx"])
        u0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        uf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(uf.view[:], inputs["uf"])
        vf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(vf.view[:], inputs["vf"])
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx.view[:], inputs["uflx"])
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx.view[:], inputs["vflx"])
        umf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(emf.view[:], inputs["emf"])
        qi0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qi0.view[:], inputs["qi0"])
        ql0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ql0.view[:], inputs["ql0"])

        # Outputs
        qrten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        slten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

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

        krel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        prel = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        pifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        ppen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(ppen.view[:], inputs["ppen"])
        thlu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(thlu_top.view[:], inputs["thlu_top"])
        qtu_top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(qtu_top.view[:], inputs["qtu_top"])
        qlubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qiubelow = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qlj_2D = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qij_2D = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        kbup = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kbup.view[:], inputs["kbup"] - 1)
        fdr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(fdr.view[:], inputs["fdr"])
        pmid0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        thlu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu_emf.view[:], inputs["thlu_emf"])
        qtu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu_emf.view[:], inputs["qtu_emf"])
        qlten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qlten_sink.view[:], inputs["qlten_sink"])
        qiten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qiten_sink.view[:], inputs["qiten_sink"])
        qvten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        testvar3D_1 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_2 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_3 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        testvar3D_4 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_thermodynamic_tendencies(
            condensation=condensation,
            kpen=kpen,
            umf_zint=umf,
            dp0=dp0,
            slflx=slflx,
            uflx=uflx,
            vflx=vflx,
            qtflx=qtflx,
            u0=u0,
            v0=v0,
            uf=uf,
            vf=vf,
            dwten=dwten,
            diten=diten,
            umf_temp=umf,
            krel=krel,
            prel=prel,
            thlu=thlu,
            qtu=qtu,
            ese=self.ese,
            esx=self.esx,
            pifc0=pifc0,
            ppen=ppen,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            qlubelow=qlubelow,
            qiubelow=qiubelow,
            qlj_2D=qlj_2D,
            qij_2D=qij_2D,
            kbup=kbup,
            fdr=fdr,
            ql0=ql0,
            qi0=qi0,
            pmid0=pmid0,
            thlu_emf=thlu_emf,
            qtu_emf=qtu_emf,
            emf=emf,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            qrten=qrten,
            qsten=qsten,
            qvten=qvten,
            qlten=qlten,
            sten=sten,
            qiten=qiten,
            qc=qc,
            qlten_det=qlten_det,
            qiten_det=qiten_det,
            slten=slten,
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
            iteration=iter_test,
            testvar3D_1=testvar3D_1,
            testvar3D_2=testvar3D_2,
            testvar3D_3=testvar3D_3,
            testvar3D_4=testvar3D_4,
        )

        return {
            "slten": slten.view[:],
            "qc": qc.view[:],
            "sten": sten.view[:],
            "qiten": qiten.view[:],
            "qlten": qlten.view[:],
            "qvten": qvten.view[:],
            "qsten": qsten.view[:],
            "qrten": qrten.view[:],
            "testvar3D_1": testvar3D_1.view[:],
            "testvar3D_2": testvar3D_2.view[:],
            "testvar3D_3": testvar3D_3.view[:],
            "testvar3D_4": testvar3D_4.view[:],
        }
