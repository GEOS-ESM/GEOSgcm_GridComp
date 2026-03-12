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
from pyMoist.UW.compute_uwshcu import buoyancy_sorting
from pyMoist.UW.config import UWConfiguration


class TranslateBuoyancySorting(TranslateFortranData2Py):
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
            "dp0": {},
            "exnifc0": {},
            "exnmid0": {},
            "krel": {},
            "pifc0": {},
            "pmid0": {},
            "qt0": {},
            "qtu": {},
            "ssqt0": {},
            "ssthl0": {},
            "sstr0": {},
            "ssu0": {},
            "ssv0": {},
            "thl0": {},
            "thlu": {},
            "thv0bot": {},
            "thv0rel": {},
            "thv0top": {},
            "tscaleh": {},
            "u0": {},
            "v0": {},
            "umf": {},
            "uu": {},
            "vu": {},
            "wlcl": {},
            "wu": {},
            "zifc0": {},
            "tr0_BuoySort": {},
            "prel": {},
            "zmid0": {},
            "qsat_pe": {},
            "thvu": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "testvar3D_1": self.grid.compute_dict(),
            "testvar3D_2": self.grid.compute_dict(),
            "testvar3D_3": self.grid.compute_dict(),
            "testvar3D_4": self.grid.compute_dict(),
            "testvar3D_5": self.grid.compute_dict(),
            # "testvar3D_6": self.grid.compute_dict(),
            # "testvar3D_7": self.grid.compute_dict(),
            # "testvar3D_8": self.grid.compute_dict(),
            # "testvar3D_9": self.grid.compute_dict(),
            # "testvar3D_10": self.grid.compute_dict(),
            # "testvar3D_11": self.grid.compute_dict(),
            # "testvar3D_12": self.grid.compute_dict(),
            # "testvar3D_13": self.grid.compute_dict(),
            # "testvar3D_14": self.grid.compute_dict(),
            # "testvar3D_15": self.grid.compute_dict(),
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

        self._buoyancy_sorting = self.stencil_factory.from_dims_halo(
            func=buoyancy_sorting,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "ncnst": config.NCNST,
                "dotransport": config.dotransport,
                "k0": config.k0,
                "niter_xc": config.niter_xc,
                "criqc": config.criqc,
                "cridist_opt": config.cridist_opt,
                "rle": config.rle,
                "rbuoy": config.rbuoy,
                "mixscale": config.mixscale,
                "rkm": config.rkm,
                "detrhgt": config.detrhgt,
                "dt": config.dt,
                "rmaxfrac": config.rmaxfrac,
                "use_self_detrain": config.use_self_detrain,
                "rdrag": config.rdrag,
                "PGFc": config.PGFc,
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

        # Outputs
        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        exnifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0.view[:], inputs["exnifc0"])
        exnmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(exnmid0.view[:], inputs["exnmid0"])
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        qt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qt0.view[:], inputs["qt0"])
        qtu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        ssqt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssqt0.view[:], inputs["ssqt0"])
        ssthl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssthl0.view[:], inputs["ssthl0"])
        sstr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(sstr0.view[:], inputs["sstr0"])
        ssu0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssu0.view[:], inputs["ssu0"])
        ssv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ssv0.view[:], inputs["ssv0"])
        thl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thl0.view[:], inputs["thl0"])
        thlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        thv0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0bot.view[:], inputs["thv0bot"])
        thv0rel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0rel.view[:], inputs["thv0rel"])
        thv0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(thv0top.view[:], inputs["thv0top"])
        tscaleh = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(tscaleh.view[:], inputs["tscaleh"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        umf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        uu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uu.view[:], inputs["uu"])
        vu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vu.view[:], inputs["vu"])
        wlcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(wlcl.view[:], inputs["wlcl"])
        wu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(wu.view[:], inputs["wu"])
        zifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(zifc0.view[:], inputs["zifc0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_BuoySort"])
        prel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        qsat_pe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qsat_pe.view[:], inputs["qsat_pe"])
        zmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(zmid0.view[:], inputs["zmid0"])
        thvu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thvu.view[:], inputs["thvu"])

        dpe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        exne = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        pe = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qte = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        thle = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        thvebot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        tre = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")
        ue = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        ve = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        tru = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        emf = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thlue = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qtue = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        wue = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        wtwb = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        rei = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ufrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        drage = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        bogbot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        bogtop = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        kpen_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)
        kbup_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)
        rhomid0j = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        fer = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        dwten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        diten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        fdr = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        dcm = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        xco = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        stop_buoyancy_sort = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        testvar3D_1 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_2 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_3 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_4 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_5 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_6 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_7 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_8 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_9 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_10 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_11 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_12 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_13 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_14 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        testvar3D_15 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._buoyancy_sorting(
            condensation=condensation,
            tscaleh=tscaleh,
            krel=krel,
            wlcl=wlcl,
            prel=prel,
            pifc0=pifc0,
            thv0rel=thv0rel,
            thl0=thl0,
            ssthl0=ssthl0,
            pmid0=pmid0,
            qt0=qt0,
            ssqt0=ssqt0,
            u0=u0,
            v0=v0,
            ssu0=ssu0,
            ssv0=ssv0,
            tre=tre,
            tr0=tr0,
            sstr0=sstr0,
            thlu=thlu,
            qtu=qtu,
            wu=wu,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            qsat_pe=qsat_pe,
            zifc0=zifc0,
            zmid0=zmid0,
            thlue=thlue,
            qtue=qtue,
            wue=wue,
            wtwb=wtwb,
            dp0=dp0,
            thv0bot=thv0bot,
            exnmid0=exnmid0,
            thv0top=thv0top,
            exnifc0=exnifc0,
            tru=tru,
            emf=emf,
            thvu=thvu,
            umf_zint=umf,
            rei=rei,
            uu=uu,
            vu=vu,
            ufrc=ufrc,
            pe=pe,
            thle=thle,
            qte=qte,
            dpe=dpe,
            exne=exne,
            thvebot=thvebot,
            ue=ue,
            ve=ve,
            drage=drage,
            bogbot=bogbot,
            bogtop=bogtop,
            kpen_IJ=kpen_IJ,
            kbup_IJ=kbup_IJ,
            rhomid0j=rhomid0j,
            fer=fer,
            dwten=dwten,
            diten=diten,
            fdr=fdr,
            dcm=dcm,
            xco=xco,
            stop_buoyancy_sort=stop_buoyancy_sort,
            iteration=iter_test,
            cush_inout=cush_inout,
            testvar3D_1=testvar3D_1,
            testvar3D_2=testvar3D_2,
            testvar3D_3=testvar3D_3,
            testvar3D_4=testvar3D_4,
            testvar3D_5=testvar3D_5,
        )

        return {
            "testvar3D_1": testvar3D_1.view[:],
            "testvar3D_2": testvar3D_2.view[:],
            "testvar3D_3": testvar3D_3.view[:],
            "testvar3D_4": testvar3D_4.view[:],
            "testvar3D_5": testvar3D_5.view[:],
            # "testvar3D_6": testvar3D_6.view[:],
            # "testvar3D_7": testvar3D_7.view[:],
            # "testvar3D_8": testvar3D_8.view[:],
            # "testvar3D_9": testvar3D_9.view[:],
            # "testvar3D_10": testvar3D_10.view[:],
            # "testvar3D_11": testvar3D_11.view[:],
            # "testvar3D_12": testvar3D_12.view[:],
            # "testvar3D_13": testvar3D_13.view[:],
            # "testvar3D_14": testvar3D_14.view[:],
            # "testvar3D_15": testvar3D_15.view[:],
        }
