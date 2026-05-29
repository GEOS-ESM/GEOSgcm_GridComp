from f90nml import Namelist
from ndsl import Quantity, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import FloatField, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import _reset_mask, compute_thermodynamic_variables, compute_thv0_thvl0, compute_uwshcu_invert_before
from pyMoist.convection.UW.config import UWConfiguration
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


# Run the following at command line before running translate test
# python UW_postprocess_netcdfs.py
# ncks -A CNV_Tracers-In.nc PrepareInputs-In.nc
# ncks -A ComputeUwshcuInv-constants2.nc ComputeUwshcuInv-constants.nc


class TranslatePrepareInputs(TranslateFortranData2Py):
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
            "pifc0_inv": {},
            "zifc0_inv": {},
            "pmid0_inv": {},
            "zmid0_inv": {},
            "kpbl_inv": {},
            "exnmid0_inv": {},
            "exnifc0_inv": {},
            "dp0_inv": {},
            "u0_inv": {},
            "v0_inv": {},
            "qv0_inv": {},
            "ql0_inv": {},
            "qi0_inv": {},
            "t0_inv": {},
            "frland": {},
            "tke_inv": {},
            "rkfre": {},
            "cush": {},
            "shfx": {},
            "evap": {},
            "cnvtr": {},
            "CNV_Tracers": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qt0": self.grid.compute_dict(),
            "s0": self.grid.compute_dict(),
            "ssqt0": self.grid.compute_dict(),
            "ssthl0": self.grid.compute_dict(),
            "sstr0": self.grid.compute_dict(),
            "sstr0_o": self.grid.compute_dict(),
            "ssu0": self.grid.compute_dict(),
            "ssv0": self.grid.compute_dict(),
            "t0": self.grid.compute_dict(),
            "thl0": self.grid.compute_dict(),
            "thv0bot": self.grid.compute_dict(),
            "thv0top": self.grid.compute_dict(),
            "thvl0": self.grid.compute_dict(),
            "thvl0bot": self.grid.compute_dict(),
            "thvl0top": self.grid.compute_dict(),
            "tr0_o": self.grid.compute_dict(),
            "trflx": self.grid.compute_dict(),
            "trten": self.grid.compute_dict(),
            "tru": self.grid.compute_dict(),
            "tru_emf": self.grid.compute_dict(),
            "tke": self.grid.compute_dict(),
        }

    def make_ntracers_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [I_DIM, J_DIM, K_DIM, "ntracers"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty

    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([I_DIM, J_DIM, K_DIM], "n/a", dtype=dtype)
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def make_ij_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([I_DIM, J_DIM], "n/a", dtype=dtype)
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")
        self.constants["JASON"] = True

    def compute(self, inputs):
        self.config = UWConfiguration(**self.constants)

        self._compute_uwshcu_invert_before = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu_invert_before,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"ncnst": self.config.NCNST},
        )

        self._compute_thermodynamic_variables = self.stencil_factory.from_dims_halo(
            func=compute_thermodynamic_variables,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"ncnst": self.config.NCNST, "dotransport": self.config.dotransport},
        )

        self._compute_thv0_thvl0 = self.stencil_factory.from_dims_halo(
            func=compute_thv0_thvl0,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "ncnst": self.config.NCNST,
                "k0": self.config.k0,
                "dotransport": self.config.dotransport,
            },
        )

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": self.config.NCNST,
            }
        )

        # Field inputs
        pifc0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0_inv.view[:, :, :], inputs["pifc0_inv"])
        zifc0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(zifc0_inv.view[:, :, :], inputs["zifc0_inv"])
        pmid0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(pmid0_inv.view[:, :, :], inputs["pmid0_inv"])
        zmid0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(zmid0_inv.view[:, :, :], inputs["zmid0_inv"])
        kpbl_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(kpbl_inv.view[:, :], inputs["kpbl_inv"])
        exnmid0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(exnmid0_inv.view[:, :, :], inputs["exnmid0_inv"])
        exnifc0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(exnifc0_inv.view[:, :, :], inputs["exnifc0_inv"])
        dp0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dp0_inv.view[:, :, :], inputs["dp0_inv"])
        u0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(u0_inv.view[:, :, :], inputs["u0_inv"])
        v0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(v0_inv.view[:, :, :], inputs["v0_inv"])
        qv0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qv0_inv.view[:, :, :], inputs["qv0_inv"])
        ql0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(ql0_inv.view[:, :, :], inputs["ql0_inv"])
        qi0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qi0_inv.view[:, :, :], inputs["qi0_inv"])
        t0_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(t0_inv.view[:, :, :], inputs["t0_inv"])
        frland = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(frland.view[:, :], inputs["frland"])
        tke_inv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(tke_inv.view[:, :, :], inputs["tke_inv"])
        rkfre = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(rkfre.view[:, :], inputs["rkfre"])
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cush.view[:, :], inputs["cush"])
        shfx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(shfx.view[:, :], inputs["shfx"])
        evap = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(evap.view[:, :], inputs["evap"])

        cnvtr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cnvtr.view[:, :], inputs["cnvtr"])

        CNV_Tracers = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(CNV_Tracers.view[:, :, :, :], inputs["CNV_Tracers"])
        tr0_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")

        # FloatFieldNTracers
        sstr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        sstr0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        tr0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        trten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")

        # FloatFields
        qt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        s0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssqt0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssthl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssu0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssv0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        t0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0top = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0bot = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # _compute_uwshcu_invert_before locals
        pmid0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        u0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        v0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        zmid0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        exnmid0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        dp0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qv0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ql0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qi0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        th0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tke_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tke_flip = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        pifc0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        zifc0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        exnifc0_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        kpbl_in = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], dtype=Int, units="n/a")
        cnvtrmax = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")

        # # Call stencils
        self._compute_uwshcu_invert_before(
            # Inputs
            pmid0_inv=pmid0_inv,
            u0_inv=u0_inv,
            v0_inv=v0_inv,
            zmid0_inv=zmid0_inv,
            exnmid0_inv=exnmid0_inv,
            dp0_inv=dp0_inv,
            qv0_inv=qv0_inv,
            ql0_inv=ql0_inv,
            qi0_inv=qi0_inv,
            t0_inv=t0_inv,
            tke_inv=tke_inv,
            tke_flip=tke_flip,
            pifc0_inv=pifc0_inv,
            zifc0_inv=zifc0_inv,
            exnifc0_inv=exnifc0_inv,
            kpbl_inv=kpbl_inv,
            cnvtr=cnvtr,
            frland=frland,
            CNV_Tracers=CNV_Tracers,
            tr0_inout=tr0_inout,
            # Outputs
            pmid0_in=pmid0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            dp0_in=dp0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tke_in=tke_in,
            pifc0_in=pifc0_in,
            zifc0_in=zifc0_in,
            exnifc0_in=exnifc0_in,
            kpbl_in=kpbl_in,
            cnvtrmax=cnvtrmax,
        )

        # _compute_thermodynamics_variables locals
        u0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        v0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tr0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qv0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qi0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        pmid0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tr0_temp = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tscaleh = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        fer_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        tpert_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qpert_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        dcm_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlsub_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qisub_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ndrop_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        nice_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        self._compute_thermodynamic_variables(
            pmid0_in=pmid0_in,
            zmid0_in=zmid0_in,
            exnmid0_in=exnmid0_in,
            dp0_in=dp0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            u0=u0,
            v0=v0,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            cush_inout=cush,
            cush=cush,
            umf_out=umf_out,
            shfx=shfx,
            evap=evap,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            qt0=qt0,
            t0=t0,
            qv0=qv0,
            qi0=qi0,
            pmid0=pmid0,
            tr0=tr0,
            tr0_temp=tr0_temp,
            sstr0=sstr0,
            thl0=thl0,
            ssthl0=ssthl0,
            ssqt0=ssqt0,
            ssu0=ssu0,
            ssv0=ssv0,
            tscaleh=tscaleh,
            fer_out=fer_out,
            fdr_out=fdr_out,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
            dcm_out=dcm_out,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qlsub_out=qlsub_out,
            qisub_out=qisub_out,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            cufrc_out=cufrc_out,
        )

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        self._reset_mask = self.stencil_factory.from_dims_halo(
            func=_reset_mask,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )
        self.condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], dtype=bool, units="n/a")
        self._reset_mask(self.condensation, False)

        # _compute_thv0_thvl0 locals
        trflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        tru = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        tru_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM, "ntracers"], units="n/a")
        umf_zint = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        thlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        wu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        thvu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        thlu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qtu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vu_emf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uemf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        ql0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        uten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qcu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        qlten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        slten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qv0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ql0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qi0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        t0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        s0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        u0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        v0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qt0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thl0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssthl0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssqt0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0bot_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thv0top_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0bot_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        thvl0top_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssu0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ssv0_o = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        qvten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        uten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        dcm = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qvten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        dwten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        diten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fer = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        xco = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qc_l = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qc_i = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qtten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        ufrclcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        wlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        self._compute_thv0_thvl0(
            pmid0_in=pmid0_in,
            exnmid0_in=exnmid0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            zmid0=zmid0_in,
            pifc0_in=pifc0_in,
            ssthl0=ssthl0,
            ssqt0=ssqt0,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            u0_in=u0_in,
            v0_in=v0_in,
            condensation=self.condensation,
            ssu0=ssu0,
            ssv0=ssv0,
            tr0=tr0,
            sstr0=sstr0,
            tr0_o=tr0_o,
            sstr0_o=sstr0_o,
            trflx=trflx,
            trten=trten,
            tru=tru,
            tru_emf=tru_emf,
            umf_zint=umf_zint,
            emf=emf,
            slflx=slflx,
            qtflx=qtflx,
            uflx=uflx,
            vflx=vflx,
            thlu=thlu,
            qtu=qtu,
            uu=uu,
            vu=vu,
            wu=wu,
            thvu=thvu,
            thlu_emf=thlu_emf,
            qtu_emf=qtu_emf,
            uu_emf=uu_emf,
            vu_emf=vu_emf,
            uemf=uemf,
            thvl0bot=thvl0bot,
            thvl0top=thvl0top,
            thvl0=thvl0,
            qt0=qt0,
            t0=t0,
            qv0=qv0,
            ql0=ql0,
            qi0=qi0,
            thl0=thl0,
            thv0bot=thv0bot,
            thv0top=thv0top,
            uten=uten,
            vten=vten,
            s0=s0,
            qcu=qcu,
            qlu=qlu,
            qiu=qiu,
            cufrc=cufrc,
            ufrc=ufrc,
            qlten_det=qlten_det,
            qiten_det=qiten_det,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            sten=sten,
            slten=slten,
            qiten=qiten,
            qv0_o=qv0_o,
            ql0_o=ql0_o,
            qi0_o=qi0_o,
            t0_o=t0_o,
            s0_o=s0_o,
            u0_o=u0_o,
            v0_o=v0_o,
            qt0_o=qt0_o,
            thl0_o=thl0_o,
            thvl0_o=thvl0_o,
            ssthl0_o=ssthl0_o,
            ssqt0_o=ssqt0_o,
            thv0bot_o=thv0bot_o,
            thv0top_o=thv0top_o,
            thvl0bot_o=thvl0bot_o,
            thvl0top_o=thvl0top_o,
            ssu0_o=ssu0_o,
            ssv0_o=ssv0_o,
            cush_inout=cush_inout,
            cush=cush,
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
            fer_out=fer_out,
            fdr_out=fdr_out,
            dcm=dcm,
            qvten=qvten,
            qrten=qrten,
            qsten=qsten,
            dwten=dwten,
            diten=diten,
            fer=fer,
            fdr=fdr,
            xco=xco,
            cin=cin,
            cinlcl=cinlcl,
            cbmf=cbmf,
            qc=qc,
            qc_l=qc_l,
            qc_i=qc_i,
            qtten=qtten,
            ufrclcl=ufrclcl,
            wlcl=wlcl,
        )

        return {
            "qt0": qt0.view[:],
            "s0": s0.view[:],
            "ssqt0": ssqt0.view[:],
            "ssthl0": ssthl0.view[:],
            "sstr0": sstr0.view[:],
            "sstr0_o": sstr0_o.view[:],
            "ssu0": ssu0.view[:],
            "ssv0": ssv0.view[:],
            "t0": t0.view[:],
            "thl0": thl0.view[:],
            "thv0bot": thv0bot.view[:],
            "thv0top": thv0top.view[:],
            "thvl0": thvl0.view[:],
            "thvl0bot": thvl0bot.view[:],
            "thvl0top": thvl0top.view[:],
            "tr0_o": tr0_o.view[:],
            "trten": trten.view[:],
            "trflx": trflx.view[:],
            "tru": tru.view[:],
            "tru_emf": tru_emf.view[:],
            "tke": tke_in.view[:],
        }
