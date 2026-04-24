from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import setup_outputs
from pyMoist.convection.UW.config import UWConfiguration


class TranslateSetupOutputs(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "CLCN": {},
            # "CUSH": {},
            "DCM_SC": {},
            "DETR_SC": {},
            "DP": {},
            "DQIDT_SC": {},
            "DQRDT_SC": {},
            "DQSDT_SC": {},
            "DQVDT_SC": {},
            "DTDT_SC": {},
            "DUDT_SC": {},
            "DVDT_SC": {},
            "MASS": {},
            "Q": {},
            "QICN": {},
            "QIDET_SC": {},
            "QILS": {},
            "QISUB_SC": {},
            "QLCN": {},
            "QLDET_SC": {},
            "QLLS": {},
            "QLSUB_SC": {},
            "T": {},
            "U": {},
            "UMF_SC": {},
            "V": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "CLCN": self.grid.compute_dict(),
            "DQADT_SC": self.grid.compute_dict(),
            "MFD_SC": self.grid.compute_dict(),
            # "PTR3D": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "QIDET_SC": self.grid.compute_dict(),
            "QIENT_SC": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "QLDET_SC": self.grid.compute_dict(),
            "QLENT_SC": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "U": self.grid.compute_dict(),
            "V": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self._setup_outputs = self.stencil_factory.from_dims_halo(
            func=setup_outputs,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "dt": config.dt,
                "SCLM_SHALLOW": config.SCLM_SHALLOW,
            },
        )

        # Inputs
        CLCN = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(CLCN.view[:, :, :], inputs["CLCN"])
        # CUSH = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM], units="n/a")
        # safe_assign_array(CUSH.view[:, :], inputs["CUSH"])
        DCM_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DCM_SC.view[:, :, :], inputs["DCM_SC"])
        DETR_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DETR_SC.view[:, :, :], inputs["DETR_SC"])
        DP = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DP.view[:, :, :], inputs["DP"])
        DQIDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DQIDT_SC.view[:, :, :], inputs["DQIDT_SC"])
        DQRDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DQRDT_SC.view[:, :, :], inputs["DQRDT_SC"])
        DQSDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DQSDT_SC.view[:, :, :], inputs["DQSDT_SC"])
        DQVDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DQVDT_SC.view[:, :, :], inputs["DQVDT_SC"])
        DTDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DTDT_SC.view[:, :, :], inputs["DTDT_SC"])
        DUDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DUDT_SC.view[:, :, :], inputs["DUDT_SC"])
        DVDT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(DVDT_SC.view[:, :, :], inputs["DVDT_SC"])
        MASS = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(MASS.view[:, :, :], inputs["MASS"])
        Q = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(Q.view[:, :, :], inputs["Q"])
        QICN = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QICN.view[:, :, :], inputs["QICN"])
        QIDET_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QIDET_SC.view[:, :, :], inputs["QIDET_SC"])
        QILS = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QILS.view[:, :, :], inputs["QILS"])
        QISUB_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QISUB_SC.view[:, :, :], inputs["QISUB_SC"])
        QLCN = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QLCN.view[:, :, :], inputs["QLCN"])
        QLDET_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QLDET_SC.view[:, :, :], inputs["QLDET_SC"])
        QLLS = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QLLS.view[:, :, :], inputs["QLLS"])
        QLSUB_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(QLSUB_SC.view[:, :, :], inputs["QLSUB_SC"])
        T = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(T.view[:, :, :], inputs["T"])
        U = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(U.view[:, :, :], inputs["U"])
        UMF_SC = QuantityFactory.zeros(
            self.quantity_factory, dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(UMF_SC.view[:, :, :], inputs["UMF_SC"])
        V = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(V.view[:, :, :], inputs["V"])

        # Outputs
        MFD_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        DQADT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        QLENT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        QIENT_SC = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        SHLW_SNO3 = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        SHLW_PRC3 = QuantityFactory.zeros(self.quantity_factory, dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        self._setup_outputs(
            Q=Q,
            T=T,
            U=U,
            V=V,
            DQVDT_SC=DQVDT_SC,
            DTDT_SC=DTDT_SC,
            DUDT_SC=DUDT_SC,
            DVDT_SC=DVDT_SC,
            MFD_SC=MFD_SC,
            DETR_SC=DETR_SC,
            UMF_SC=UMF_SC,
            DP=DP,
            DQADT_SC=DQADT_SC,
            MASS=MASS,
            CLCN=CLCN,
            QLENT_SC=QLENT_SC,
            QIENT_SC=QIENT_SC,
            QLDET_SC=QLDET_SC,
            QIDET_SC=QIDET_SC,
            QLCN=QLCN,
            QICN=QICN,
            QLLS=QLLS,
            QILS=QILS,
            QLSUB_SC=QLSUB_SC,
            QISUB_SC=QISUB_SC,
            DQRDT_SC=DQRDT_SC,
            DQSDT_SC=DQSDT_SC,
            DQIDT_SC=DQIDT_SC,
            # CUSH=CUSH,
            SHLW_PRC3=SHLW_PRC3,
            SHLW_SNO3=SHLW_SNO3,
        )

        return {
            "CLCN": CLCN.view[:],
            "DQADT_SC": DQADT_SC.view[:],
            "MFD_SC": MFD_SC.view[:],
            "Q": Q.view[:],
            "QIDET_SC": QIDET_SC.view[:],
            "QIENT_SC": QIENT_SC.view[:],
            "QILS": QILS.view[:],
            "QLDET_SC": QLDET_SC.view[:],
            "QLENT_SC": QLENT_SC.view[:],
            "QLLS": QLLS.view[:],
            "T": T.view[:],
            "U": U.view[:],
            "V": V.view[:],
        }
