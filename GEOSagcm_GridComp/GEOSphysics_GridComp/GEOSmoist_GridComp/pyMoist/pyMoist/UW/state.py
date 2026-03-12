import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Float, Z_INTERFACE_DIM


@dataclasses.dataclass
class UWState(State):
    @dataclasses.dataclass
    class Input:
        """
        ComputeUwshcuInv inputs
        """

        PLE: Quantity = dataclasses.field(
            metadata={
                "name": "PLE",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "Pa",
                "intent": "?",
                "dtype": Float,
            }
        )
        ZLE: Quantity = dataclasses.field(
            metadata={
                "name": "ZLE",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )
        QLLS: Quantity = dataclasses.field(
            metadata={
                "name": "QLLS",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        QILS: Quantity = dataclasses.field(
            metadata={
                "name": "QILS",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        QLCN: Quantity = dataclasses.field(
            metadata={
                "name": "QLCN",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        QICN: Quantity = dataclasses.field(
            metadata={
                "name": "QICN",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        kpbl_inv: Quantity = dataclasses.field(
            metadata={
                "name": "kpbl_inv",
                "dims": [X_DIM, Y_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )
        frland: Quantity = dataclasses.field(
            metadata={
                "name": "frland",
                "dims": [X_DIM, Y_DIM],
                "units": "fraction",
                "intent": "?",
                "dtype": Float,
            }
        )
        tke_inv: Quantity = dataclasses.field(
            metadata={
                "name": "tke_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "m2/s2",
                "intent": "?",
                "dtype": Float,
            }
        )
        shfx: Quantity = dataclasses.field(
            metadata={
                "name": "skfx",
                "dims": [X_DIM, Y_DIM],
                "units": "J",
                "intent": "?",
                "dtype": Float,
            }
        )
        evap: Quantity = dataclasses.field(
            metadata={
                "name": "evap",
                "dims": [X_DIM, Y_DIM],
                "units": "kg/m^2/s",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Input_Output:
        """
        ComputeUwshcuInv inouts
        """

        u0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "u0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        v0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "v0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        qv0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qv0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg",
                "intent": "?",
                "dtype": Float,
            }
        )
        t0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "t0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "K",
                "intent": "?",
                "dtype": Float,
            }
        )
        cush: Quantity = dataclasses.field(
            metadata={
                "name": "cush",
                "dims": [X_DIM, Y_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )
        CNV_Tracers: Quantity = dataclasses.field(
            metadata={
                "name": "CNV_Tracers",
                "dims": [X_DIM, Y_DIM, Z_DIM, "ntracers"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        cnvtr: Quantity = dataclasses.field(
            metadata={
                "name": "cnvtr",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        CLCN: Quantity = dataclasses.field(
            metadata={
                "name": "CLCN",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Output:
        """
        ComputeUwshcuInv outputs
        """

        RKFRE: Quantity = dataclasses.field(
            metadata={
                "name": "RKFRE",
                "dims": [X_DIM, Y_DIM],
                "units": "fraction",
                "intent": "?",
                "dtype": Float,
            }
        )
        MFD_SC: Quantity = dataclasses.field(
            metadata={
                "name": "MFD_SC",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        DQADT_SC: Quantity = dataclasses.field(
            metadata={
                "name": "DQADT_SC",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        QLENT_SC: Quantity = dataclasses.field(
            metadata={
                "name": "QLENT_SC",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        QIENT_SC: Quantity = dataclasses.field(
            metadata={
                "name": "QIENT_SC",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        umf_inv: Quantity = dataclasses.field(
            metadata={
                "name": "umf_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "kg/m^2/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcm_inv: Quantity = dataclasses.field(
            metadata={
                "name": "dcm_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qtflx_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qtflx_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        slflx_inv: Quantity = dataclasses.field(
            metadata={
                "name": "slflx_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        uflx_inv: Quantity = dataclasses.field(
            metadata={
                "name": "uflx_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vflx_inv: Quantity = dataclasses.field(
            metadata={
                "name": "vflx_inv",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qvten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qvten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        qlten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qlten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        qiten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qiten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        tten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "tten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "K/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        uten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "uten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m/s^2",
                "intent": "?",
                "dtype": Float,
            }
        )
        vten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "vten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m/s^2",
                "intent": "?",
                "dtype": Float,
            }
        )
        qrten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qrten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        qsten_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qsten_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg/kg/s",
                "intent": "?",
                "dtype": Float,
            }
        )
        cufrc_inv: Quantity = dataclasses.field(
            metadata={
                "name": "cufrc_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "fraction",
                "intent": "?",
                "dtype": Float,
            }
        )
        fer_inv: Quantity = dataclasses.field(
            metadata={
                "name": "fer_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1/Pa",
                "intent": "?",
                "dtype": Float,
            }
        )
        fdr_inv: Quantity = dataclasses.field(
            metadata={
                "name": "fdr_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1/Pa",
                "intent": "?",
                "dtype": Float,
            }
        )
        ndrop_inv: Quantity = dataclasses.field(
            metadata={
                "name": "ndrop_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        nice_inv: Quantity = dataclasses.field(
            metadata={
                "name": "nice_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qldet_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qldet_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qlsub_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qlsub_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qidet_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qidet_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qisub_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qisub_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        tpert_out: Quantity = dataclasses.field(
            metadata={
                "name": "tpert_out",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qpert_out: Quantity = dataclasses.field(
            metadata={
                "name": "qpert_out",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        CNV_MFC: Quantity = dataclasses.field(
            metadata={
                "name": "CNV_MFC",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        CNV_MFD: Quantity = dataclasses.field(
            metadata={
                "name": "CNV_MFD",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        SHLW_PRC3: Quantity = dataclasses.field(
            metadata={
                "name": "SHLW_PRC3",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        SHLW_SNO3: Quantity = dataclasses.field(
            metadata={
                "name": "SHLW_SNO3",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        SC_QT: Quantity | None = dataclasses.field(
            metadata={
                "name": "SC_QT",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        SC_MSE: Quantity | None = dataclasses.field(
            metadata={
                "name": "SC_MSE",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        CUSH_SC: Quantity | None = dataclasses.field(
            metadata={
                "name": "CUSH_SC",
                "dims": [X_DIM, Y_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )
        ql0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "ql0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        qi0_inv: Quantity = dataclasses.field(
            metadata={
                "name": "qi0_inv",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    input: Input
    input_output: Input_Output
    output: Output
