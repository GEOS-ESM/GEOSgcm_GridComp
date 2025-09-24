from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.driver_interface.driver_interface import DriverInterface
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios
from pyMoist.convective_parameterization.GF_2020.GF_2020 import GF2020

# required only for workarounds
from ndsl.dsl.typing import Int
import numpy as np


class TranslateGF_2020_driver_interface_setup(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self.quantity_factory.set_extra_dim_lengths(
            **{
                "nmp": 2,
                "maxiens": 3,
            }
        )

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {
            "PLE": grid.compute_dict(),
            "PLO": grid.compute_dict(),
            "ZLE": grid.compute_dict(),
            "ZLO": grid.compute_dict(),
            "PK": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "KH": grid.compute_dict(),
            "T1": grid.compute_dict(),
            "Q1": grid.compute_dict(),
            "U1": grid.compute_dict(),
            "V1": grid.compute_dict(),
            "W1": grid.compute_dict(),
            "BYNCY": grid.compute_dict(),
            "QLCN": grid.compute_dict(),
            "QICN": grid.compute_dict(),
            "QLLS": grid.compute_dict(),
            "QILS": grid.compute_dict(),
            "CLCN": grid.compute_dict(),
            "CLLS": grid.compute_dict(),
            "QV_DYN_IN": grid.compute_dict(),
            "PLE_DYN_IN": grid.compute_dict(),
            "U_DYN_IN": grid.compute_dict(),
            "V_DYN_IN": grid.compute_dict(),
            "T_DYN_IN": grid.compute_dict(),
            "RADSW": grid.compute_dict(),
            "RADLW": grid.compute_dict(),
            "DQDT_BL": grid.compute_dict(),
            "DTDT_BL": grid.compute_dict(),
            "FRLAND": grid.compute_dict(),
            "AREA": grid.compute_dict(),
            "T2M": grid.compute_dict(),
            "SH": grid.compute_dict(),
            "EVAP": grid.compute_dict(),
            "PHIS": grid.compute_dict(),
            "KPBLIN": grid.compute_dict(),
            "DTDTDYN": grid.compute_dict(),
            "DQVDTDYN": grid.compute_dict(),
            "TPWI": grid.compute_dict(),
            "TPWI_star": grid.compute_dict(),
            "CNV_TR": grid.compute_dict(),
            "WQT_DC": grid.compute_dict(),
            "CNV_MFC": grid.compute_dict(),
            "PRFIL": grid.compute_dict(),
            "CNV_MF0": grid.compute_dict(),
            "CNV_PRC3": grid.compute_dict(),
            "CNV_MFD": grid.compute_dict(),
            "CNV_DQCDT": grid.compute_dict(),
            "CNV_UPDF": grid.compute_dict(),
            "CNV_CVW": grid.compute_dict(),
            "CNV_QC": grid.compute_dict(),
            "ENTLAM": grid.compute_dict(),
            "REVSU": grid.compute_dict(),
            "DQDT_GF": grid.compute_dict(),
            "DTDT_GF": grid.compute_dict(),
            "DUDT_GF": grid.compute_dict(),
            "DVDT_GF": grid.compute_dict(),
            "MUPDP": grid.compute_dict(),
            "MDNDP": grid.compute_dict(),
            "MUPSH": grid.compute_dict(),
            "MUPMD": grid.compute_dict(),
            "CNPCPRATE": grid.compute_dict(),
            "LIGHTN_DENS": grid.compute_dict(),
            "SIGMA_DEEP": grid.compute_dict(),
            "SIGMA_MID": grid.compute_dict(),
            "MFDP": grid.compute_dict(),
            "MFSH": grid.compute_dict(),
            "MFMD": grid.compute_dict(),
            "ERRDP": grid.compute_dict(),
            "ERRSH": grid.compute_dict(),
            "ERRMD": grid.compute_dict(),
            "AA0": grid.compute_dict(),
            "AA1": grid.compute_dict(),
            "AA2": grid.compute_dict(),
            "AA3": grid.compute_dict(),
            "AA1_BL": grid.compute_dict(),
            "AA1_CIN": grid.compute_dict(),
            "TAU_BL": grid.compute_dict(),
            "TAU_EC": grid.compute_dict(),
        }

        self.out_vars = {
            "aot500": grid.compute_dict(),
            "temp2m": grid.compute_dict(),
            "sflux_r": grid.compute_dict(),
            "sflux_t": grid.compute_dict(),
            "topt": grid.compute_dict(),
            "xland": grid.compute_dict(),
            "dx2d": grid.compute_dict(),
            "kpbl": grid.compute_dict(),
            "temp": grid.compute_dict(),
            "press": grid.compute_dict(),
            "rvap": grid.compute_dict(),
            "up": grid.compute_dict(),
            "vp": grid.compute_dict(),
            "wp": grid.compute_dict(),
            "zt3d": grid.compute_dict(),
            "zm3d": grid.compute_dict(),
            "dm3d": grid.compute_dict(),
            "khloc": grid.compute_dict(),
            "curr_rvap": grid.compute_dict(),
            "mp_ice": grid.compute_dict(),
            "mp_liq": grid.compute_dict(),
            "mp_cf": grid.compute_dict(),
            "buoy_exc": grid.compute_dict(),
            "DZ": grid.compute_dict(),
            "AIR_DEN": grid.compute_dict(),
            "entr3d": grid.compute_dict(),
            # debug
            "MASS_N": grid.compute_dict(),
            "ZLE_N": grid.compute_dict(),
            "ZLO_N": grid.compute_dict(),
            "PLO_N": grid.compute_dict(),
            "PK_N": grid.compute_dict(),
            "ec3d": grid.compute_dict(),
        }

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        self.temporaries = Temporaries.make(self.quantity_factory)

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface is True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    def make_ij_quantity(self, data) -> Quantity:
        quantity = self.quantity_factory.empty([X_DIM, Y_DIM], "n/a")
        quantity.view[:, :] = quantity.np.asarray(data[:, :])
        return quantity

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF_2020-constants")

    def compute(self, inputs):
        GF_2020_config = GF2020Config(
            STOCHASTIC_CONVECTION=bool(self.constants["STOCHASTIC_CNV"]),
            STOCH_TOP=self.constants["STOCH_TOP"],
            STOCH_BOT=self.constants["STOCH_BOT"],
            GF_MIN_AREA=self.constants["GF_MIN_AREA"],
            GF_ENV_SETTING=int(self.constants["GF_ENV_SETTING"]),
            ENTRVERSION=self.constants["ENTRVERSION"],
            CONVECTION_TRACER=self.constants["CONVECTION_TRACER"],
            C1=self.constants["C1"],
            ADV_TRIGGER=self.constants["ADV_TRIGGER"],
        )

        PLE = self.make_ijk_quantity(inputs.pop("PLE"), interface=True)
        PLO = self.make_ijk_quantity(inputs.pop("PLO"))
        ZLE = self.make_ijk_quantity(inputs.pop("ZLE"), interface=True)
        ZLO = self.make_ijk_quantity(inputs.pop("ZLO"))
        PK = self.make_ijk_quantity(inputs.pop("PK"))
        MASS = self.make_ijk_quantity(inputs.pop("MASS"))
        KH = self.make_ijk_quantity(inputs.pop("KH"), interface=True)
        T1 = self.make_ijk_quantity(inputs.pop("T1"))
        Q1 = self.make_ijk_quantity(inputs.pop("Q1"))
        U1 = self.make_ijk_quantity(inputs.pop("U1"))
        V1 = self.make_ijk_quantity(inputs.pop("V1"))
        W1 = self.make_ijk_quantity(inputs.pop("W1"))
        BYNCY = self.make_ijk_quantity(inputs.pop("BYNCY"))
        QLCN = self.make_ijk_quantity(inputs.pop("QLCN"))
        QICN = self.make_ijk_quantity(inputs.pop("QICN"))
        QLLS = self.make_ijk_quantity(inputs.pop("QLLS"))
        QILS = self.make_ijk_quantity(inputs.pop("QILS"))
        CLCN = self.make_ijk_quantity(inputs.pop("CLCN"))
        CLLS = self.make_ijk_quantity(inputs.pop("CLLS"))
        QV_DYN_IN = self.make_ijk_quantity(inputs.pop("QV_DYN_IN"))
        PLE_DYN_IN = self.make_ijk_quantity(inputs.pop("PLE_DYN_IN"), interface=True)
        U_DYN_IN = self.make_ijk_quantity(inputs.pop("U_DYN_IN"))
        V_DYN_IN = self.make_ijk_quantity(inputs.pop("V_DYN_IN"))
        T_DYN_IN = self.make_ijk_quantity(inputs.pop("T_DYN_IN"))
        RADSW = self.make_ijk_quantity(inputs.pop("RADSW"))
        RADLW = self.make_ijk_quantity(inputs.pop("RADLW"))
        DQDT_BL = self.make_ijk_quantity(inputs.pop("DQDT_BL"))
        DTDT_BL = self.make_ijk_quantity(inputs.pop("DTDT_BL"))
        FRLAND = self.make_ij_quantity(inputs.pop("FRLAND"))
        AREA = self.make_ij_quantity(inputs.pop("AREA"))
        T2M = self.make_ij_quantity(inputs.pop("T2M"))
        SH = self.make_ij_quantity(inputs.pop("SH"))
        EVAP = self.make_ij_quantity(inputs.pop("EVAP"))
        PHIS = self.make_ij_quantity(inputs.pop("PHIS"))
        KPBLIN = self.make_ij_quantity(inputs.pop("KPBLIN"))
        DTDTDYN = self.make_ijk_quantity(inputs.pop("DTDTDYN"))
        DQVDTDYN = self.make_ijk_quantity(inputs.pop("DQVDTDYN"))
        TPWI = self.make_ij_quantity(inputs.pop("TPWI"))
        TPWI_star = self.make_ij_quantity(inputs.pop("TPWI_star"))
        CNV_TR = self.make_ijk_quantity(inputs.pop("CNV_TR"))

        # inputs just getting reset to zero (not explicitly tested)
        WQT_DC = self.make_ijk_quantity(inputs.pop("WQT_DC"), interface=True)
        CNV_MFC = self.make_ijk_quantity(inputs.pop("CNV_MFC"), interface=True)
        PRFIL = self.make_ijk_quantity(inputs.pop("PRFIL"), interface=True)
        CNV_MF0 = self.make_ijk_quantity(inputs.pop("CNV_MF0"))
        CNV_PRC3 = self.make_ijk_quantity(inputs.pop("CNV_PRC3"))
        CNV_MFD = self.make_ijk_quantity(inputs.pop("CNV_MFD"))
        CNV_DQCDT = self.make_ijk_quantity(inputs.pop("CNV_DQCDT"))
        CNV_UPDF = self.make_ijk_quantity(inputs.pop("CNV_UPDF"))
        CNV_CVW = self.make_ijk_quantity(inputs.pop("CNV_CVW"))
        CNV_QC = self.make_ijk_quantity(inputs.pop("CNV_QC"))
        ENTLAM = self.make_ijk_quantity(inputs.pop("ENTLAM"))
        REVSU = self.make_ijk_quantity(inputs.pop("REVSU"))
        DQDT_GF = self.make_ijk_quantity(inputs.pop("DQDT_GF"))
        DTDT_GF = self.make_ijk_quantity(inputs.pop("DTDT_GF"))
        DUDT_GF = self.make_ijk_quantity(inputs.pop("DUDT_GF"))
        DVDT_GF = self.make_ijk_quantity(inputs.pop("DVDT_GF"))
        MUPDP = self.make_ijk_quantity(inputs.pop("MUPDP"))
        MDNDP = self.make_ijk_quantity(inputs.pop("MDNDP"))
        MUPSH = self.make_ijk_quantity(inputs.pop("MUPSH"))
        MUPMD = self.make_ijk_quantity(inputs.pop("MUPMD"))
        CNPCPRATE = self.make_ij_quantity(inputs.pop("CNPCPRATE"))
        LIGHTN_DENS = self.make_ij_quantity(inputs.pop("LIGHTN_DENS"))
        SIGMA_DEEP = self.make_ij_quantity(inputs.pop("SIGMA_DEEP"))
        SIGMA_MID = self.make_ij_quantity(inputs.pop("SIGMA_MID"))
        MFDP = self.make_ij_quantity(inputs.pop("MFDP"))
        MFSH = self.make_ij_quantity(inputs.pop("MFSH"))
        MFMD = self.make_ij_quantity(inputs.pop("MFMD"))
        ERRDP = self.make_ij_quantity(inputs.pop("ERRDP"))
        ERRSH = self.make_ij_quantity(inputs.pop("ERRSH"))
        ERRMD = self.make_ij_quantity(inputs.pop("ERRMD"))
        AA0 = self.make_ij_quantity(inputs.pop("AA0"))
        AA1 = self.make_ij_quantity(inputs.pop("AA1"))
        AA2 = self.make_ij_quantity(inputs.pop("AA2"))
        AA3 = self.make_ij_quantity(inputs.pop("AA3"))
        AA1_BL = self.make_ij_quantity(inputs.pop("AA1_BL"))
        AA1_CIN = self.make_ij_quantity(inputs.pop("AA1_CIN"))
        TAU_BL = self.make_ij_quantity(inputs.pop("TAU_BL"))
        TAU_EC = self.make_ij_quantity(inputs.pop("TAU_EC"))

        # Construct stencils
        driver_interface = DriverInterface(
            stencil_factory=self.stencil_factory,
            GF_2020_config=GF_2020_config,
        )

        maximum_t2m = np.max(T2M.field)

        # debug stuff
        MASS_N = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ZLE_N = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        ZLE_N_surface = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ZLO_N = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PLO_N = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PK_N = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ec3d = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        driver_interface(
            PLE,
            PLO,
            ZLE,
            ZLO,
            PK,
            MASS,
            KH,
            T1,
            Q1,
            U1,
            V1,
            W1,
            BYNCY,
            QLCN,
            QICN,
            QLLS,
            QILS,
            CLCN,
            CLLS,
            QV_DYN_IN,
            PLE_DYN_IN,
            U_DYN_IN,
            V_DYN_IN,
            T_DYN_IN,
            RADSW,
            RADLW,
            DQDT_BL,
            DTDT_BL,
            FRLAND,
            AREA,
            T2M,
            SH,
            EVAP,
            PHIS,
            KPBLIN,
            DTDTDYN,
            DQVDTDYN,
            TPWI,
            TPWI_star,
            CNV_TR,
            # fields computed here and passed to the driver
            self.temporaries.aot500,
            self.temporaries.temp2m,
            self.temporaries.sflux_r,
            self.temporaries.sflux_t,
            self.temporaries.topt,
            self.temporaries.xland,
            self.temporaries.dx2d,
            self.temporaries.kpbl,
            self.temporaries.temp,
            self.temporaries.press,
            self.temporaries.rvap,
            self.temporaries.up,
            self.temporaries.vp,
            self.temporaries.wp,
            self.temporaries.zt3d,
            self.temporaries.zm3d,
            self.temporaries.dm3d,
            self.temporaries.khloc,
            self.temporaries.curr_rvap,
            self.temporaries.mp_ice,
            self.temporaries.mp_liq,
            self.temporaries.mp_cf,
            self.temporaries.buoy_exc,
            # fields computed here and used immediately after driver conclusion
            self.temporaries.DZ,
            self.temporaries.AIR_DEN,
            # fields computed here and used elsewhere in the model
            self.temporaries.entr3d,
            # fields just getting reset to zero
            WQT_DC,
            CNV_MFC,
            PRFIL,
            CNV_MF0,
            CNV_PRC3,
            CNV_MFD,
            CNV_DQCDT,
            CNV_UPDF,
            CNV_CVW,
            CNV_QC,
            ENTLAM,
            REVSU,
            DQDT_GF,
            DTDT_GF,
            DUDT_GF,
            DVDT_GF,
            MUPDP,
            MDNDP,
            MUPSH,
            MUPMD,
            CNPCPRATE,
            LIGHTN_DENS,
            SIGMA_DEEP,
            SIGMA_MID,
            MFDP,
            MFSH,
            MFMD,
            ERRDP,
            ERRSH,
            ERRMD,
            AA0,
            AA1,
            AA2,
            AA3,
            AA1_BL,
            AA1_CIN,
            TAU_BL,
            TAU_EC,
            # workarounds
            maximum_t2m,
            Int(PLE.field.shape[0]),
            Int(PLE.field.shape[1]),
            # debug stuff
            MASS_N,
            ZLE_N,
            ZLE_N_surface,
            ZLO_N,
            PLO_N,
            PK_N,
            ec3d,
        )

        ZLE_N.field[:, :, -1] = ZLE_N_surface.field

        return {
            "aot500": self.temporaries.aot500.field,
            "temp2m": self.temporaries.temp2m.field,
            "sflux_r": self.temporaries.sflux_r.field,
            "sflux_t": self.temporaries.sflux_t.field,
            "topt": self.temporaries.topt.field,
            "xland": self.temporaries.xland.field,
            "dx2d": self.temporaries.dx2d.field,
            "kpbl": self.temporaries.kpbl.field,
            "temp": np.moveaxis(self.temporaries.temp.field, 2, 0),
            "press": np.moveaxis(self.temporaries.press.field, 2, 0),
            "rvap": np.moveaxis(self.temporaries.rvap.field, 2, 0),
            "up": np.moveaxis(self.temporaries.up.field, 2, 0),
            "vp": np.moveaxis(self.temporaries.vp.field, 2, 0),
            "wp": np.moveaxis(self.temporaries.wp.field, 2, 0),
            "zt3d": np.moveaxis(self.temporaries.zt3d.field, 2, 0),
            "zm3d": np.moveaxis(self.temporaries.zm3d.field, 2, 0),
            "dm3d": np.moveaxis(self.temporaries.dm3d.field, 2, 0),
            "khloc": np.moveaxis(self.temporaries.khloc.field, 2, 0),
            "curr_rvap": np.moveaxis(self.temporaries.curr_rvap.field, 2, 0),
            "mp_ice": np.transpose(self.temporaries.mp_ice.field, (3, 2, 0, 1)),
            "mp_liq": np.transpose(self.temporaries.mp_liq.field, (3, 2, 0, 1)),
            "mp_cf": np.transpose(self.temporaries.mp_cf.field, (3, 2, 0, 1)),
            "buoy_exc": np.moveaxis(self.temporaries.buoy_exc.field, 2, 0),
            "DZ": self.temporaries.DZ.field,
            "AIR_DEN": self.temporaries.AIR_DEN.field,
            "entr3d": self.temporaries.entr3d.field,
            "MASS_N": MASS_N.field,
            "ZLE_N": ZLE_N.field,
            "ZLO_N": ZLO_N.field,
            "PLO_N": PLO_N.field,
            "PK_N": PK_N.field,
            "ec3d": np.moveaxis(ec3d.field, 2, 0),
        }
