from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.driver.config import config
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver


class TranslateGFDL_1M_driver(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "RAD_QV": {},
            "RAD_QL": {},
            "RAD_QR": {},
            "RAD_QI": {},
            "RAD_QS": {},
            "RAD_QG": {},
            "RAD_CF": {},
            "NACTLI": {},
            "DQVDTmic": {},
            "DQLDTmic": {},
            "DQRDTmic": {},
            "DQIDTmic": {},
            "DQSDTmic": {},
            "DQGDTmic": {},
            "DQADTmic": {},
            "DTDTmic": {},
            "T": {},
            "W": {},
            "U": {},
            "V": {},
            "DUDTmic": {},
            "DVDTmic": {},
            "DZ": {},
            "DP": {},
            "AREA": {},
            "FRLAND": {},
            "CNV_FRC": {},
            "SRF_TYPE": {},
            "EIS": {},
            "RHCRIT3D": {},
            "REV_LS": {},
            "RSU_LS": {},
            "PRCP_RAIN": {},
            "PRCP_SNOW": {},
            "PRCP_ICE": {},
            "PRCP_GRAUPEL": {},
            "PFL_LS": {},
            "PFI_LS": {},
            "KBOT": {},
            "KTOP": {},
        }

        # Float Inputs
        self.in_vars["parameters"] = [
            "DT_MOIST",
            "ANV_ICEFALL",
            "LS_ICEFALL",
            "LHYDROSTATIC",
            "LPHYS_HYDROSTATIC",
            # Namelist options
            "mp_time",
            "t_min",
            "t_sub",
            "tau_r2g",
            "tau_smlt",
            "tau_g2r",
            "dw_land",
            "dw_ocean",
            "vi_fac",
            "vr_fac",
            "vs_fac",
            "vg_fac",
            "ql_mlt",
            "do_qa",
            "fix_negative",
            "vi_max",
            "vs_max",
            "vg_max",
            "vr_max",
            "qs_mlt",
            "qs0_crt",
            "qi_gen",
            "ql0_max",
            "qi0_max",
            "qi0_crt",
            "qr0_crt",
            "fast_sat_adj",
            "rh_inc",
            "rh_ins",
            "rh_inr",
            "const_vi",
            "const_vs",
            "const_vg",
            "const_vr",
            "use_ccn",
            "rthreshu",
            "rthreshs",
            "ccn_l",
            "ccn_o",
            "qc_crt",
            "tau_g2v",
            "tau_v2g",
            "tau_s2v",
            "tau_v2s",
            "tau_revp",
            "tau_frz",
            "do_bigg",
            "do_evap",
            "do_subl",
            "sat_adj0",
            "c_piacr",
            "tau_imlt",
            "tau_v2l",
            "tau_l2v",
            "tau_i2v",
            "tau_i2s",
            "tau_l2r",
            "qi_lim",
            "ql_gen",
            "c_paut",
            "c_psaci",
            "c_pgacs",
            "c_pgaci",
            "z_slope_liq",
            "z_slope_ice",
            "prog_ccn",
            "c_cracw",
            "alin",
            "clin",
            "preciprad",
            "cld_min",
            "use_ppm",
            "mono_prof",
            "do_sedi_heat",
            "sedi_transport",
            "do_sedi_w",
            "de_ice",
            "icloud_f",
            "irain_f",
            "mp_print",
        ]

        # FloatField Outputs
        self.out_vars = {
            "DQADTmic": self.grid.compute_dict(),
            "DTDTmic": self.grid.compute_dict(),
            "DUDTmic": self.grid.compute_dict(),
            "DVDTmic": self.grid.compute_dict(),
            "W": self.grid.compute_dict(),
            "DQVDTmic": self.grid.compute_dict(),
            "DQLDTmic": self.grid.compute_dict(),
            "DQRDTmic": self.grid.compute_dict(),
            "DQIDTmic": self.grid.compute_dict(),
            "DQSDTmic": self.grid.compute_dict(),
            "DQGDTmic": self.grid.compute_dict(),
            "PRCP_RAIN": self.grid.compute_dict(),
            "PRCP_SNOW": self.grid.compute_dict(),
            "PRCP_ICE": self.grid.compute_dict(),
            "PRCP_GRAUPEL": self.grid.compute_dict(),
            "PFL_LS": self.grid.compute_dict(),
            "PFI_LS": self.grid.compute_dict(),
            "REV_LS": self.grid.compute_dict(),
            "RSU_LS": self.grid.compute_dict(),
        }

    def make_ij_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM],
            "n/a",
        )
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def make_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def compute(self, inputs):
        # FloatField Variables
        RAD_QV = self.make_ijk_field(inputs["RAD_QV"])
        RAD_QL = self.make_ijk_field(inputs["RAD_QL"])
        RAD_QR = self.make_ijk_field(inputs["RAD_QR"])
        RAD_QI = self.make_ijk_field(inputs["RAD_QI"])
        RAD_QS = self.make_ijk_field(inputs["RAD_QS"])
        RAD_QG = self.make_ijk_field(inputs["RAD_QG"])
        RAD_CF = self.make_ijk_field(inputs["RAD_CF"])
        NACTLI = self.make_ijk_field(inputs["NACTLI"])
        DQVDTmic = self.make_ijk_field(inputs["DQVDTmic"])
        DQLDTmic = self.make_ijk_field(inputs["DQLDTmic"])
        DQRDTmic = self.make_ijk_field(inputs["DQRDTmic"])
        DQIDTmic = self.make_ijk_field(inputs["DQIDTmic"])
        DQSDTmic = self.make_ijk_field(inputs["DQSDTmic"])
        DQGDTmic = self.make_ijk_field(inputs["DQGDTmic"])
        DQADTmic = self.make_ijk_field(inputs["DQADTmic"])
        DTDTmic = self.make_ijk_field(inputs["DTDTmic"])
        T = self.make_ijk_field(inputs["T"])
        W = self.make_ijk_field(inputs["W"])
        U = self.make_ijk_field(inputs["U"])
        V = self.make_ijk_field(inputs["V"])
        DUDTmic = self.make_ijk_field(inputs["DUDTmic"])
        DVDTmic = self.make_ijk_field(inputs["DVDTmic"])
        DZ = self.make_ijk_field(inputs["DZ"])
        DP = self.make_ijk_field(inputs["DP"])
        AREA = self.make_ij_field(inputs["AREA"])
        FRLAND = self.make_ij_field(inputs["FRLAND"])
        CNV_FRC = self.make_ij_field(inputs["CNV_FRC"])
        SRF_TYPE = self.make_ij_field(inputs["SRF_TYPE"])
        EIS = self.make_ij_field(inputs["EIS"])
        RHCRIT3D = self.make_ijk_field(inputs["RHCRIT3D"])
        REV_LS = self.make_ijk_field(inputs["REV_LS"])
        RSU_LS = self.make_ijk_field(inputs["RSU_LS"])
        PRCP_RAIN = self.make_ij_field(inputs["PRCP_RAIN"])
        PRCP_SNOW = self.make_ij_field(inputs["PRCP_SNOW"])
        PRCP_ICE = self.make_ij_field(inputs["PRCP_ICE"])
        PRCP_GRAUPEL = self.make_ij_field(inputs["PRCP_GRAUPEL"])
        PFL_LS = self.make_ijk_field(inputs["PFL_LS"])
        PFI_LS = self.make_ijk_field(inputs["PFI_LS"])

        # Point Variables
        DT_MOIST = Float(inputs["DT_MOIST"])
        ANV_ICEFALL = Float(inputs["ANV_ICEFALL"])
        LS_ICEFALL = Float(inputs["LS_ICEFALL"])
        LHYDROSTATIC = bool(inputs["LHYDROSTATIC"])
        LPHYS_HYDROSTATIC = bool(inputs["LPHYS_HYDROSTATIC"])
        KBOT = Int(inputs["KBOT"])
        KTOP = Int(inputs["KTOP"])
        # Namelist options
        mp_time = Float(inputs["mp_time"])
        t_min = Float(inputs["t_min"])
        t_sub = Float(inputs["t_sub"])
        tau_r2g = Float(inputs["tau_r2g"])
        tau_smlt = Float(inputs["tau_smlt"])
        tau_g2r = Float(inputs["tau_g2r"])
        dw_land = Float(inputs["dw_land"])
        dw_ocean = Float(inputs["dw_ocean"])
        vi_fac = Float(inputs["vi_fac"])
        vr_fac = Float(inputs["vr_fac"])
        vs_fac = Float(inputs["vs_fac"])
        vg_fac = Float(inputs["vg_fac"])
        ql_mlt = Float(inputs["ql_mlt"])
        do_qa = bool(inputs["do_qa"])
        fix_negative = bool(inputs["fix_negative"])
        vi_max = Float(inputs["vi_max"])
        vs_max = Float(inputs["vs_max"])
        vg_max = Float(inputs["vg_max"])
        vr_max = Float(inputs["vr_max"])
        qs_mlt = Float(inputs["qs_mlt"])
        qs0_crt = Float(inputs["qs0_crt"])
        qi_gen = Float(inputs["qi_gen"])
        ql0_max = Float(inputs["ql0_max"])
        qi0_max = Float(inputs["qi0_max"])
        qi0_crt = Float(inputs["qi0_crt"])
        qr0_crt = Float(inputs["qr0_crt"])
        fast_sat_adj = bool(inputs["fast_sat_adj"])
        rh_inc = Float(inputs["rh_inc"])
        rh_ins = Float(inputs["rh_ins"])
        rh_inr = Float(inputs["rh_inr"])
        const_vi = bool(inputs["const_vi"])
        const_vs = bool(inputs["const_vs"])
        const_vg = bool(inputs["const_vg"])
        const_vr = bool(inputs["const_vr"])
        use_ccn = Float(inputs["use_ccn"])
        rthreshu = Float(inputs["rthreshu"])
        rthreshs = Float(inputs["rthreshs"])
        ccn_l = Float(inputs["ccn_l"])
        ccn_o = Float(inputs["ccn_o"])
        qc_crt = Float(inputs["qc_crt"])
        tau_g2v = Float(inputs["tau_g2v"])
        tau_v2g = Float(inputs["tau_v2g"])
        tau_s2v = Float(inputs["tau_s2v"])
        tau_v2s = Float(inputs["tau_v2s"])
        tau_revp = Float(inputs["tau_revp"])
        tau_frz = Float(inputs["tau_frz"])
        do_bigg = Float(inputs["do_bigg"])
        do_evap = Float(inputs["do_evap"])
        do_subl = Float(inputs["do_subl"])
        sat_adj0 = Float(inputs["sat_adj0"])
        c_piacr = Float(inputs["c_piacr"])
        tau_imlt = Float(inputs["tau_imlt"])
        tau_v2l = Float(inputs["tau_v2l"])
        tau_l2v = Float(inputs["tau_l2v"])
        tau_i2v = Float(inputs["tau_i2v"])
        tau_i2s = Float(inputs["tau_i2s"])
        tau_l2r = Float(inputs["tau_l2r"])
        qi_lim = Float(inputs["qi_lim"])
        ql_gen = Float(inputs["ql_gen"])
        c_paut = Float(inputs["c_paut"])
        c_psaci = Float(inputs["c_psaci"])
        c_pgacs = Float(inputs["c_pgacs"])
        c_pgaci = Float(inputs["c_pgaci"])
        z_slope_liq = Float(inputs["z_slope_liq"])
        z_slope_ice = Float(inputs["z_slope_ice"])
        prog_ccn = Float(inputs["prog_ccn"])
        c_cracw = Float(inputs["c_cracw"])
        alin = Float(inputs["alin"])
        clin = Float(inputs["clin"])
        preciprad = Float(inputs["preciprad"])
        cld_min = Float(inputs["cld_min"])
        use_ppm = bool(inputs["use_ppm"])
        mono_prof = Float(inputs["mono_prof"])
        do_sedi_heat = Float(inputs["do_sedi_heat"])
        sedi_transport = bool(inputs["sedi_transport"])
        do_sedi_w = bool(inputs["do_sedi_w"])
        de_ice = Float(inputs["de_ice"])
        icloud_f = Float(inputs["icloud_f"])
        irain_f = Float(inputs["irain_f"])
        mp_print = Float(inputs["mp_print"])

        # Create config class, to be replaced by a proper feature w/pyMoist integration
        namelist = config(
            LPHYS_HYDROSTATIC,
            LHYDROSTATIC,
            DT_MOIST,
            mp_time,
            t_min,
            t_sub,
            tau_r2g,
            tau_smlt,
            tau_g2r,
            dw_land,
            dw_ocean,
            vi_fac,
            vr_fac,
            vs_fac,
            vg_fac,
            ql_mlt,
            do_qa,
            fix_negative,
            vi_max,
            vs_max,
            vg_max,
            vr_max,
            qs_mlt,
            qs0_crt,
            qi_gen,
            ql0_max,
            qi0_max,
            qi0_crt,
            qr0_crt,
            fast_sat_adj,
            rh_inc,
            rh_ins,
            rh_inr,
            const_vi,
            const_vs,
            const_vg,
            const_vr,
            use_ccn,
            rthreshu,
            rthreshs,
            ccn_l,
            ccn_o,
            qc_crt,
            tau_g2v,
            tau_v2g,
            tau_s2v,
            tau_v2s,
            tau_revp,
            tau_frz,
            do_bigg,
            do_evap,
            do_subl,
            sat_adj0,
            c_piacr,
            tau_imlt,
            tau_v2l,
            tau_l2v,
            tau_i2v,
            tau_i2s,
            tau_l2r,
            qi_lim,
            ql_gen,
            c_paut,
            c_psaci,
            c_pgacs,
            c_pgaci,
            z_slope_liq,
            z_slope_ice,
            prog_ccn,
            c_cracw,
            alin,
            clin,
            preciprad,
            cld_min,
            use_ppm,
            mono_prof,
            do_sedi_heat,
            sedi_transport,
            do_sedi_w,
            de_ice,
            icloud_f,
            irain_f,
            mp_print,
        )

        stencil = MicrophysicsDriver(
            self.stencil_factory,
            self.quantity_factory,
            namelist,
        )

        stencil(
            namelist,
            RAD_QV,
            RAD_QL,
            RAD_QR,
            RAD_QI,
            RAD_QS,
            RAD_QG,
            RAD_CF,
            NACTLI,
            DQVDTmic,
            DQLDTmic,
            DQRDTmic,
            DQIDTmic,
            DQSDTmic,
            DQGDTmic,
            DQADTmic,
            DTDTmic,
            T,
            W,
            U,
            V,
            DUDTmic,
            DVDTmic,
            DZ,
            DP,
            AREA,
            FRLAND,
            CNV_FRC,
            SRF_TYPE,
            EIS,
            RHCRIT3D,
            ANV_ICEFALL,
            LS_ICEFALL,
        )

        return {
            "DQADTmic": DQADTmic.view[:],
            "DTDTmic": DTDTmic.view[:],
            "DUDTmic": DUDTmic.view[:],
            "DVDTmic": DVDTmic.view[:],
            "W": W.view[:],
            "DQVDTmic": DQVDTmic.view[:],
            "DQLDTmic": DQLDTmic.view[:],
            "DQRDTmic": DQRDTmic.view[:],
            "DQIDTmic": DQIDTmic.view[:],
            "DQSDTmic": DQSDTmic.view[:],
            "DQGDTmic": DQGDTmic.view[:],
            "PRCP_RAIN": stencil.outputs.rain.view[:],
            "PRCP_SNOW": stencil.outputs.snow.view[:],
            "PRCP_ICE": stencil.outputs.ice.view[:],
            "PRCP_GRAUPEL": stencil.outputs.graupel.view[:],
            "PFL_LS": stencil.outputs.m2_rain.view[:],
            "PFI_LS": stencil.outputs.m2_sol.view[:],
            "REV_LS": stencil.outputs.revap.view[:],
            "RSU_LS": stencil.outputs.isubl.view[:],
        }
