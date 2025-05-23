from __future__ import annotations

import os
from typing import Dict, List

import cffi
import numpy as np
from mpi4py import MPI

from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import Float
from ndsl.optional_imports import cupy as cp
from pyMoist.interface.cuda_profiler import CUDAProfiler, TimedCUDAProfiler
from pyMoist.interface.f_py_conversion import FortranPythonConversion
from pyMoist.interface.flags import gfdl_1m_flags_f_to_python, moist_flags_f_to_python
from pyMoist.interface.wrapper import GEOSPyMoistWrapper, MemorySpace


class PYMOIST_WRAPPER:
    def __init__(self) -> None:
        self.ready = False

    def init(
        self,
        pyMoist_flags: cffi.FFI.CData,
        backend: str = "dace:cpu",
    ) -> None:
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.backend = backend
        self.flags = moist_flags_f_to_python(pyMoist_flags)
        print(f"Moist Flags:\n{self.flags}")
        # For Fortran<->NumPy conversion
        if is_gpu_backend(self.backend):
            numpy_module = cp
            self.fortran_mem_space = MemorySpace.DEVICE
        else:
            numpy_module = np
            self.fortran_mem_space = MemorySpace.HOST
        self.f_py = FortranPythonConversion(
            self.flags.npx,
            self.flags.npy,
            self.flags.npz,
            numpy_module,
        )

        # Initalize pyMoist
        self.pymoist = GEOSPyMoistWrapper(self.flags, backend)

        self._timings: Dict[str, List[float]] = {}
        self.ready = True

    def finalize(self):
        import json

        with open("pymoist_timings.json", "w") as f:
            json.dump(self._timings, f, indent=4)

    def aer_activation(
        self,
        f_aero_dgn: cffi.FFI.CData,
        f_aero_num: cffi.FFI.CData,
        f_aero_hygroscopicity: cffi.FFI.CData,
        f_aero_sigma: cffi.FFI.CData,
        f_frland: cffi.FFI.CData,
        f_nn_ocean: np.float32,
        f_nn_land: np.float32,
        f_t: cffi.FFI.CData,
        f_plo: cffi.FFI.CData,
        f_qicn: cffi.FFI.CData,
        f_qils: cffi.FFI.CData,
        f_qlcn: cffi.FFI.CData,
        f_qlls: cffi.FFI.CData,
        f_vvel: cffi.FFI.CData,
        f_tke: cffi.FFI.CData,
        f_nacti: cffi.FFI.CData,
        f_nwfa: cffi.FFI.CData,
        f_nactl: cffi.FFI.CData,
    ):
        CUDAProfiler.start_cuda_profiler()
        with TimedCUDAProfiler("[AER] Fortran -> Python", self._timings):
            aero_dgn = self.f_py.fortran_to_python(
                f_aero_dgn,
                [
                    self.flags.npx,
                    self.flags.npy,
                    self.flags.npz,
                    self.flags.n_modes,
                ],
            )
            aero_num = self.f_py.fortran_to_python(
                f_aero_num,
                [
                    self.flags.npx,
                    self.flags.npy,
                    self.flags.npz,
                    self.flags.n_modes,
                ],
            )
            aero_hygroscopicity = self.f_py.fortran_to_python(
                f_aero_hygroscopicity,
                [
                    self.flags.npx,
                    self.flags.npy,
                    self.flags.npz,
                    self.flags.n_modes,
                ],
            )
            aero_sigma = self.f_py.fortran_to_python(
                f_aero_sigma,
                [
                    self.flags.npx,
                    self.flags.npy,
                    self.flags.npz,
                    self.flags.n_modes,
                ],
            )

            frland = self.f_py.fortran_to_python(
                f_frland, [self.flags.npx, self.flags.npy]
            )

            t = self.f_py.fortran_to_python(f_t)
            plo = self.f_py.fortran_to_python(f_plo)
            qicn = self.f_py.fortran_to_python(f_qicn)
            qils = self.f_py.fortran_to_python(f_qils)
            qlcn = self.f_py.fortran_to_python(f_qlcn)
            qlls = self.f_py.fortran_to_python(f_qlls)
            vvel = self.f_py.fortran_to_python(f_vvel)
            tke = self.f_py.fortran_to_python(f_tke)
            nacti = self.f_py.fortran_to_python(f_nacti)
            nwfa = self.f_py.fortran_to_python(f_nwfa)
            nactl = self.f_py.fortran_to_python(f_nactl)
            self.f_py.device_sync()

        # Run Aer Activation
        with TimedCUDAProfiler("Aer Activation numerics", self._timings):
            self.pymoist.aer_activation(
                aero_dgn=aero_dgn,
                aero_num=aero_num,
                aero_hygroscopicity=aero_hygroscopicity,
                aero_sigma=aero_sigma,
                frland=frland,
                nn_ocean=Float(f_nn_ocean),
                nn_land=Float(f_nn_land),
                t=t,
                plo=plo,
                qicn=qicn,
                qils=qils,
                qlcn=qlcn,
                qlls=qlls,
                vvel=vvel,
                tke=tke,
                nwfa=nwfa,
                nacti=nacti,
                nactl=nactl,
            )

        # Convert NumPy arrays back to Fortran
        with TimedCUDAProfiler("Python -> Fortran", self._timings):
            self.f_py.python_to_fortran(aero_dgn, f_aero_dgn)
            self.f_py.python_to_fortran(aero_num, f_aero_num)
            self.f_py.python_to_fortran(aero_hygroscopicity, f_aero_hygroscopicity)
            self.f_py.python_to_fortran(aero_sigma, f_aero_sigma)
            self.f_py.python_to_fortran(frland, f_frland)
            self.f_py.python_to_fortran(t, f_t)
            self.f_py.python_to_fortran(plo, f_plo)
            self.f_py.python_to_fortran(qicn, f_qicn)
            self.f_py.python_to_fortran(qils, f_qils)
            self.f_py.python_to_fortran(qlcn, f_qlcn)
            self.f_py.python_to_fortran(qlls, f_qlls)
            self.f_py.python_to_fortran(vvel, f_vvel)
            self.f_py.python_to_fortran(tke, f_tke)
            self.f_py.python_to_fortran(nacti, f_nacti)
            self.f_py.python_to_fortran(nwfa, f_nwfa)
            self.f_py.python_to_fortran(nactl, f_nactl)

    def GFDL_1M_evap_subl_hystpdf(
        self,
        f_EIS: cffi.FFI.CData,
        f_dw_land: np.float32,
        f_dw_ocean: np.float32,
        f_PDF_shape: np.int32,
        f_TurnRHCrit_Param: cffi.FFI.CData,
        f_PLmb: cffi.FFI.CData,
        f_KLCL: cffi.FFI.CData,
        f_PLEmb: cffi.FFI.CData,
        f_Area: cffi.FFI.CData,
        f_dt_moist: np.float32,
        f_CNV_FRC: cffi.FFI.CData,
        f_SRF_TYPE: cffi.FFI.CData,
        f_T: cffi.FFI.CData,
        f_QLCN: cffi.FFI.CData,
        f_QICN: cffi.FFI.CData,
        f_QLLS: cffi.FFI.CData,
        f_QILS: cffi.FFI.CData,
        f_CCW_EVAP_EFF: np.float32,
        f_CCI_EVAP_EFF: np.float32,
        f_Q: cffi.FFI.CData,
        f_CLLS: cffi.FFI.CData,
        f_CLCN: cffi.FFI.CData,
        f_NACTL: cffi.FFI.CData,
        f_NACTI: cffi.FFI.CData,
        f_QST: cffi.FFI.CData,
        f_SUBLC: cffi.FFI.CData,
        f_EVAPC: cffi.FFI.CData,
        f_RHX: cffi.FFI.CData,
        f_LMELTFRZ: bool = True,
    ):
        CUDAProfiler.start_cuda_profiler()

        with TimedCUDAProfiler("[GFDL_1M] Fortran -> Python", self._timings):
            eis = self.f_py.fortran_to_python(
                f_EIS,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            klcl = self.f_py.fortran_to_python(
                f_KLCL,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            aera = self.f_py.fortran_to_python(
                f_Area,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            cnv_fraction = self.f_py.fortran_to_python(
                f_CNV_FRC,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            surface_type = self.f_py.fortran_to_python(
                f_SRF_TYPE,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )

            plmb = self.f_py.fortran_to_python(f_PLmb)
            plemb = self.f_py.fortran_to_python(f_PLEmb)
            T = self.f_py.fortran_to_python(f_T)
            QLCN = self.f_py.fortran_to_python(f_QLCN)
            QICN = self.f_py.fortran_to_python(f_QICN)
            QLLS = self.f_py.fortran_to_python(f_QLLS)
            QILS = self.f_py.fortran_to_python(f_QILS)
            Q = self.f_py.fortran_to_python(f_Q)
            CLLS = self.f_py.fortran_to_python(f_CLLS)
            CLCN = self.f_py.fortran_to_python(f_CLCN)
            NACTL = self.f_py.fortran_to_python(f_NACTL)
            NACTI = self.f_py.fortran_to_python(f_NACTI)
            QST = self.f_py.fortran_to_python(f_QST)

            self.f_py.device_sync()

        with TimedCUDAProfiler("[GFDL_1M] Run", self._timings):
            self.pymoist.GFDL_1M_evap(
                EIS=eis,
                dw_land=Float(f_dw_land),
                dw_ocean=Float(f_dw_ocean),
                PDFSHAPE=Float(f_PDF_shape),
                TURNRHCRIT_PARAM=Float(f_TurnRHCrit_Param),
                PLmb=plmb,
                KLCL=klcl,
                PLEmb=plemb,
                AREA=aera,
                DT_MOIST=Float(f_dt_moist),
                CNV_FRC=cnv_fraction,
                SRF_TYPE=surface_type,
                T=T,
                QLCN=QLCN,
                QICN=QICN,
                QLLS=QLLS,
                QILS=QILS,
                CCW_EVAP_EFF=Float(f_CCW_EVAP_EFF),
                CCI_EVAP_EFF=Float(f_CCI_EVAP_EFF),
                Q=Q,
                CLLS=CLLS,
                CLCN=CLCN,
                NACTL=NACTL,
                NACTI=NACTI,
                QST=QST,
                LMELTFRZ=f_LMELTFRZ,
            )

        # Convert NumPy arrays back to Fortran
        with TimedCUDAProfiler("[GFDL_1M] Python -> Fortran", self._timings):
            self.f_py.python_to_fortran(eis, f_EIS)
            self.f_py.python_to_fortran(klcl, f_KLCL)
            self.f_py.python_to_fortran(aera, f_Area)
            self.f_py.python_to_fortran(cnv_fraction, f_CNV_FRC)
            self.f_py.python_to_fortran(surface_type, f_SRF_TYPE)
            self.f_py.python_to_fortran(plmb, f_PLmb)
            self.f_py.python_to_fortran(plemb, f_PLEmb)
            self.f_py.python_to_fortran(T, f_T)
            self.f_py.python_to_fortran(QLCN, f_QLCN)
            self.f_py.python_to_fortran(QICN, f_QICN)
            self.f_py.python_to_fortran(QLLS, f_QLLS)
            self.f_py.python_to_fortran(QILS, f_QILS)
            self.f_py.python_to_fortran(Q, f_Q)
            self.f_py.python_to_fortran(CLLS, f_CLLS)
            self.f_py.python_to_fortran(CLCN, f_CLCN)
            self.f_py.python_to_fortran(NACTL, f_NACTL)
            self.f_py.python_to_fortran(NACTI, f_NACTI)
            self.f_py.python_to_fortran(QST, f_QST)
            self.f_py.python_to_fortran(self.pymoist.gfdl_1M.sublc.view[:], f_SUBLC)
            self.f_py.python_to_fortran(self.pymoist.gfdl_1M.evapc.view[:], f_EVAPC)
            self.f_py.python_to_fortran(self.pymoist.gfdl_1M.rhx.view[:], f_RHX)

    def driver(
        self,
        f_qv: cffi.FFI.CData,
        f_ql: cffi.FFI.CData,
        f_qr: cffi.FFI.CData,
        f_qi: cffi.FFI.CData,
        f_qs: cffi.FFI.CData,
        f_qg: cffi.FFI.CData,
        f_qa: cffi.FFI.CData,
        f_qn: cffi.FFI.CData,
        f_qv_dt: cffi.FFI.CData,
        f_ql_dt: cffi.FFI.CData,
        f_qr_dt: cffi.FFI.CData,
        f_qi_dt: cffi.FFI.CData,
        f_qs_dt: cffi.FFI.CData,
        f_qg_dt: cffi.FFI.CData,
        f_qa_dt: cffi.FFI.CData,
        f_t_dt: cffi.FFI.CData,
        f_t: cffi.FFI.CData,
        f_w: cffi.FFI.CData,
        f_uin: cffi.FFI.CData,
        f_vin: cffi.FFI.CData,
        f_udt: cffi.FFI.CData,
        f_vdt: cffi.FFI.CData,
        f_dz: cffi.FFI.CData,
        f_dp: cffi.FFI.CData,
        f_area: cffi.FFI.CData,
        f_fr_land: cffi.FFI.CData,
        f_cnv_frc: cffi.FFI.CData,
        f_srf_type: cffi.FFI.CData,
        f_eis: cffi.FFI.CData,
        f_rhcrit3d: cffi.FFI.CData,
        f_anv_icefall: np.float32,
        f_ls_icefall: np.float32,
        f_hydrostatic: bool,
        f_phys_hydrostatic: bool,
        f_rain: cffi.FFI.CData,
        f_snow: cffi.FFI.CData,
        f_ice: cffi.FFI.CData,
        f_graupel: cffi.FFI.CData,
        f_m2_rain: cffi.FFI.CData,
        f_m2_sol: cffi.FFI.CData,
        f_revap: cffi.FFI.CData,
        f_isubl: cffi.FFI.CData,
    ):
        CUDAProfiler.start_cuda_profiler()

        with TimedCUDAProfiler("[GFDL_1M] Fortran -> Python", self._timings):
            qv = self.f_py.fortran_to_python(f_qv)
            ql = self.f_py.fortran_to_python(f_ql)
            qr = self.f_py.fortran_to_python(f_qr)
            qi = self.f_py.fortran_to_python(f_qi)
            qs = self.f_py.fortran_to_python(f_qs)
            qg = self.f_py.fortran_to_python(f_qg)
            qa = self.f_py.fortran_to_python(f_qa)
            qn = self.f_py.fortran_to_python(f_qn)
            qv_dt = self.f_py.fortran_to_python(f_qv_dt)
            ql_dt = self.f_py.fortran_to_python(f_ql_dt)
            qr_dt = self.f_py.fortran_to_python(f_qr_dt)
            qi_dt = self.f_py.fortran_to_python(f_qi_dt)
            qs_dt = self.f_py.fortran_to_python(f_qs_dt)
            qg_dt = self.f_py.fortran_to_python(f_qg_dt)
            qa_dt = self.f_py.fortran_to_python(f_qa_dt)
            t_dt = self.f_py.fortran_to_python(f_t_dt)
            t = self.f_py.fortran_to_python(f_t)
            w = self.f_py.fortran_to_python(f_w)
            uin = self.f_py.fortran_to_python(f_uin)
            vin = self.f_py.fortran_to_python(f_vin)
            udt = self.f_py.fortran_to_python(f_udt)
            vdt = self.f_py.fortran_to_python(f_vdt)
            dz = self.f_py.fortran_to_python(f_dz)
            dp = self.f_py.fortran_to_python(f_dp)
            area = self.f_py.fortran_to_python(
                f_area,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            fr_land = self.f_py.fortran_to_python(
                f_fr_land,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            cnv_frc = self.f_py.fortran_to_python(
                f_cnv_frc,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            srf_type = self.f_py.fortran_to_python(
                f_srf_type,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            eis = self.f_py.fortran_to_python(
                f_eis,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            rhcrit3d = self.f_py.fortran_to_python(f_rhcrit3d)

            self.f_py.device_sync()

        with TimedCUDAProfiler("[GFDL_1M] Run", self._timings):
            self.pymoist.driver(
                qv=qv,
                ql=ql,
                qr=qr,
                qi=qi,
                qs=qs,
                qg=qg,
                qa=qa,
                qn=qn,  # NACTL + NACTI
                qv_dt=qv_dt,
                ql_dt=ql_dt,
                qr_dt=qr_dt,
                qi_dt=qi_dt,
                qs_dt=qs_dt,
                qg_dt=qg_dt,
                qa_dt=qa_dt,
                t_dt=t_dt,
                t=t,
                w=w,
                u=uin,
                v=vin,
                u_dt=udt,
                v_dt=vdt,
                dz=dz,
                dp=dp,
                area=area,
                fr_land=fr_land,
                cnv_frc=cnv_frc,
                srf_type=srf_type,
                eis=eis,
                rhcrit3d=rhcrit3d,
                anv_icefall=Float(f_anv_icefall),
                ls_icefall=Float(f_ls_icefall),
            )

        # Convert NumPy arrays back to Fortran
        with TimedCUDAProfiler("[GFDL_1M] Python -> Fortran", self._timings):
            self.f_py.python_to_fortran(qv, f_qv)
            self.f_py.python_to_fortran(ql, f_ql)
            self.f_py.python_to_fortran(qr, f_qr)
            self.f_py.python_to_fortran(qi, f_qi)
            self.f_py.python_to_fortran(qs, f_qs)
            self.f_py.python_to_fortran(qg, f_qg)
            self.f_py.python_to_fortran(qa, f_qa)
            self.f_py.python_to_fortran(qn, f_qn)
            self.f_py.python_to_fortran(qv_dt, f_qv_dt)
            self.f_py.python_to_fortran(ql_dt, f_ql_dt)
            self.f_py.python_to_fortran(qr_dt, f_qr_dt)
            self.f_py.python_to_fortran(qi_dt, f_qi_dt)
            self.f_py.python_to_fortran(qs_dt, f_qs_dt)
            self.f_py.python_to_fortran(qg_dt, f_qg_dt)
            self.f_py.python_to_fortran(qa_dt, f_qa_dt)
            self.f_py.python_to_fortran(t_dt, f_t_dt)
            self.f_py.python_to_fortran(t, f_t)
            self.f_py.python_to_fortran(w, f_w)
            self.f_py.python_to_fortran(uin, f_uin)
            self.f_py.python_to_fortran(vin, f_vin)
            self.f_py.python_to_fortran(udt, f_udt)
            self.f_py.python_to_fortran(vdt, f_vdt)
            self.f_py.python_to_fortran(dz, f_dz)
            self.f_py.python_to_fortran(dp, f_dp)
            self.f_py.python_to_fortran(area, f_area)
            self.f_py.python_to_fortran(fr_land, f_fr_land)
            self.f_py.python_to_fortran(cnv_frc, f_cnv_frc)
            self.f_py.python_to_fortran(srf_type, f_srf_type)
            self.f_py.python_to_fortran(eis, f_eis)
            self.f_py.python_to_fortran(rhcrit3d, f_rhcrit3d)
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.rain.view[:], f_rain
            )
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.snow.view[:], f_snow
            )
            self.f_py.python_to_fortran(self.pymoist.driver.outputs.ice.view[:], f_ice)
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.graupel.view[:], f_graupel
            )
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.m2_rain.view[:], f_m2_rain
            )
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.m2_sol.view[:], f_m2_sol
            )
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.revap.view[:], f_revap
            )
            self.f_py.python_to_fortran(
                self.pymoist.driver.outputs.isubl.view[:], f_isubl
            )

    def UW_shallow_convection(
        self,
        dotransport: int,
        k0: int,
        windsrcavg: int,
        qtsrchgt: float,
        qtsrc_fac: float,
        thlsrc_fac: float,
        frc_rasn: float,
        rbuoy: float,
        epsvarw: float,
        use_CINcin: int,
        mumin1: float,
        rmaxfrac: float,
        PGFc: float,
        dt: float,
        niter_xc: int,
        criqc: float,
        rle: float,
        cridist_opt: int,
        mixscale: float,
        rdrag: float,
        rkm: float,
        use_self_detrain: int,
        detrhgt: float,
        use_cumpenent: int,
        rpen: float,
        use_momenflx: int,
        rdrop: float,
        iter_cin: int,
        f_pifc0_inv: cffi.FFI.CData,
        f_zifc0_inv: cffi.FFI.CData,
        f_pmid0_inv: cffi.FFI.CData,
        f_zmid0_inv: cffi.FFI.CData,
        f_kpbl_inv: cffi.FFI.CData,
        f_exnmid0_inv: cffi.FFI.CData,
        f_exnifc0_inv: cffi.FFI.CData,
        f_dp0_inv: cffi.FFI.CData,
        f_u0_inv: cffi.FFI.CData,
        f_v0_inv: cffi.FFI.CData,
        f_qv0_inv: cffi.FFI.CData,
        f_ql0_inv: cffi.FFI.CData,
        f_qi0_inv: cffi.FFI.CData,
        f_t0_inv: cffi.FFI.CData,
        f_frland: cffi.FFI.CData,
        f_tke_inv: cffi.FFI.CData,
        f_rkfre: cffi.FFI.CData,
        f_cush: cffi.FFI.CData,
        f_shfx: cffi.FFI.CData,
        f_evap: cffi.FFI.CData,
        f_cnvtr: cffi.FFI.CData,
        f_CNV_Tracers: cffi.FFI.CData,
        f_umf_inv: cffi.FFI.CData,
        f_dcm_inv: cffi.FFI.CData,
        f_qtflx_inv: cffi.FFI.CData,
        f_slflx_inv: cffi.FFI.CData,
        f_uflx_inv: cffi.FFI.CData,
        f_vflx_inv: cffi.FFI.CData,
        f_qvten_inv: cffi.FFI.CData,
        f_qlten_inv: cffi.FFI.CData,
        f_qiten_inv: cffi.FFI.CData,
        f_tten_inv: cffi.FFI.CData,
        f_uten_inv: cffi.FFI.CData,
        f_vten_inv: cffi.FFI.CData,
        f_qrten_inv: cffi.FFI.CData,
        f_qsten_inv: cffi.FFI.CData,
        f_cufrc_inv: cffi.FFI.CData,
        f_fer_inv: cffi.FFI.CData,
        f_fdr_inv: cffi.FFI.CData,
        f_ndrop_inv: cffi.FFI.CData,
        f_nice_inv: cffi.FFI.CData,
        f_qldet_inv: cffi.FFI.CData,
        f_qlsub_inv: cffi.FFI.CData,
        f_qidet_inv: cffi.FFI.CData,
        f_qisub_inv: cffi.FFI.CData,
        f_tpert_out: cffi.FFI.CData,
        f_qpert_out: cffi.FFI.CData,
    ):
        CUDAProfiler.start_cuda_profiler()

        with TimedCUDAProfiler("[UW] Fortran -> Python", self._timings):
            in_pifc0_inv = self.f_py.fortran_to_python(
                f_pifc0_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            in_zifc0_inv = self.f_py.fortran_to_python(
                f_zifc0_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            in_pmid0_inv = self.f_py.fortran_to_python(f_pmid0_inv)
            in_zmid0_inv = self.f_py.fortran_to_python(f_zmid0_inv)
            in_kpbl_inv = self.f_py.fortran_to_python(
                f_kpbl_inv,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_exnmid0_inv = self.f_py.fortran_to_python(f_exnmid0_inv)
            in_exnifc0_inv = self.f_py.fortran_to_python(
                f_exnifc0_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            in_dp0_inv = self.f_py.fortran_to_python(f_dp0_inv)
            in_u0_inv = self.f_py.fortran_to_python(f_u0_inv)
            in_v0_inv = self.f_py.fortran_to_python(f_v0_inv)
            in_qv0_inv = self.f_py.fortran_to_python(f_qv0_inv)
            in_ql0_inv = self.f_py.fortran_to_python(f_ql0_inv)
            in_qi0_inv = self.f_py.fortran_to_python(f_qi0_inv)
            in_t0_inv = self.f_py.fortran_to_python(f_t0_inv)
            in_frland = self.f_py.fortran_to_python(
                f_frland,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_tke_inv = self.f_py.fortran_to_python(
                f_tke_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            in_rkfre = self.f_py.fortran_to_python(
                f_rkfre,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_cush = self.f_py.fortran_to_python(
                f_cush,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_shfx = self.f_py.fortran_to_python(
                f_shfx,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_evap = self.f_py.fortran_to_python(
                f_evap,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_cnvtr = self.f_py.fortran_to_python(
                f_cnvtr,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            in_CNV_Tracers = self.f_py.fortran_to_python(
                f_CNV_Tracers,
                [
                    self.flags.npx,
                    self.flags.npy,
                    self.flags.npz,
                    self.pymoist.UW_config.NCNST,
                ],
            )
            self.f_py.fortran_to_python(f_CNV_Tracers)

            out_umf_inv = self.f_py.fortran_to_python(
                f_umf_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            out_dcm_inv = self.f_py.fortran_to_python(f_dcm_inv)
            out_qtflx_inv = self.f_py.fortran_to_python(
                f_qtflx_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            out_slflx_inv = self.f_py.fortran_to_python(
                f_slflx_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            out_uflx_inv = self.f_py.fortran_to_python(
                f_uflx_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            out_vflx_inv = self.f_py.fortran_to_python(
                f_vflx_inv,
                [self.flags.npx, self.flags.npy, self.flags.npz + 1],
            )
            out_qvten_inv = self.f_py.fortran_to_python(f_qvten_inv)
            out_qlten_inv = self.f_py.fortran_to_python(f_qlten_inv)
            out_qiten_inv = self.f_py.fortran_to_python(f_qiten_inv)
            out_tten_inv = self.f_py.fortran_to_python(f_tten_inv)
            out_uten_inv = self.f_py.fortran_to_python(f_uten_inv)
            out_vten_inv = self.f_py.fortran_to_python(f_vten_inv)
            out_qrten_inv = self.f_py.fortran_to_python(f_qrten_inv)
            out_qsten_inv = self.f_py.fortran_to_python(f_qsten_inv)
            out_cufrc_inv = self.f_py.fortran_to_python(f_cufrc_inv)
            out_fer_inv = self.f_py.fortran_to_python(f_fer_inv)
            out_fdr_inv = self.f_py.fortran_to_python(f_fdr_inv)
            out_ndrop_inv = self.f_py.fortran_to_python(f_ndrop_inv)
            out_nice_inv = self.f_py.fortran_to_python(f_nice_inv)
            out_qldet_inv = self.f_py.fortran_to_python(f_qldet_inv)
            out_qlsub_inv = self.f_py.fortran_to_python(f_qlsub_inv)
            out_qidet_inv = self.f_py.fortran_to_python(f_qidet_inv)
            out_qisub_inv = self.f_py.fortran_to_python(f_qisub_inv)
            out_tpert_out = self.f_py.fortran_to_python(
                f_tpert_out,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )
            out_qpert_out = self.f_py.fortran_to_python(
                f_qpert_out,
                [
                    self.flags.npx,
                    self.flags.npy,
                ],
            )

        with TimedCUDAProfiler("[UW] Run", self._timings):
            self.pymoist.UW_shallow_convection(
                # Inputs
                pifc0_inv=in_pifc0_inv,
                zifc0_inv=in_zifc0_inv,
                pmid0_inv=in_pmid0_inv,
                zmid0_inv=in_zmid0_inv,
                kpbl_inv=in_kpbl_inv,
                exnmid0_inv=in_exnmid0_inv,
                exnifc0_inv=in_exnifc0_inv,
                dp0_inv=in_dp0_inv,
                u0_inv=in_u0_inv,
                v0_inv=in_v0_inv,
                qv0_inv=in_qv0_inv,
                ql0_inv=in_ql0_inv,
                qi0_inv=in_qi0_inv,
                t0_inv=in_t0_inv,
                frland=in_frland,
                tke_inv=in_tke_inv,
                rkfre=in_rkfre,
                cush=in_cush,
                shfx=in_shfx,
                evap=in_evap,
                cnvtr=in_cnvtr,
                CNV_Tracers=in_CNV_Tracers,
                # Parameters
                dotransport=dotransport,
                k0=k0,
                windsrcavg=windsrcavg,
                qtsrchgt=qtsrchgt,
                qtsrc_fac=qtsrc_fac,
                thlsrc_fac=thlsrc_fac,
                frc_rasn=frc_rasn,
                rbuoy=rbuoy,
                epsvarw=epsvarw,
                use_CINcin=use_CINcin,
                mumin1=mumin1,
                rmaxfrac=rmaxfrac,
                PGFc=PGFc,
                dt=dt,
                niter_xc=niter_xc,
                criqc=criqc,
                rle=rle,
                cridist_opt=cridist_opt,
                mixscale=mixscale,
                rdrag=rdrag,
                rkm=rkm,
                use_self_detrain=use_self_detrain,
                detrhgt=detrhgt,
                use_cumpenent=use_cumpenent,
                rpen=rpen,
                use_momenflx=use_momenflx,
                rdrop=rdrop,
                iter_cin=iter_cin,
                # Outputs
                umf_inv=out_umf_inv,
                dcm_inv=out_dcm_inv,
                qtflx_inv=out_qtflx_inv,
                slflx_inv=out_slflx_inv,
                uflx_inv=out_uflx_inv,
                vflx_inv=out_vflx_inv,
                qvten_inv=out_qvten_inv,
                qlten_inv=out_qlten_inv,
                qiten_inv=out_qiten_inv,
                tten_inv=out_tten_inv,
                uten_inv=out_uten_inv,
                vten_inv=out_vten_inv,
                qrten_inv=out_qrten_inv,
                qsten_inv=out_qsten_inv,
                cufrc_inv=out_cufrc_inv,
                fer_inv=out_fer_inv,
                fdr_inv=out_fdr_inv,
                ndrop_inv=out_ndrop_inv,
                nice_inv=out_nice_inv,
                qldet_inv=out_qldet_inv,
                qlsub_inv=out_qlsub_inv,
                qidet_inv=out_qidet_inv,
                qisub_inv=out_qisub_inv,
                tpert_out=out_tpert_out,
                qpert_out=out_qpert_out,
            )

        with TimedCUDAProfiler("[GFDL_1M] Python -> Fortran", self._timings):
            self.f_py.python_to_fortran(out_umf_inv, f_umf_inv)
            self.f_py.python_to_fortran(out_dcm_inv, f_dcm_inv)
            self.f_py.python_to_fortran(out_qtflx_inv, f_qtflx_inv)
            self.f_py.python_to_fortran(out_slflx_inv, f_slflx_inv)
            self.f_py.python_to_fortran(out_uflx_inv, f_uflx_inv)
            self.f_py.python_to_fortran(out_vflx_inv, f_vflx_inv)
            self.f_py.python_to_fortran(out_qvten_inv, f_qvten_inv)
            self.f_py.python_to_fortran(out_qlten_inv, f_qlten_inv)
            self.f_py.python_to_fortran(out_qiten_inv, f_qiten_inv)
            self.f_py.python_to_fortran(out_tten_inv, f_tten_inv)
            self.f_py.python_to_fortran(out_uten_inv, f_uten_inv)
            self.f_py.python_to_fortran(out_vten_inv, f_vten_inv)
            self.f_py.python_to_fortran(out_qrten_inv, f_qrten_inv)
            self.f_py.python_to_fortran(out_qsten_inv, f_qsten_inv)
            self.f_py.python_to_fortran(out_cufrc_inv, f_cufrc_inv)
            self.f_py.python_to_fortran(out_fer_inv, f_fer_inv)
            self.f_py.python_to_fortran(out_fdr_inv, f_fdr_inv)
            self.f_py.python_to_fortran(out_ndrop_inv, f_ndrop_inv)
            self.f_py.python_to_fortran(out_nice_inv, f_nice_inv)
            self.f_py.python_to_fortran(out_qldet_inv, f_qldet_inv)
            self.f_py.python_to_fortran(out_qlsub_inv, f_qlsub_inv)
            self.f_py.python_to_fortran(out_qidet_inv, f_qidet_inv)
            self.f_py.python_to_fortran(out_qisub_inv, f_qisub_inv)
            self.f_py.python_to_fortran(out_tpert_out, f_tpert_out)
            self.f_py.python_to_fortran(out_qpert_out, f_qpert_out)


WRAPPER = PYMOIST_WRAPPER()


def pyMoist_run_GFDL_1M_evap_subl_hystpdf(
    dw_land: np.float32,
    dw_ocean: np.float32,
    PDFSHAPE: np.int32,
    TURNRHCRIT_PARAM: np.float32,
    DT_MOIST: np.float32,
    CCW_EVAP_EFF: np.float32,
    CCI_EVAP_EFF: np.float32,
    LMELTFRZ: bool,
    AREA: cffi.FFI.CData,
    CNV_FRC: cffi.FFI.CData,
    SRF_TYPE: cffi.FFI.CData,
    KLCL: cffi.FFI.CData,
    EIS: cffi.FFI.CData,
    PLmb: cffi.FFI.CData,
    PLEmb: cffi.FFI.CData,
    NACTL: cffi.FFI.CData,
    NACTI: cffi.FFI.CData,
    QST: cffi.FFI.CData,
    T: cffi.FFI.CData,
    Q: cffi.FFI.CData,
    QLCN: cffi.FFI.CData,
    QICN: cffi.FFI.CData,
    QLLS: cffi.FFI.CData,
    QILS: cffi.FFI.CData,
    CLLS: cffi.FFI.CData,
    CLCN: cffi.FFI.CData,
    SUBLC: cffi.FFI.CData,
    EVAPC: cffi.FFI.CData,
    RHX: cffi.FFI.CData,
):
    if not WRAPPER.ready:
        raise RuntimeError("[pyMoist] Bad init, did you call init?")
    WRAPPER.GFDL_1M_evap_subl_hystpdf(
        f_EIS=EIS,
        f_dw_land=dw_land,
        f_dw_ocean=dw_ocean,
        f_PDF_shape=PDFSHAPE,
        f_TurnRHCrit_Param=TURNRHCRIT_PARAM,
        f_PLmb=PLmb,
        f_PLEmb=PLEmb,
        f_Area=AREA,
        f_dt_moist=DT_MOIST,
        f_CNV_FRC=CNV_FRC,
        f_SRF_TYPE=SRF_TYPE,
        f_KLCL=KLCL,
        f_T=T,
        f_QLCN=QLCN,
        f_QICN=QICN,
        f_QLLS=QLLS,
        f_QILS=QILS,
        f_CCI_EVAP_EFF=CCI_EVAP_EFF,
        f_CCW_EVAP_EFF=CCW_EVAP_EFF,
        f_Q=Q,
        f_CLLS=CLLS,
        f_CLCN=CLCN,
        f_NACTL=NACTL,
        f_NACTI=NACTI,
        f_QST=QST,
        f_LMELTFRZ=LMELTFRZ,
        f_SUBLC=SUBLC,
        f_EVAPC=EVAPC,
        f_RHX=RHX,
    )


def pymoist_run_GFDL_1M_driver(
    RAD_QV: cffi.FFI.CData,
    RAD_QL: cffi.FFI.CData,
    RAD_QR: cffi.FFI.CData,
    RAD_QI: cffi.FFI.CData,
    RAD_QS: cffi.FFI.CData,
    RAD_QG: cffi.FFI.CData,
    RAD_CF: cffi.FFI.CData,
    NACTAll: cffi.FFI.CData,
    DQVDTmic: cffi.FFI.CData,
    DQLDTmic: cffi.FFI.CData,
    DQRDTmic: cffi.FFI.CData,
    DQIDTmic: cffi.FFI.CData,
    DQSDTmic: cffi.FFI.CData,
    DQGDTmic: cffi.FFI.CData,
    DQADTmic: cffi.FFI.CData,
    DTDTmic: cffi.FFI.CData,
    T: cffi.FFI.CData,
    W: cffi.FFI.CData,
    U: cffi.FFI.CData,
    V: cffi.FFI.CData,
    DUDTmic: cffi.FFI.CData,
    DVDTmic: cffi.FFI.CData,
    DZ: cffi.FFI.CData,
    DP: cffi.FFI.CData,
    AREA: cffi.FFI.CData,
    FRLAND: cffi.FFI.CData,
    CNV_FRC: cffi.FFI.CData,
    SRF_TYPE: cffi.FFI.CData,
    EIS: cffi.FFI.CData,
    RHCRIT3D: cffi.FFI.CData,
    DT_MOIST: np.float32,
    ANV_ICEFALL: np.float32,
    LS_ICEFALL: np.float32,
    REV_LS: cffi.FFI.CData,
    RSU_LS: cffi.FFI.CData,
    PRCP_RAIN: cffi.FFI.CData,
    PRCP_SNOW: cffi.FFI.CData,
    PRCP_ICE: cffi.FFI.CData,
    PRCP_GRAUPEL: cffi.FFI.CData,
    PFL_LS: cffi.FFI.CData,
    PFI_LS: cffi.FFI.CData,
    LHYDROSTATIC: bool,
    LPHYS_HYDROSTATIC: bool,
):
    if not WRAPPER.ready:
        raise RuntimeError("[pyMoist] Bad init, did you call init?")

    WRAPPER.driver(
        f_qv=RAD_QV,
        f_ql=RAD_QL,
        f_qr=RAD_QR,
        f_qi=RAD_QI,
        f_qs=RAD_QS,
        f_qg=RAD_QG,
        f_qa=RAD_CF,
        f_qn=NACTAll,
        f_qv_dt=DQVDTmic,
        f_ql_dt=DQLDTmic,
        f_qr_dt=DQRDTmic,
        f_qi_dt=DQIDTmic,
        f_qs_dt=DQSDTmic,
        f_qg_dt=DQGDTmic,
        f_qa_dt=DQADTmic,
        f_t_dt=DTDTmic,
        f_t=T,
        f_w=W,
        f_uin=U,
        f_vin=V,
        f_udt=DUDTmic,
        f_vdt=DVDTmic,
        f_dz=DZ,
        f_dp=DP,
        f_area=AREA,
        f_fr_land=FRLAND,
        f_cnv_frc=CNV_FRC,
        f_srf_type=SRF_TYPE,
        f_eis=EIS,
        f_rhcrit3d=RHCRIT3D,
        f_anv_icefall=ANV_ICEFALL,
        f_ls_icefall=LS_ICEFALL,
        f_hydrostatic=LHYDROSTATIC,
        f_phys_hydrostatic=LPHYS_HYDROSTATIC,
        f_rain=PRCP_RAIN,
        f_snow=PRCP_SNOW,
        f_ice=PRCP_ICE,
        f_graupel=PRCP_GRAUPEL,
        f_m2_rain=PFL_LS,
        f_m2_sol=PFI_LS,
        f_revap=REV_LS,
        f_isubl=RSU_LS,
    )


def pyMoist_run_AerActivation(
    aero_dgn: cffi.FFI.CData,
    aero_num: cffi.FFI.CData,
    aero_hygroscopicity: cffi.FFI.CData,
    aero_sigma: cffi.FFI.CData,
    frland: cffi.FFI.CData,
    nn_ocean: np.float32,
    nn_land: np.float32,
    t: cffi.FFI.CData,
    plo: cffi.FFI.CData,
    qicn: cffi.FFI.CData,
    qils: cffi.FFI.CData,
    qlcn: cffi.FFI.CData,
    qlls: cffi.FFI.CData,
    vvel: cffi.FFI.CData,
    tke: cffi.FFI.CData,
    nacti: cffi.FFI.CData,
    nwfa: cffi.FFI.CData,
    nactl: cffi.FFI.CData,
):
    if not WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.aer_activation(
        aero_dgn,
        aero_num,
        aero_hygroscopicity,
        aero_sigma,
        frland,
        nn_ocean,
        nn_land,
        t,
        plo,
        qicn,
        qils,
        qlcn,
        qlls,
        vvel,
        tke,
        nacti,
        nwfa,
        nactl,
    )


def pyMoist_finalize():
    if WRAPPER.ready:
        WRAPPER.finalize()


def pyMoist_init(pyMoist_flags: cffi.FFI.CData):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "dace:cpu")
    if WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.init(
        pyMoist_flags=pyMoist_flags,
        backend=BACKEND,
    )


def gfdl_1m_init(gfdl_1m_flags: cffi.FFI.CData) -> None:
    if not WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.pymoist.init_gfdl_1m_configuration(
        flags=gfdl_1m_flags_f_to_python(gfdl_1m_flags),
    )


def compute_uwshcu_init(
    ncnst: int,
    k0: int,
    windsrcavg: int,
):
    print("Reached compute_uwschu_init")
    if not WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.pymoist.init_UW_configuration(
        NCNST=ncnst,
        k0=k0,
        windsrcavg=windsrcavg,
    )


def compute_uwshcu_run(
    # inputs
    dotransport: int,
    k0: int,
    windsrcavg: int,
    qtsrchgt: float,
    qtsrc_fac: float,
    thlsrc_fac: float,
    frc_rasn: float,
    rbuoy: float,
    epsvarw: float,
    use_CINcin: int,
    mumin1: float,
    rmaxfrac: float,
    PGFc: float,
    dt: float,
    niter_xc: int,
    criqc: float,
    rle: float,
    cridist_opt: int,
    mixscale: float,
    rdrag: float,
    rkm: float,
    use_self_detrain: int,
    detrhgt: float,
    use_cumpenent: int,
    rpen: float,
    use_momenflx: int,
    rdrop: float,
    iter_cin: int,
    pifc0_inv: "cffi.FFI.CData",
    pifc0_inv_dim_sizes: "cffi.FFI.CData",
    pifc0_inv_rank: int,
    zifc0_inv: "cffi.FFI.CData",
    zifc0_inv_dim_sizes: "cffi.FFI.CData",
    zifc0_inv_rank: int,
    pmid0_inv: "cffi.FFI.CData",
    pmid0_inv_dim_sizes: "cffi.FFI.CData",
    pmid0_inv_rank: int,
    zmid0_inv: "cffi.FFI.CData",
    zmid0_inv_dim_sizes: "cffi.FFI.CData",
    zmid0_inv_rank: int,
    kpbl_inv: "cffi.FFI.CData",
    kpbl_inv_dim_sizes: "cffi.FFI.CData",
    kpbl_inv_rank: int,
    exnmid0_inv: "cffi.FFI.CData",
    exnmid0_inv_dim_sizes: "cffi.FFI.CData",
    exnmid0_inv_rank: int,
    exnifc0_inv: "cffi.FFI.CData",
    exnifc0_inv_dim_sizes: "cffi.FFI.CData",
    exnifc0_inv_rank: int,
    dp0_inv: "cffi.FFI.CData",
    dp0_inv_dim_sizes: "cffi.FFI.CData",
    dp0_inv_rank: int,
    u0_inv: "cffi.FFI.CData",
    u0_inv_dim_sizes: "cffi.FFI.CData",
    u0_inv_rank: int,
    v0_inv: "cffi.FFI.CData",
    v0_inv_dim_sizes: "cffi.FFI.CData",
    v0_inv_rank: int,
    qv0_inv: "cffi.FFI.CData",
    qv0_inv_dim_sizes: "cffi.FFI.CData",
    qv0_inv_rank: int,
    ql0_inv: "cffi.FFI.CData",
    ql0_inv_dim_sizes: "cffi.FFI.CData",
    ql0_inv_rank: int,
    qi0_inv: "cffi.FFI.CData",
    qi0_inv_dim_sizes: "cffi.FFI.CData",
    qi0_inv_rank: int,
    t0_inv: "cffi.FFI.CData",
    t0_inv_dim_sizes: "cffi.FFI.CData",
    t0_inv_rank: int,
    frland_in: "cffi.FFI.CData",
    frland_in_dim_sizes: "cffi.FFI.CData",
    frland_in_rank: int,
    tke_inv: "cffi.FFI.CData",
    tke_inv_dim_sizes: "cffi.FFI.CData",
    tke_inv_rank: int,
    rkfre: "cffi.FFI.CData",
    rkfre_dim_sizes: "cffi.FFI.CData",
    rkfre_rank: int,
    cush: "cffi.FFI.CData",
    cush_dim_sizes: "cffi.FFI.CData",
    cush_rank: int,
    shfx: "cffi.FFI.CData",
    shfx_dim_sizes: "cffi.FFI.CData",
    shfx_rank: int,
    evap: "cffi.FFI.CData",
    evap_dim_sizes: "cffi.FFI.CData",
    evap_rank: int,
    cnvtr: "cffi.FFI.CData",
    cnvtr_dim_sizes: "cffi.FFI.CData",
    cnvtr_rank: int,
    CNV_Tracers: "cffi.FFI.CData",
    CNV_Tracers_dim_sizes: "cffi.FFI.CData",
    CNV_Tracers_rank: int,
    # inputs-outputs
    # outputs
    umf_inv: "cffi.FFI.CData",
    umf_inv_dim_sizes: "cffi.FFI.CData",
    umf_inv_rank: int,
    dcm_inv: "cffi.FFI.CData",
    dcm_inv_dim_sizes: "cffi.FFI.CData",
    dcm_inv_rank: int,
    qtflx_inv: "cffi.FFI.CData",
    qtflx_inv_dim_sizes: "cffi.FFI.CData",
    qtflx_inv_rank: int,
    slflx_inv: "cffi.FFI.CData",
    slflx_inv_dim_sizes: "cffi.FFI.CData",
    slflx_inv_rank: int,
    uflx_inv: "cffi.FFI.CData",
    uflx_inv_dim_sizes: "cffi.FFI.CData",
    uflx_inv_rank: int,
    vflx_inv: "cffi.FFI.CData",
    vflx_inv_dim_sizes: "cffi.FFI.CData",
    vflx_inv_rank: int,
    qvten_inv: "cffi.FFI.CData",
    qvten_inv_dim_sizes: "cffi.FFI.CData",
    qvten_inv_rank: int,
    qlten_inv: "cffi.FFI.CData",
    qlten_inv_dim_sizes: "cffi.FFI.CData",
    qlten_inv_rank: int,
    qiten_inv: "cffi.FFI.CData",
    qiten_inv_dim_sizes: "cffi.FFI.CData",
    qiten_inv_rank: int,
    tten_inv: "cffi.FFI.CData",
    tten_inv_dim_sizes: "cffi.FFI.CData",
    tten_inv_rank: int,
    uten_inv: "cffi.FFI.CData",
    uten_inv_dim_sizes: "cffi.FFI.CData",
    uten_inv_rank: int,
    vten_inv: "cffi.FFI.CData",
    vten_inv_dim_sizes: "cffi.FFI.CData",
    vten_inv_rank: int,
    qrten_inv: "cffi.FFI.CData",
    qrten_inv_dim_sizes: "cffi.FFI.CData",
    qrten_inv_rank: int,
    qsten_inv: "cffi.FFI.CData",
    qsten_inv_dim_sizes: "cffi.FFI.CData",
    qsten_inv_rank: int,
    cufrc_inv: "cffi.FFI.CData",
    cufrc_inv_dim_sizes: "cffi.FFI.CData",
    cufrc_inv_rank: int,
    fer_inv: "cffi.FFI.CData",
    fer_inv_dim_sizes: "cffi.FFI.CData",
    fer_inv_rank: int,
    fdr_inv: "cffi.FFI.CData",
    fdr_inv_dim_sizes: "cffi.FFI.CData",
    fdr_inv_rank: int,
    ndrop_inv: "cffi.FFI.CData",
    ndrop_inv_dim_sizes: "cffi.FFI.CData",
    ndrop_inv_rank: int,
    nice_inv: "cffi.FFI.CData",
    nice_inv_dim_sizes: "cffi.FFI.CData",
    nice_inv_rank: int,
    qldet_inv: "cffi.FFI.CData",
    qldet_inv_dim_sizes: "cffi.FFI.CData",
    qldet_inv_rank: int,
    qlsub_inv: "cffi.FFI.CData",
    qlsub_inv_dim_sizes: "cffi.FFI.CData",
    qlsub_inv_rank: int,
    qidet_inv: "cffi.FFI.CData",
    qidet_inv_dim_sizes: "cffi.FFI.CData",
    qidet_inv_rank: int,
    qisub_inv: "cffi.FFI.CData",
    qisub_inv_dim_sizes: "cffi.FFI.CData",
    qisub_inv_rank: int,
    tpert_out: "cffi.FFI.CData",
    tpert_out_dim_sizes: "cffi.FFI.CData",
    tpert_out_rank: int,
    qpert_out: "cffi.FFI.CData",
    qpert_out_dim_sizes: "cffi.FFI.CData",
    qpert_out_rank: int,
):
    if not WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.UW_shallow_convection(
        dotransport=dotransport,
        k0=k0,
        windsrcavg=windsrcavg,
        qtsrchgt=qtsrchgt,
        qtsrc_fac=qtsrc_fac,
        thlsrc_fac=thlsrc_fac,
        frc_rasn=frc_rasn,
        rbuoy=rbuoy,
        epsvarw=epsvarw,
        use_CINcin=use_CINcin,
        mumin1=mumin1,
        rmaxfrac=rmaxfrac,
        PGFc=PGFc,
        dt=dt,
        niter_xc=niter_xc,
        criqc=criqc,
        rle=rle,
        cridist_opt=cridist_opt,
        mixscale=mixscale,
        rdrag=rdrag,
        rkm=rkm,
        use_self_detrain=use_self_detrain,
        detrhgt=detrhgt,
        use_cumpenent=use_cumpenent,
        rpen=rpen,
        use_momenflx=use_momenflx,
        rdrop=rdrop,
        iter_cin=iter_cin,
        f_pifc0_inv=pifc0_inv,
        f_zifc0_inv=zifc0_inv,
        f_pmid0_inv=pmid0_inv,
        f_zmid0_inv=zmid0_inv,
        f_kpbl_inv=kpbl_inv,
        f_exnmid0_inv=exnmid0_inv,
        f_exnifc0_inv=exnifc0_inv,
        f_dp0_inv=dp0_inv,
        f_u0_inv=u0_inv,
        f_v0_inv=v0_inv,
        f_qv0_inv=qv0_inv,
        f_ql0_inv=ql0_inv,
        f_qi0_inv=qi0_inv,
        f_t0_inv=t0_inv,
        f_frland=frland_in,
        f_tke_inv=tke_inv,
        f_rkfre=rkfre,
        f_cush=cush,
        f_shfx=shfx,
        f_evap=evap,
        f_cnvtr=cnvtr,
        f_CNV_Tracers=CNV_Tracers,
        f_umf_inv=umf_inv,
        f_dcm_inv=dcm_inv,
        f_qtflx_inv=qtflx_inv,
        f_slflx_inv=slflx_inv,
        f_uflx_inv=uflx_inv,
        f_vflx_inv=vflx_inv,
        f_qvten_inv=qvten_inv,
        f_qlten_inv=qlten_inv,
        f_qiten_inv=qiten_inv,
        f_tten_inv=tten_inv,
        f_uten_inv=uten_inv,
        f_vten_inv=vten_inv,
        f_qrten_inv=qrten_inv,
        f_qsten_inv=qsten_inv,
        f_cufrc_inv=cufrc_inv,
        f_fer_inv=fer_inv,
        f_fdr_inv=fdr_inv,
        f_ndrop_inv=ndrop_inv,
        f_nice_inv=nice_inv,
        f_qldet_inv=qldet_inv,
        f_qlsub_inv=qlsub_inv,
        f_qidet_inv=qidet_inv,
        f_qisub_inv=qisub_inv,
        f_tpert_out=tpert_out,
        f_qpert_out=qpert_out,
    )
