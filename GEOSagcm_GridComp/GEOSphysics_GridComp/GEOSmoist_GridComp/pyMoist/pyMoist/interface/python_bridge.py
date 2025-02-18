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

    def GFDL_1M_driver(
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
            self.pymoist.GFDL_1M_driver(
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
                uin=uin,
                vin=vin,
                udt=udt,
                vdt=vdt,
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
                self.pymoist.GFDL_1M_driver.rain.view[:], f_rain
            )
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.snow.view[:], f_snow
            )
            self.f_py.python_to_fortran(self.pymoist.GFDL_1M_driver.ice.view[:], f_ice)
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.graupel.view[:], f_graupel
            )
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.m2_rain.view[:], f_m2_rain
            )
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.m2_sol.view[:], f_m2_sol
            )
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.revap.view[:], f_revap
            )
            self.f_py.python_to_fortran(
                self.pymoist.GFDL_1M_driver.isubl.view[:], f_isubl
            )


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

    WRAPPER.GFDL_1M_driver(
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
        raise RuntimeError("[PYMOIST WRAPPER] Double init")
    WRAPPER.init(
        pyMoist_flags=pyMoist_flags,
        backend=BACKEND,
    )


def gfdl_1m_init(gfdl_1m_flags: cffi.FFI.CData):
    if not WRAPPER.ready:
        raise RuntimeError("[GFDL_1M WRAPPER] pyMoist_init needs to be called first")
    WRAPPER.pymoist.init_gfdl_1m_flags(
        flags=gfdl_1m_flags_f_to_python(gfdl_1m_flags),
    )
