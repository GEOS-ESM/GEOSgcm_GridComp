import os
from typing import TYPE_CHECKING, Dict, List

import numpy as np
from mpi4py import MPI

from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import Float
from ndsl.optional_imports import cupy as cp
from pyMoist.interface.cuda_profiler import CUDAProfiler, TimedCUDAProfiler
from pyMoist.interface.f_py_conversion import FortranPythonConversion
from pyMoist.interface.flags import flags_fv_to_python
from pyMoist.interface.wrapper import GEOSPyMoistWrapper, MemorySpace


if TYPE_CHECKING:
    import cffi


class PYMOIST_WRAPPER:
    def __init__(self) -> None:
        self.ready = False

    def init(
        self,
        fv_flags: "cffi.FFI.CData",
        backend: str = "dace:cpu",
    ) -> None:
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.backend = backend
        self.flags = flags_fv_to_python(fv_flags)
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

        # Setup pyFV3's dynamical core
        self.pymoist = GEOSPyMoistWrapper(self.flags, backend)

        self._timings: Dict[str, List[float]] = {}
        self.ready = True

    def finalize(self):
        import json

        with open("pymoist_timings.json", "w") as f:
            json.dump(self._timings, f, indent=4)

    def aer_activation(
        self,
        f_aero_dgn: "cffi.FFI.CData",
        f_aero_num: "cffi.FFI.CData",
        f_aero_hygroscopicity: "cffi.FFI.CData",
        f_aero_sigma: "cffi.FFI.CData",
        f_frland: "cffi.FFI.CData",
        f_nn_ocean: np.float32,
        f_nn_land: np.float32,
        f_t: "cffi.FFI.CData",
        f_plo: "cffi.FFI.CData",
        f_qicn: "cffi.FFI.CData",
        f_qils: "cffi.FFI.CData",
        f_qlcn: "cffi.FFI.CData",
        f_qlls: "cffi.FFI.CData",
        f_vvel: "cffi.FFI.CData",
        f_tke: "cffi.FFI.CData",
        f_nacti: "cffi.FFI.CData",
        f_nwfa: "cffi.FFI.CData",
        f_nactl: "cffi.FFI.CData",
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

    def GFDL_1M__EVAP_SUBL_HYSTPDF(
        self,
        f_EIS: "cffi.FFI.CData",
        f_dw_land: np.float32,
        f_dw_ocean: np.float32,
        f_PDF_shape: np.int32,
        f_TurnRHCrit_Param: "cffi.FFI.CData",
        f_PLmb: "cffi.FFI.CData",
        f_KLCL: "cffi.FFI.CData",
        f_PLEmb: "cffi.FFI.CData",
        f_Area: "cffi.FFI.CData",
        f_dt_moist: np.float32,
        f_CNV_FRC: "cffi.FFI.CData",
        f_SRF_TYPE: "cffi.FFI.CData",
        f_T: "cffi.FFI.CData",
        f_QLCN: "cffi.FFI.CData",
        f_QICN: "cffi.FFI.CData",
        f_QLLS: "cffi.FFI.CData",
        f_QILS: "cffi.FFI.CData",
        f_CCW_EVAP_EFF: np.float32,
        f_CCI_EVAP_EFF: np.float32,
        f_Q: "cffi.FFI.CData",
        f_CLLS: "cffi.FFI.CData",
        f_CLCN: "cffi.FFI.CData",
        f_NACTL: "cffi.FFI.CData",
        f_NACTI: "cffi.FFI.CData",
        f_QST: "cffi.FFI.CData",
        f_SUBLC: "cffi.FFI.CData",
        f_EVAPC: "cffi.FFI.CData",
        f_RHX: "cffi.FFI.CData",
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
            self.pymoist.gfdl_1M(
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


WRAPPER = PYMOIST_WRAPPER()


def pyMoist_run_GFDL1M(
    dw_land: np.float32,
    dw_ocean: np.float32,
    PDFSHAPE: np.int32,
    TURNRHCRIT_PARAM: np.float32,
    DT_MOIST: np.float32,
    CCW_EVAP_EFF: np.float32,
    CCI_EVAP_EFF: np.float32,
    LMELTFRZ: bool,
    AREA: "cffi.FFI.CData",
    CNV_FRC: "cffi.FFI.CData",
    SRF_TYPE: "cffi.FFI.CData",
    KLCL: "cffi.FFI.CData",
    EIS: "cffi.FFI.CData",
    PLmb: "cffi.FFI.CData",
    PLEmb: "cffi.FFI.CData",
    NACTL: "cffi.FFI.CData",
    NACTI: "cffi.FFI.CData",
    QST: "cffi.FFI.CData",
    T: "cffi.FFI.CData",
    Q: "cffi.FFI.CData",
    QLCN: "cffi.FFI.CData",
    QICN: "cffi.FFI.CData",
    QLLS: "cffi.FFI.CData",
    QILS: "cffi.FFI.CData",
    CLLS: "cffi.FFI.CData",
    CLCN: "cffi.FFI.CData",
    SUBLC: "cffi.FFI.CData",
    EVAPC: "cffi.FFI.CData",
    RHX: "cffi.FFI.CData",
):
    if not WRAPPER.ready:
        raise RuntimeError("[pyMoist] Bad init, did you call init?")
    WRAPPER.GFDL_1M__EVAP_SUBL_HYSTPDF(
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


def pyMoist_run_AerActivation(
    aero_dgn: "cffi.FFI.CData",
    aero_num: "cffi.FFI.CData",
    aero_hygroscopicity: "cffi.FFI.CData",
    aero_sigma: "cffi.FFI.CData",
    frland: "cffi.FFI.CData",
    nn_ocean: np.float32,
    nn_land: np.float32,
    t: "cffi.FFI.CData",
    plo: "cffi.FFI.CData",
    qicn: "cffi.FFI.CData",
    qils: "cffi.FFI.CData",
    qlcn: "cffi.FFI.CData",
    qlls: "cffi.FFI.CData",
    vvel: "cffi.FFI.CData",
    tke: "cffi.FFI.CData",
    nacti: "cffi.FFI.CData",
    nwfa: "cffi.FFI.CData",
    nactl: "cffi.FFI.CData",
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


def pyMoist_init(fv_flags: "cffi.FFI.CData"):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "dace:cpu")
    if WRAPPER.ready:
        raise RuntimeError("[PYMOIST WRAPPER] Double init")
    WRAPPER.init(
        fv_flags=fv_flags,
        backend=BACKEND,
    )
