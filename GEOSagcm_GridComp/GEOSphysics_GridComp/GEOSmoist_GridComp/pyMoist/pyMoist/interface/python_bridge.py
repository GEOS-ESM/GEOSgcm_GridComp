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
from pyMoist.interface.wrapper import GEOSPyMoistWrapper, MAPLStates, MemorySpace


class PYMOIST_WRAPPER:
    def __init__(self) -> None:
        self.ready = False

    def init(
        self,
        mapl_states: MAPLStates,
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
        self.pymoist = GEOSPyMoistWrapper(mapl_states, self.flags, backend)

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

            frland = self.f_py.fortran_to_python(f_frland, [self.flags.npx, self.flags.npy])

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


WRAPPER = PYMOIST_WRAPPER()


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


def pymoist_interface_GFDL_1M():
    if WRAPPER.ready:
        WRAPPER.pymoist.GFDL_1M_Microphysics()


def pyMoist_finalize():
    if WRAPPER.ready:
        WRAPPER.finalize()


def pyMoist_init(
    import_state: cffi.FFI.CData,
    export_state: cffi.FFI.CData,
    internal_state: cffi.FFI.CData,
    mapl_comp: cffi.FFI.CData,
    pyMoist_flags: cffi.FFI.CData,
):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "dace:cpu")
    if WRAPPER.ready:
        raise RuntimeError("[PYMOIST WRAPPER] Double init")
    WRAPPER.init(
        mapl_states=MAPLStates(import_state, export_state, internal_state, mapl_comp),
        pyMoist_flags=pyMoist_flags,
        backend=BACKEND,
    )


def gfdl_1m_init(gfdl_1m_flags: cffi.FFI.CData) -> None:
    if not WRAPPER.ready:
        raise RuntimeError("[GFDL_1M WRAPPER] pyMoist_init needs to be called first")
    if not WRAPPER.pymoist._GFDL_1M_ready:
        WRAPPER.pymoist.init_gfdl_1m_configuration(
            flags=gfdl_1m_flags_f_to_python(gfdl_1m_flags),
        )
        WRAPPER.pymoist._GFDL_1M_ready = True
