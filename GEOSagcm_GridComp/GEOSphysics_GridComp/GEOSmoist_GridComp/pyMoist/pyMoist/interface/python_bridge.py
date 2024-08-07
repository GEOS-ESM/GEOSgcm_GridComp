import os

# from f_py_conversion import FortranPythonConversion
# from cuda_profiler import CUDAProfiler, TimedCUDAProfiler
from mpi4py import MPI
from ndsl.optional_imports import cupy as cp
import numpy as np
from ndsl.dsl.gt4py_utils import is_gpu_backend
from typing import TYPE_CHECKING
from pyMoist.wrapper import GEOSPyMoistWrapper, MoistFlags


if TYPE_CHECKING:
    import cffi


class PYMOIST_WRAPPER:
    def __init__(
        self,
        flags: MoistFlags,
        backend: str = "dace:gpu",
    ) -> None:
        print("Wrapper.init")
        # self.rank = comm.Get_rank()
        # self.backend = backend
        # # For Fortran<->NumPy conversion
        # if is_gpu_backend(self.backend):
        #     numpy_module = cp
        #     fortran_mem_space = MemorySpace.DEVICE
        # else:
        #     numpy_module = np
        #     fortran_mem_space = MemorySpace.HOST
        # self.f_py = FortranPythonConversion(
        #     npx,
        #     npy,
        #     npz,
        #     is_,
        #     ie,
        #     js,
        #     je,
        #     isd,
        #     ied,
        #     jsd,
        #     jed,
        #     tracer_count,
        #     numpy_module,
        # )

        # # Input pressure levels
        # ak = self.f_py._fortran_to_numpy(ak_cdata, [npz + 1])
        # bk = self.f_py._fortran_to_numpy(bk_cdata, [npz + 1])

        # # Setup pyFV3's dynamical core
        # self.dycore = GeosDycoreWrapper(
        #     fv_flags=fv_flags,
        #     bdt=bdt,
        #     comm=comm,
        #     ak=ak,
        #     bk=bk,
        #     backend=self.backend,
        #     tracer_count=tracer_count,
        #     fortran_mem_space=fortran_mem_space,
        # )

        # self._timings = {}

    def finalize(self):
        # import json

        # with open("pyfv3_timings.json", "w") as f:
        #     json.dump(self._timings, f, indent=4)
        print("Wrapper.finalize")

    def __call__(
        self,
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
        print("Wrapper.__call__")
        # CUDAProfiler.start_cuda_profiler()
        # with TimedCUDAProfiler("Fortran -> Python", self._timings):
        #     # Convert Fortran arrays to NumPy
        #     state_in = self.f_py.fortran_to_python(
        #         # input
        #         u,
        #         v,
        #         w,
        #         delz,
        #         pt,
        #         delp,
        #         q,
        #         ps,
        #         pe,
        #         pk,
        #         peln,
        #         pkz,
        #         phis,
        #         q_con,
        #         omga,
        #         ua,
        #         va,
        #         uc,
        #         vc,
        #         mfx,
        #         mfy,
        #         cx,
        #         cy,
        #         diss_est,
        #     )

        # # Run pyFV3
        # with TimedCUDAProfiler("Numerics", self._timings):
        #     state_out, self._timings = self.dycore(
        #         self._timings,
        #         state_in["u"],
        #         state_in["v"],
        #         state_in["w"],
        #         state_in["delz"],
        #         state_in["pt"],
        #         state_in["delp"],
        #         state_in["q"],
        #         state_in["ps"],
        #         state_in["pe"],
        #         state_in["pk"],
        #         state_in["peln"],
        #         state_in["pkz"],
        #         state_in["phis"],
        #         state_in["q_con"],
        #         state_in["omga"],
        #         state_in["ua"],
        #         state_in["va"],
        #         state_in["uc"],
        #         state_in["vc"],
        #         state_in["mfxd"],
        #         state_in["mfyd"],
        #         state_in["cxd"],
        #         state_in["cyd"],
        #         state_in["diss_estd"],
        #     )

        # # Convert NumPy arrays back to Fortran
        # with TimedCUDAProfiler("Python -> Fortran", self._timings):
        #     self.f_py.python_to_fortran(
        #         # input
        #         state_out,
        #         # output
        #         u,
        #         v,
        #         w,
        #         delz,
        #         pt,
        #         delp,
        #         q,
        #         ps,
        #         pe,
        #         pk,
        #         peln,
        #         pkz,
        #         phis,
        #         q_con,
        #         omga,
        #         ua,
        #         va,
        #         uc,
        #         vc,
        #         mfx,
        #         mfy,
        #         cx,
        #         cy,
        #         diss_est,
        #     )


# Below is the entry point to the interface
# ToDo: we should build the object outside of the sim loop from fortran
# potentially by writing a pyfv3_interface_setup and caching the ptr Fortran side
# or by having a central python interpreter object handled by CFFI to register against
WRAPPER = None


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
    global WRAPPER
    if not WRAPPER:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER(
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
    if WRAPPER is not None:
        WRAPPER.finalize()


def pyMoist_init(flags: MoistFlags):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "gt:gpu")

    global WRAPPER
    if WRAPPER is not None:
        raise RuntimeError("[PYMOIST WRAPPER] Double init")
    WRAPPER = PYMOIST_WRAPPER(
        flags=flags,
        backend=BACKEND,
    )
