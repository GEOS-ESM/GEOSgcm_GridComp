from __future__ import annotations

import os
from typing import Dict, List

import cffi
import numpy as np
from mpi4py import MPI

from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import Float, Int
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
                dotransport=Int(dotransport),
                k0=Int(k0),
                windsrcavg=Int(windsrcavg),
                qtsrchgt=Float(qtsrchgt),
                qtsrc_fac=Float(qtsrc_fac),
                thlsrc_fac=Float(thlsrc_fac),
                frc_rasn=Float(frc_rasn),
                rbuoy=Float(rbuoy),
                epsvarw=Float(epsvarw),
                use_CINcin=Int(use_CINcin),
                mumin1=Float(mumin1),
                rmaxfrac=Float(rmaxfrac),
                PGFc=Float(PGFc),
                dt=Float(dt),
                niter_xc=Int(niter_xc),
                criqc=Float(criqc),
                rle=Float(rle),
                cridist_opt=Int(cridist_opt),
                mixscale=Float(mixscale),
                rdrag=Float(rdrag),
                rkm=Float(rkm),
                use_self_detrain=Int(use_self_detrain),
                detrhgt=Float(detrhgt),
                use_cumpenent=Int(use_cumpenent),
                rpen=Float(rpen),
                use_momenflx=Int(use_momenflx),
                rdrop=Float(rdrop),
                iter_cin=Int(iter_cin),
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
    mapl_comp: cffi.FFI.CData,
    pyMoist_flags: cffi.FFI.CData,
):
    # Read in the backend
    BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "dace:cpu")
    if WRAPPER.ready:
        raise RuntimeError("[GEOS WRAPPER] Bad init, did you call init?")
    WRAPPER.init(
        mapl_states=MAPLStates(
            import_=import_state,
            export=export_state,
            mapl_comp=mapl_comp,
        ),
        pyMoist_flags=pyMoist_flags,
        backend=BACKEND,
    )


def gfdl_1m_init(
    gfdl_1m_flags: cffi.FFI.CData,
    internal_state: cffi.FFI.CData,
) -> None:
    if not WRAPPER.ready:
        raise RuntimeError("[GFDL_1M WRAPPER] pyMoist_init needs to be called first")
    if not WRAPPER.pymoist._GFDL_1M_ready:
        WRAPPER.pymoist.init_gfdl_1m_configuration(
            flags=gfdl_1m_flags_f_to_python(gfdl_1m_flags),
            internal_state=internal_state,
        )
        WRAPPER.pymoist._GFDL_1M_ready = True


def compute_uwshcu_init(
    ncnst: int,
    k0: int,
    windsrcavg: int,
):
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
