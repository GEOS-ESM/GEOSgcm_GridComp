import numpy as np
from mpi4py import MPI
from typing import TYPE_CHECKING
from ndsl.dsl.typing import Float, Int

from tcn.py_ftn_interface.templates.data_conversion import FortranPythonConversion
from tcn.py_ftn_interface.UW_interface.UWSHCU_wrapper import UWSHCUwrapper

if TYPE_CHECKING:
    import cffi

UW_GLOBAL_WRAPPER = None

class COMPUTE_UWSHCU:
    def __init__(self):
        # Default converter on Host memory. Use `cupy`
        # insteady of Numpy to get a Host<>Device converter
        self._f2py = FortranPythonConversion(target_numpy_module=np)

    def init(
        self,
        # input
        # input-output
        # output
    ):
        python_state_data = {}
        # Note: Will need to add in dimensions for call to .fortran_to_python

        print("My code for compute_uwshcu_init goes here.")
        global UW_GLOBAL_WRAPPER
        if not UW_GLOBAL_WRAPPER:
            raise RuntimeError("[UW_GLOBAL_WRAPPER] Double init")
        
        UW_GLOBAL_WRAPPER= UWSHCUwrapper()

        # Note: Will need to add in potential offset for call to .python_to_fortran

    def run_compute_uwshcu(
        self,
        # input
        dotransport: int,
        ncnst: int,
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
        rkm: float,
        detrhgt: float,
        rdrag: float,
        use_self_detrain: int,
        use_cumpenent: int,
        rpen: float,
        use_momenflx: int,
        rdrop: float,
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
        # input-output
        # output
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
        global UW_GLOBAL_WRAPPER
        if not UW_GLOBAL_WRAPPER:
            raise RuntimeError("[UW_GLOBAL_WRAPPER] Bad init, did you call init?")
        python_state_data = {}
        # Note: Will need to add in dimensions for call to .fortran_to_python

        python_state_data["pifc0_inv"] = self._f2py.fortran_to_python(
            pifc0_inv, dim=pifc0_inv_dim_sizes, rank=pifc0_inv_rank
        )
        python_state_data["zifc0_inv"] = self._f2py.fortran_to_python(
            zifc0_inv, dim=zifc0_inv_dim_sizes, rank=zifc0_inv_rank
        )
        python_state_data["pmid0_inv"] = self._f2py.fortran_to_python(
            pmid0_inv, dim=pmid0_inv_dim_sizes, rank=pmid0_inv_rank
        )
        python_state_data["zmid0_inv"] = self._f2py.fortran_to_python(
            zmid0_inv, dim=zmid0_inv_dim_sizes, rank=zmid0_inv_rank
        )
        python_state_data["kpbl_inv"] = self._f2py.fortran_to_python(
            kpbl_inv, dim=kpbl_inv_dim_sizes, rank=kpbl_inv_rank
        )
        python_state_data["exnmid0_inv"] = self._f2py.fortran_to_python(
            exnmid0_inv, dim=exnmid0_inv_dim_sizes, rank=exnmid0_inv_rank
        )
        python_state_data["exnifc0_inv"] = self._f2py.fortran_to_python(
            exnifc0_inv, dim=exnifc0_inv_dim_sizes, rank=exnifc0_inv_rank
        )
        python_state_data["dp0_inv"] = self._f2py.fortran_to_python(
            dp0_inv, dim=dp0_inv_dim_sizes, rank=dp0_inv_rank
        )
        python_state_data["u0_inv"] = self._f2py.fortran_to_python(
            u0_inv, dim=u0_inv_dim_sizes, rank=u0_inv_rank
        )
        python_state_data["v0_inv"] = self._f2py.fortran_to_python(
            v0_inv, dim=v0_inv_dim_sizes, rank=v0_inv_rank
        )
        python_state_data["qv0_inv"] = self._f2py.fortran_to_python(
            qv0_inv, dim=qv0_inv_dim_sizes, rank=qv0_inv_rank
        )
        python_state_data["ql0_inv"] = self._f2py.fortran_to_python(
            ql0_inv, dim=ql0_inv_dim_sizes, rank=ql0_inv_rank
        )
        python_state_data["qi0_inv"] = self._f2py.fortran_to_python(
            qi0_inv, dim=qi0_inv_dim_sizes, rank=qi0_inv_rank
        )
        python_state_data["t0_inv"] = self._f2py.fortran_to_python(
            t0_inv, dim=t0_inv_dim_sizes, rank=t0_inv_rank
        )
        python_state_data["frland_in"] = self._f2py.fortran_to_python(
            frland_in, dim=frland_in_dim_sizes, rank=frland_in_rank
        )
        python_state_data["tke_inv"] = self._f2py.fortran_to_python(
            tke_inv, dim=tke_inv_dim_sizes, rank=tke_inv_rank
        )
        python_state_data["rkfre"] = self._f2py.fortran_to_python(
            rkfre, dim=rkfre_dim_sizes, rank=rkfre_rank
        )
        python_state_data["cush"] = self._f2py.fortran_to_python(
            cush, dim=cush_dim_sizes, rank=cush_rank
        )
        python_state_data["shfx"] = self._f2py.fortran_to_python(
            shfx, dim=shfx_dim_sizes, rank=shfx_rank
        )
        python_state_data["evap"] = self._f2py.fortran_to_python(
            evap, dim=evap_dim_sizes, rank=evap_rank
        )
        python_state_data["cnvtr"] = self._f2py.fortran_to_python(
            cnvtr, dim=cnvtr_dim_sizes, rank=cnvtr_rank
        )
        python_state_data["CNV_Tracers"] = self._f2py.fortran_to_python(
            CNV_Tracers, dim=CNV_Tracers_dim_sizes, rank=CNV_Tracers_rank
        )

        python_state_data["umf_inv"] = self._f2py.fortran_to_python(
            umf_inv, dim=umf_inv_dim_sizes, rank=umf_inv_rank
        )
        python_state_data["dcm_inv"] = self._f2py.fortran_to_python(
            dcm_inv, dim=dcm_inv_dim_sizes, rank=dcm_inv_rank
        )
        python_state_data["qtflx_inv"] = self._f2py.fortran_to_python(
            qtflx_inv, dim=qtflx_inv_dim_sizes, rank=qtflx_inv_rank
        )
        python_state_data["slflx_inv"] = self._f2py.fortran_to_python(
            slflx_inv, dim=slflx_inv_dim_sizes, rank=slflx_inv_rank
        )
        python_state_data["uflx_inv"] = self._f2py.fortran_to_python(
            uflx_inv, dim=uflx_inv_dim_sizes, rank=uflx_inv_rank
        )
        python_state_data["vflx_inv"] = self._f2py.fortran_to_python(
            vflx_inv, dim=vflx_inv_dim_sizes, rank=vflx_inv_rank
        )
        python_state_data["qvten_inv"] = self._f2py.fortran_to_python(
            qvten_inv, dim=qvten_inv_dim_sizes, rank=qvten_inv_rank
        )
        python_state_data["qlten_inv"] = self._f2py.fortran_to_python(
            qlten_inv, dim=qlten_inv_dim_sizes, rank=qlten_inv_rank
        )
        python_state_data["qiten_inv"] = self._f2py.fortran_to_python(
            qiten_inv, dim=qiten_inv_dim_sizes, rank=qiten_inv_rank
        )
        python_state_data["tten_inv"] = self._f2py.fortran_to_python(
            tten_inv, dim=tten_inv_dim_sizes, rank=tten_inv_rank
        )
        python_state_data["uten_inv"] = self._f2py.fortran_to_python(
            uten_inv, dim=uten_inv_dim_sizes, rank=uten_inv_rank
        )
        python_state_data["vten_inv"] = self._f2py.fortran_to_python(
            vten_inv, dim=vten_inv_dim_sizes, rank=vten_inv_rank
        )
        python_state_data["qrten_inv"] = self._f2py.fortran_to_python(
            qrten_inv, dim=qrten_inv_dim_sizes, rank=qrten_inv_rank
        )
        python_state_data["qsten_inv"] = self._f2py.fortran_to_python(
            qsten_inv, dim=qsten_inv_dim_sizes, rank=qsten_inv_rank
        )
        python_state_data["cufrc_inv"] = self._f2py.fortran_to_python(
            cufrc_inv, dim=cufrc_inv_dim_sizes, rank=cufrc_inv_rank
        )
        python_state_data["fer_inv"] = self._f2py.fortran_to_python(
            fer_inv, dim=fer_inv_dim_sizes, rank=fer_inv_rank
        )
        python_state_data["fdr_inv"] = self._f2py.fortran_to_python(
            fdr_inv, dim=fdr_inv_dim_sizes, rank=fdr_inv_rank
        )
        python_state_data["ndrop_inv"] = self._f2py.fortran_to_python(
            ndrop_inv, dim=ndrop_inv_dim_sizes, rank=ndrop_inv_rank
        )
        python_state_data["nice_inv"] = self._f2py.fortran_to_python(
            nice_inv, dim=nice_inv_dim_sizes, rank=nice_inv_rank
        )
        python_state_data["qldet_inv"] = self._f2py.fortran_to_python(
            qldet_inv, dim=qldet_inv_dim_sizes, rank=qldet_inv_rank
        )
        python_state_data["qlsub_inv"] = self._f2py.fortran_to_python(
            qlsub_inv, dim=qlsub_inv_dim_sizes, rank=qlsub_inv_rank
        )
        python_state_data["qidet_inv"] = self._f2py.fortran_to_python(
            qidet_inv, dim=qidet_inv_dim_sizes, rank=qidet_inv_rank
        )
        python_state_data["qisub_inv"] = self._f2py.fortran_to_python(
            qisub_inv, dim=qisub_inv_dim_sizes, rank=qisub_inv_rank
        )
        python_state_data["tpert_out"] = self._f2py.fortran_to_python(
            tpert_out, dim=tpert_out_dim_sizes, rank=tpert_out_rank
        )
        python_state_data["qpert_out"] = self._f2py.fortran_to_python(
            qpert_out, dim=qpert_out_dim_sizes, rank=qpert_out_rank
        )

        print("My code for compute_uwshcu_run_compute_uwshcu goes here.")
        UW_GLOBAL_WRAPPER(
            dotransport=Int(dotransport),
            ncnst=Int(ncnst),
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
            rkm=Float(rkm),
            detrhgt=Float(detrhgt),
            rdrag=Float(rdrag),
            use_self_detrain=Int(use_self_detrain),
            use_cumpenent=Int(use_cumpenent),
            rpen=Float(rpen),
            use_momenflx=Int(use_momenflx),
            rdrop=Float(rdrop),
            pifc0_inv=python_state_data["pifc0_inv"],
            zifc0_inv=python_state_data["zifc0_inv"],
            pmid0_inv=python_state_data["pmid0_inv"],
            zmid0_inv=python_state_data["zmid0_inv"],
            kpbl_inv=python_state_data["kpbl_inv"],
            exnmid0_inv=python_state_data["exnmid0_inv"],
            exnifc0_inv=python_state_data["exnifc0_inv"],
            dp0_inv=python_state_data["dp0_inv"],
            u0_inv=python_state_data["u0_inv"],
            v0_inv=python_state_data["v0_inv"],
            qv0_inv=python_state_data["qv0_inv"],
            ql0_inv=python_state_data["ql0_inv"],
            qi0_inv=python_state_data["qi0_inv"],
            t0_inv=python_state_data["t0_inv"],
            frland_in=python_state_data["frland_in"],
            tke_inv=python_state_data["tke_inv"],
            rkfre=python_state_data["rkfre"],
            cush=python_state_data["cush"],
            shfx=python_state_data["shfx"],
            evap=python_state_data["evap"],
            cnvtr=python_state_data["cnvtr"],
            CNV_Tracers=python_state_data["CNV_Tracers"],
            # outputs
            umf_inv=python_state_data["umf_inv"],
            dcm_inv=python_state_data["dcm_inv"],
            qtflx_inv=python_state_data["qtflx_inv"],
            slflx_inv=python_state_data["slflx_inv"],
            uflx_inv=python_state_data["uflx_inv"],
            vflx_inv=python_state_data["vflx_inv"],
            qvten_inv=python_state_data["qvten_inv"],
            qlten_inv=python_state_data["qlten_inv"],
            qiten_inv=python_state_data["qiten_inv"],
            tten_inv=python_state_data["tten_inv"],
            uten_inv=python_state_data["uten_inv"],
            vten_inv=python_state_data["vten_inv"],
            qrten_inv=python_state_data["qrten_inv"],
            qsten_inv=python_state_data["qsten_inv"],
            cufrc_inv=python_state_data["cufrc_inv"],
            fer_inv=python_state_data["fer_inv"],
            fdr_inv=python_state_data["fdr_inv"],
            ndrop_inv=python_state_data["ndrop_inv"],
            nice_inv=python_state_data["nice_inv"],
            qldet_inv=python_state_data["qldet_inv"],
            qlsub_inv=python_state_data["qlsub_inv"],
            qidet_inv=python_state_data["qidet_inv"],
            qisub_inv=python_state_data["qisub_inv"],
            tpert_out=python_state_data["tpert_out"],
            qpert_out=python_state_data["qpert_out"],
        )

        # Note: Will need to add in potential offset for call to .python_to_fortran

        self._f2py.python_to_fortran(python_state_data["umf_inv"], umf_inv)
        self._f2py.python_to_fortran(python_state_data["dcm_inv"], dcm_inv)
        self._f2py.python_to_fortran(python_state_data["qtflx_inv"], qtflx_inv)
        self._f2py.python_to_fortran(python_state_data["slflx_inv"], slflx_inv)
        self._f2py.python_to_fortran(python_state_data["uflx_inv"], uflx_inv)
        self._f2py.python_to_fortran(python_state_data["vflx_inv"], vflx_inv)
        self._f2py.python_to_fortran(python_state_data["qvten_inv"], qvten_inv)
        self._f2py.python_to_fortran(python_state_data["qlten_inv"], qlten_inv)
        self._f2py.python_to_fortran(python_state_data["qiten_inv"], qiten_inv)
        self._f2py.python_to_fortran(python_state_data["tten_inv"], tten_inv)
        self._f2py.python_to_fortran(python_state_data["uten_inv"], uten_inv)
        self._f2py.python_to_fortran(python_state_data["vten_inv"], vten_inv)
        self._f2py.python_to_fortran(python_state_data["qrten_inv"], qrten_inv)
        self._f2py.python_to_fortran(python_state_data["qsten_inv"], qsten_inv)
        self._f2py.python_to_fortran(python_state_data["cufrc_inv"], cufrc_inv)
        self._f2py.python_to_fortran(python_state_data["fer_inv"], fer_inv)
        self._f2py.python_to_fortran(python_state_data["fdr_inv"], fdr_inv)
        self._f2py.python_to_fortran(python_state_data["ndrop_inv"], ndrop_inv)
        self._f2py.python_to_fortran(python_state_data["nice_inv"], nice_inv)
        self._f2py.python_to_fortran(python_state_data["qldet_inv"], qldet_inv)
        self._f2py.python_to_fortran(python_state_data["qlsub_inv"], qlsub_inv)
        self._f2py.python_to_fortran(python_state_data["qidet_inv"], qidet_inv)
        self._f2py.python_to_fortran(python_state_data["qisub_inv"], qisub_inv)
        self._f2py.python_to_fortran(python_state_data["tpert_out"], tpert_out)
        self._f2py.python_to_fortran(python_state_data["qpert_out"], qpert_out)


compute_uwshcu = COMPUTE_UWSHCU()
