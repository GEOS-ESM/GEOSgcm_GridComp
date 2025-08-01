from distutils.sysconfig import get_config_var

import cffi
from mpi4py import MPI


TMPFILEBASE = "pyMoist_interface_py"

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof("int"):
    _mpi_comm_t = "int"
else:
    _mpi_comm_t = "void*"

source = """
from {} import ffi
from datetime import datetime
from mpi4py import MPI
from pyMoist.interface.python_bridge import (
    pyMoist_init,
    pyMoist_run_AerActivation,
    pyMoist_run_GFDL_1M_evap_subl_hystpdf,
    pymoist_run_GFDL_1M_driver,
    pyMoist_finalize,
    gfdl_1m_init,
    compute_uwshcu_init,
    compute_uwshcu_run,
)
import traceback

@ffi.def_extern()
def pymoist_interface_py_init(flags) -> int:

    try:
        pyMoist_init(flags)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def gfdl_1m_interface_py_init(flags) -> int:

    try:
        gfdl_1m_init(flags)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_run_AerActivation(
    aero_dgn, aero_num, aero_hygroscopicity, aero_sigma,
    frland, nn_ocean, nn_land,
    t, plo,
    qicn, qils, qlcn, qlls,
    vvel, tke,
    nacti, nwfa, nactl) -> int:

    try:
        pyMoist_run_AerActivation(
            aero_dgn, aero_num, aero_hygroscopicity, aero_sigma,
            frland, nn_ocean, nn_land,
            t, plo,
            qicn, qils, qlcn, qlls,
            vvel, tke,
            nacti, nwfa, nactl)
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_run_GFDL1M(
    dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
    DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
    LMELTFRZ,
    AREA, CNV_FRC, SRF_TYPE,
    KLCL,
    EIS, PLmb, PLEmb, NACTL, NACTI, QST,
    T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
    SUBLC, EVAPC, RHX):

    try:
        pyMoist_run_GFDL_1M_evap_subl_hystpdf(
            dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
            DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
            LMELTFRZ,
            AREA, CNV_FRC, SRF_TYPE,
            KLCL,
            EIS, PLmb, PLEmb, NACTL, NACTI, QST,
            T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
            SUBLC, EVAPC, RHX
        )
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def pymoist_interface_py_run_GFDL_1M_driver(
        RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll,
        DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic,
        DQSDTmic, DQGDTmic, DQADTmic, DTDTmic,
        T, W, U, V, DUDTmic, DVDTmic, DZ, DP,
        AREA, FRLAND, CNV_FRC, SRF_TYPE, EIS, RHCRIT3D,
        DT_MOIST, ANV_ICEFALL, LS_ICEFALL,
        REV_LS, RSU_LS,
        PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, PFL_LS, PFI_LS,
        LHYDROSTATIC, LPHYS_HYDROSTATIC
    ):

        try:
            pymoist_run_GFDL_1M_driver(
                RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF,
                NACTAll,
                DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic,
                DQSDTmic, DQGDTmic, DQADTmic, DTDTmic,
                T, W, U, V, DUDTmic, DVDTmic, DZ, DP,
                AREA, FRLAND, CNV_FRC, SRF_TYPE, EIS, RHCRIT3D,
                DT_MOIST, ANV_ICEFALL, LS_ICEFALL,
                REV_LS, RSU_LS,
                PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, PFL_LS, PFI_LS,
                LHYDROSTATIC, LPHYS_HYDROSTATIC
            )
        except Exception as err:
            print("Error in Python:")
            print(traceback.format_exc())
            return -1
        return 0

@ffi.def_extern()
def pymoist_interface_py_finalize() -> int:
    try:
        pyMoist_finalize()
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def compute_uwshcu_py_init(
    # inputs
    ncnst:int,
    k0:int,
    windsrcavg:int,
    ) -> int:

    # Transform init code
    try:
        compute_uwshcu_init(
            ncnst=ncnst,
            k0=k0,
            windsrcavg=windsrcavg,
        )
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0

@ffi.def_extern()
def compute_uwshcu_py_run(
        # inputs
        dotransport:int,
        k0:int,
        windsrcavg:int,
        qtsrchgt:float,
        qtsrc_fac:float,
        thlsrc_fac:float,
        frc_rasn:float,
        rbuoy:float,
        epsvarw:float,
        use_CINcin:int,
        mumin1:float,
        rmaxfrac:float,
        PGFc:float,
        dt:float,
        niter_xc:int,
        criqc:float,
        rle:float,
        cridist_opt:int,
        mixscale:float,
        rdrag:float,
        rkm:float,
        use_self_detrain:int,
        detrhgt:float,
        use_cumpenent:int,
        rpen:float,
        use_momenflx:int,
        rdrop:float,
        iter_cin:int,
        pifc0_inv:'cffi.FFI.CData',
        pifc0_inv_dim_sizes:'cffi.FFI.CData',
        pifc0_inv_rank:int,
        zifc0_inv:'cffi.FFI.CData',
        zifc0_inv_dim_sizes:'cffi.FFI.CData',
        zifc0_inv_rank:int,
        pmid0_inv:'cffi.FFI.CData',
        pmid0_inv_dim_sizes:'cffi.FFI.CData',
        pmid0_inv_rank:int,
        zmid0_inv:'cffi.FFI.CData',
        zmid0_inv_dim_sizes:'cffi.FFI.CData',
        zmid0_inv_rank:int,
        kpbl_inv:'cffi.FFI.CData',
        kpbl_inv_dim_sizes:'cffi.FFI.CData',
        kpbl_inv_rank:int,
        exnmid0_inv:'cffi.FFI.CData',
        exnmid0_inv_dim_sizes:'cffi.FFI.CData',
        exnmid0_inv_rank:int,
        exnifc0_inv:'cffi.FFI.CData',
        exnifc0_inv_dim_sizes:'cffi.FFI.CData',
        exnifc0_inv_rank:int,
        dp0_inv:'cffi.FFI.CData',
        dp0_inv_dim_sizes:'cffi.FFI.CData',
        dp0_inv_rank:int,
        u0_inv:'cffi.FFI.CData',
        u0_inv_dim_sizes:'cffi.FFI.CData',
        u0_inv_rank:int,
        v0_inv:'cffi.FFI.CData',
        v0_inv_dim_sizes:'cffi.FFI.CData',
        v0_inv_rank:int,
        qv0_inv:'cffi.FFI.CData',
        qv0_inv_dim_sizes:'cffi.FFI.CData',
        qv0_inv_rank:int,
        ql0_inv:'cffi.FFI.CData',
        ql0_inv_dim_sizes:'cffi.FFI.CData',
        ql0_inv_rank:int,
        qi0_inv:'cffi.FFI.CData',
        qi0_inv_dim_sizes:'cffi.FFI.CData',
        qi0_inv_rank:int,
        t0_inv:'cffi.FFI.CData',
        t0_inv_dim_sizes:'cffi.FFI.CData',
        t0_inv_rank:int,
        frland_in:'cffi.FFI.CData',
        frland_in_dim_sizes:'cffi.FFI.CData',
        frland_in_rank:int,
        tke_inv:'cffi.FFI.CData',
        tke_inv_dim_sizes:'cffi.FFI.CData',
        tke_inv_rank:int,
        rkfre:'cffi.FFI.CData',
        rkfre_dim_sizes:'cffi.FFI.CData',
        rkfre_rank:int,
        cush:'cffi.FFI.CData',
        cush_dim_sizes:'cffi.FFI.CData',
        cush_rank:int,
        shfx:'cffi.FFI.CData',
        shfx_dim_sizes:'cffi.FFI.CData',
        shfx_rank:int,
        evap:'cffi.FFI.CData',
        evap_dim_sizes:'cffi.FFI.CData',
        evap_rank:int,
        cnvtr:'cffi.FFI.CData',
        cnvtr_dim_sizes:'cffi.FFI.CData',
        cnvtr_rank:int,
        CNV_Tracers:'cffi.FFI.CData',
        CNV_Tracers_dim_sizes:'cffi.FFI.CData',
        CNV_Tracers_rank:int,

        # inputs-outputs

        # outputs
        umf_inv:'cffi.FFI.CData',
        umf_inv_dim_sizes:'cffi.FFI.CData',
        umf_inv_rank:int,
        dcm_inv:'cffi.FFI.CData',
        dcm_inv_dim_sizes:'cffi.FFI.CData',
        dcm_inv_rank:int,
        qtflx_inv:'cffi.FFI.CData',
        qtflx_inv_dim_sizes:'cffi.FFI.CData',
        qtflx_inv_rank:int,
        slflx_inv:'cffi.FFI.CData',
        slflx_inv_dim_sizes:'cffi.FFI.CData',
        slflx_inv_rank:int,
        uflx_inv:'cffi.FFI.CData',
        uflx_inv_dim_sizes:'cffi.FFI.CData',
        uflx_inv_rank:int,
        vflx_inv:'cffi.FFI.CData',
        vflx_inv_dim_sizes:'cffi.FFI.CData',
        vflx_inv_rank:int,
        qvten_inv:'cffi.FFI.CData',
        qvten_inv_dim_sizes:'cffi.FFI.CData',
        qvten_inv_rank:int,
        qlten_inv:'cffi.FFI.CData',
        qlten_inv_dim_sizes:'cffi.FFI.CData',
        qlten_inv_rank:int,
        qiten_inv:'cffi.FFI.CData',
        qiten_inv_dim_sizes:'cffi.FFI.CData',
        qiten_inv_rank:int,
        tten_inv:'cffi.FFI.CData',
        tten_inv_dim_sizes:'cffi.FFI.CData',
        tten_inv_rank:int,
        uten_inv:'cffi.FFI.CData',
        uten_inv_dim_sizes:'cffi.FFI.CData',
        uten_inv_rank:int,
        vten_inv:'cffi.FFI.CData',
        vten_inv_dim_sizes:'cffi.FFI.CData',
        vten_inv_rank:int,
        qrten_inv:'cffi.FFI.CData',
        qrten_inv_dim_sizes:'cffi.FFI.CData',
        qrten_inv_rank:int,
        qsten_inv:'cffi.FFI.CData',
        qsten_inv_dim_sizes:'cffi.FFI.CData',
        qsten_inv_rank:int,
        cufrc_inv:'cffi.FFI.CData',
        cufrc_inv_dim_sizes:'cffi.FFI.CData',
        cufrc_inv_rank:int,
        fer_inv:'cffi.FFI.CData',
        fer_inv_dim_sizes:'cffi.FFI.CData',
        fer_inv_rank:int,
        fdr_inv:'cffi.FFI.CData',
        fdr_inv_dim_sizes:'cffi.FFI.CData',
        fdr_inv_rank:int,
        ndrop_inv:'cffi.FFI.CData',
        ndrop_inv_dim_sizes:'cffi.FFI.CData',
        ndrop_inv_rank:int,
        nice_inv:'cffi.FFI.CData',
        nice_inv_dim_sizes:'cffi.FFI.CData',
        nice_inv_rank:int,
        qldet_inv:'cffi.FFI.CData',
        qldet_inv_dim_sizes:'cffi.FFI.CData',
        qldet_inv_rank:int,
        qlsub_inv:'cffi.FFI.CData',
        qlsub_inv_dim_sizes:'cffi.FFI.CData',
        qlsub_inv_rank:int,
        qidet_inv:'cffi.FFI.CData',
        qidet_inv_dim_sizes:'cffi.FFI.CData',
        qidet_inv_rank:int,
        qisub_inv:'cffi.FFI.CData',
        qisub_inv_dim_sizes:'cffi.FFI.CData',
        qisub_inv_rank:int,
        tpert_out:'cffi.FFI.CData',
        tpert_out_dim_sizes:'cffi.FFI.CData',
        tpert_out_rank:int,
        qpert_out:'cffi.FFI.CData',
        qpert_out_dim_sizes:'cffi.FFI.CData',
        qpert_out_rank:int,
    ) -> int:
    # Transform init code
    try:
        compute_uwshcu_run(
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
            pifc0_inv=pifc0_inv,
            pifc0_inv_dim_sizes=pifc0_inv_dim_sizes,
            pifc0_inv_rank=pifc0_inv_rank,
            zifc0_inv=zifc0_inv,
            zifc0_inv_dim_sizes=zifc0_inv_dim_sizes,
            zifc0_inv_rank=zifc0_inv_rank,
            pmid0_inv=pmid0_inv,
            pmid0_inv_dim_sizes=pmid0_inv_dim_sizes,
            pmid0_inv_rank=pmid0_inv_rank,
            zmid0_inv=zmid0_inv,
            zmid0_inv_dim_sizes=zmid0_inv_dim_sizes,
            zmid0_inv_rank=zmid0_inv_rank,
            kpbl_inv=kpbl_inv,
            kpbl_inv_dim_sizes=kpbl_inv_dim_sizes,
            kpbl_inv_rank=kpbl_inv_rank,
            exnmid0_inv=exnmid0_inv,
            exnmid0_inv_dim_sizes=exnmid0_inv_dim_sizes,
            exnmid0_inv_rank=exnmid0_inv_rank,
            exnifc0_inv=exnifc0_inv,
            exnifc0_inv_dim_sizes=exnifc0_inv_dim_sizes,
            exnifc0_inv_rank=exnifc0_inv_rank,
            dp0_inv=dp0_inv,
            dp0_inv_dim_sizes=dp0_inv_dim_sizes,
            dp0_inv_rank=dp0_inv_rank,
            u0_inv=u0_inv,
            u0_inv_dim_sizes=u0_inv_dim_sizes,
            u0_inv_rank=u0_inv_rank,
            v0_inv=v0_inv,
            v0_inv_dim_sizes=v0_inv_dim_sizes,
            v0_inv_rank=v0_inv_rank,
            qv0_inv=qv0_inv,
            qv0_inv_dim_sizes=qv0_inv_dim_sizes,
            qv0_inv_rank=qv0_inv_rank,
            ql0_inv=ql0_inv,
            ql0_inv_dim_sizes=ql0_inv_dim_sizes,
            ql0_inv_rank=ql0_inv_rank,
            qi0_inv=qi0_inv,
            qi0_inv_dim_sizes=qi0_inv_dim_sizes,
            qi0_inv_rank=qi0_inv_rank,
            t0_inv=t0_inv,
            t0_inv_dim_sizes=t0_inv_dim_sizes,
            t0_inv_rank=t0_inv_rank,
            frland_in=frland_in,
            frland_in_dim_sizes=frland_in_dim_sizes,
            frland_in_rank=frland_in_rank,
            tke_inv=tke_inv,
            tke_inv_dim_sizes=tke_inv_dim_sizes,
            tke_inv_rank=tke_inv_rank,
            rkfre=rkfre,
            rkfre_dim_sizes=rkfre_dim_sizes,
            rkfre_rank=rkfre_rank,
            cush=cush,
            cush_dim_sizes=cush_dim_sizes,
            cush_rank=cush_rank,
            shfx=shfx,
            shfx_dim_sizes=shfx_dim_sizes,
            shfx_rank=shfx_rank,
            evap=evap,
            evap_dim_sizes=evap_dim_sizes,
            evap_rank=evap_rank,
            cnvtr=cnvtr,
            cnvtr_dim_sizes=cnvtr_dim_sizes,
            cnvtr_rank=cnvtr_rank,
            CNV_Tracers=CNV_Tracers,
            CNV_Tracers_dim_sizes=CNV_Tracers_dim_sizes,
            CNV_Tracers_rank=CNV_Tracers_rank,


            umf_inv=umf_inv,
            umf_inv_dim_sizes=umf_inv_dim_sizes,
            umf_inv_rank=umf_inv_rank,
            dcm_inv=dcm_inv,
            dcm_inv_dim_sizes=dcm_inv_dim_sizes,
            dcm_inv_rank=dcm_inv_rank,
            qtflx_inv=qtflx_inv,
            qtflx_inv_dim_sizes=qtflx_inv_dim_sizes,
            qtflx_inv_rank=qtflx_inv_rank,
            slflx_inv=slflx_inv,
            slflx_inv_dim_sizes=slflx_inv_dim_sizes,
            slflx_inv_rank=slflx_inv_rank,
            uflx_inv=uflx_inv,
            uflx_inv_dim_sizes=uflx_inv_dim_sizes,
            uflx_inv_rank=uflx_inv_rank,
            vflx_inv=vflx_inv,
            vflx_inv_dim_sizes=vflx_inv_dim_sizes,
            vflx_inv_rank=vflx_inv_rank,
            qvten_inv=qvten_inv,
            qvten_inv_dim_sizes=qvten_inv_dim_sizes,
            qvten_inv_rank=qvten_inv_rank,
            qlten_inv=qlten_inv,
            qlten_inv_dim_sizes=qlten_inv_dim_sizes,
            qlten_inv_rank=qlten_inv_rank,
            qiten_inv=qiten_inv,
            qiten_inv_dim_sizes=qiten_inv_dim_sizes,
            qiten_inv_rank=qiten_inv_rank,
            tten_inv=tten_inv,
            tten_inv_dim_sizes=tten_inv_dim_sizes,
            tten_inv_rank=tten_inv_rank,
            uten_inv=uten_inv,
            uten_inv_dim_sizes=uten_inv_dim_sizes,
            uten_inv_rank=uten_inv_rank,
            vten_inv=vten_inv,
            vten_inv_dim_sizes=vten_inv_dim_sizes,
            vten_inv_rank=vten_inv_rank,
            qrten_inv=qrten_inv,
            qrten_inv_dim_sizes=qrten_inv_dim_sizes,
            qrten_inv_rank=qrten_inv_rank,
            qsten_inv=qsten_inv,
            qsten_inv_dim_sizes=qsten_inv_dim_sizes,
            qsten_inv_rank=qsten_inv_rank,
            cufrc_inv=cufrc_inv,
            cufrc_inv_dim_sizes=cufrc_inv_dim_sizes,
            cufrc_inv_rank=cufrc_inv_rank,
            fer_inv=fer_inv,
            fer_inv_dim_sizes=fer_inv_dim_sizes,
            fer_inv_rank=fer_inv_rank,
            fdr_inv=fdr_inv,
            fdr_inv_dim_sizes=fdr_inv_dim_sizes,
            fdr_inv_rank=fdr_inv_rank,
            ndrop_inv=ndrop_inv,
            ndrop_inv_dim_sizes=ndrop_inv_dim_sizes,
            ndrop_inv_rank=ndrop_inv_rank,
            nice_inv=nice_inv,
            nice_inv_dim_sizes=nice_inv_dim_sizes,
            nice_inv_rank=nice_inv_rank,
            qldet_inv=qldet_inv,
            qldet_inv_dim_sizes=qldet_inv_dim_sizes,
            qldet_inv_rank=qldet_inv_rank,
            qlsub_inv=qlsub_inv,
            qlsub_inv_dim_sizes=qlsub_inv_dim_sizes,
            qlsub_inv_rank=qlsub_inv_rank,
            qidet_inv=qidet_inv,
            qidet_inv_dim_sizes=qidet_inv_dim_sizes,
            qidet_inv_rank=qidet_inv_rank,
            qisub_inv=qisub_inv,
            qisub_inv_dim_sizes=qisub_inv_dim_sizes,
            qisub_inv_rank=qisub_inv_rank,
            tpert_out=tpert_out,
            tpert_out_dim_sizes=tpert_out_dim_sizes,
            tpert_out_rank=tpert_out_rank,
            qpert_out=qpert_out,
            qpert_out_dim_sizes=qpert_out_dim_sizes,
            qpert_out_rank=qpert_out_rank,
        )
    except Exception as err:
        print("Error in Python:")
        print(traceback.format_exc())
        return -1
    return 0


""".format(
    TMPFILEBASE
)

with open("moist.h") as f:
    data = "".join([line for line in f if not line.startswith("#")])
    data = data.replace("CFFI_DLLEXPORT", "")
    ffi.embedding_api(data)

ffi.set_source(
    TMPFILEBASE,
    '#include "moist.h"',
    library_dirs=[get_config_var("LIBDIR")],
)

ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".*", verbose=True)
