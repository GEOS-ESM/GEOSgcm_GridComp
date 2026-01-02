import xarray as xr
import numpy as np
import os
import re

# Path to netcdf files
moist_path = "/Users/kfandric/netcdf/"

output_dir = moist_path
os.makedirs(output_dir, exist_ok=True)

# Savepoint you are testing
test_name = "PrepareInputs"

# Merge ComputeUwshcu-In.nc with infiles needed for translate test
merge_in_files = "ncks -A "+moist_path+"ComputeUwshcu-In.nc "+moist_path+test_name+"-In.nc"
os.system(merge_in_files)

# Define which dimensions to keep for each variable
keep_dim_maps = {
    "PrepareInputs-Out.nc": {
        "qt0": [0, 1],
        "s0": [0, 1],
        "ssqt0": [0,1],
        "ssthl0": [0,1],
        "sstr0": [0,1,2],
        "sstr0_o": [0,1,2],
        "ssu0": [0,1],
        "ssv0":[0,1],
        "t0": [0,1],
        "thl0": [0,1],
        "thv0bot": [0,1],
        "thv0top": [0,1],
        "thvl0": [0,1],
        "thvl0bot": [0,1],
        "thvl0top": [0,1],
        "tr0_o":[0,1,2],
        "trflx": [0,1,2],
        "trten": [0,1,2],
        "tru": [0,1,2],
        "tru_emf": [0,1,2],
        "tke": [0,1],
    },
    # "FindPbl-In.nc": {
    #     "condensation": [0],
    #     "cush": [0],
    #     "kpbl_in": [0],
    #     "pifc0": [0,1],
    #     "qt0": [0,1],
    #     "thvl0": [0,1],
    #     "thvl0bot": [0,1],
    #     "thvl0top": [0,1],
    #     "tke": [0,1],
    #     "u0": [0,1],
    #     "v0": [0,1],
    #     "zmid0": [0,1],
    # },
    # "FindPbl-Out.nc": {
    #     "kinv": [0,1],
    #     "qtavg": [0,1],
    #     "thvlavg": [0,1],
    #     "thvlmin": [0,1],
    #     "tkeavg": [0,1],
    #     "uavg": [0,1],
    #     "vavg": [0,1],
    #     "tscaleh": [0],
    # },
    # "FindKlcl-In.nc": {
    #     "condensation": [0],
    #     "pifc0": [0,1],
    #     "pmid0": [0,1],
    #     "qt0": [0,1],
    #     "ssu0": [0,1],
    #     "ssv0": [0,1],
    #     "thvlmin": [0,1],
    #     "u0": [0,1],
    #     "v0": [0,1],
    #     "t0": [0,1],
    #     "qv0": [0,1],
    #     "qtavg": [0,1],
    #     "uavg": [0,1],
    #     "vavg":[0,1],
    #     "kinv": [0,1],
    #     "thl0":[0,1],
    #     "ssthl0":[0,1],
    #     "ssqt0":[0,1],
    #     "tr0_FindKlcl": [0,1,2],
    # },
    # "FindKlcl-Out.nc": {
    #     "qtsrc": [0,1],
    #     "thlsrc": [0,1],
    #     "thvlsrc": [0,1],
    #     "usrc": [0,1],
    #     "vsrc": [0,1],
    #     "trsrc": [0,2],
    #     "klcl": [0,1],
    #     "plcl": [0,1],
    #     "qt0lcl": [0,1],
    #     "thl0lcl": [0,1],
    #     "thv0lcl": [0,1],
    # },
    # "ComputeCinCinlcl-In.nc": {
    #     "condensation": [0],
    #     "kinv":[0,1],
    #     "klcl":[0,1],
    #     "pifc0":[0,1],
    #     "plcl":[0,1],
    #     "qtsrc":[0,1],
    #     "thlsrc":[0,1],
    #     "thv0bot":[0,1],
    #     "thv0lcl":[0,1],
    #     "thv0top":[0,1],
    #     "thvlmin":[0,1],
    #     "thvlsrc":[0,1],
    #     "tkeavg":[0,1],
    #     "trsrc":[0,2],
    #     "usrc":[0,1],
    #     "vsrc":[0,1],
    #     "rkfre": [0],
    # },
    # "ComputeCinCinlcl-Out.nc": {
    #     "cin_i": [0],
    #     "cinlcl_i": [0],
    #     "ke": [0],
    #     "kinv_o": [0,1],
    #     "klfc_o": [0,1],
    #     "plcl_o": [0,1],
    #     "plfc_o": [0,1],
    #     "qtsrc_o": [0,1],
    #     "thlsrc_o": [0,1],
    #     "thvlmin_o": [0,1],
    #     "thv0lcl_o": [0,1],
    #     "tkeavg_o": [0,1],
    #     "trsrc_o": [0,2],
    #     "usrc_o": [0,1],
    #     "vsrc_o": [0,1],
    # },
    # "ComputeDelCIN-In.nc": {
    #     "condensation":[0],
    #     "cin_i": [0],
    #     "cinlcl_i": [0],
    #     "cin":[0],
    #     "cinlcl":[0],
    # },
    # "ComputeDelCIN-Out.nc": {
    #     "del_CIN":[0],
    # },
    # "AverageInitialFinalCIN1-In.nc": {
    #     "condensation":[0],
    #     "cin_i": [0],
    #     "cinlcl_i": [0],
    #     "cin":[0],
    #     "cinlcl":[0],
    #     "del_CIN":[0],
    #     "ke":[0],
    # },
    # "AverageInitialFinalCIN1-Out.nc": {
    #     "cin":[0],
    #     "cinlcl":[0],
    # },
    # "AverageInitialFinalCIN3-In.nc": {
    #     "condensation":[0],
    #     "qtflx_s": [0,1],
    #     "slflx_s": [0,1],
    #     "uflx_s":[0,1],
    #     "vflx_s":[0,1],
    #     "umf_s":[0,1],
    #     "zifc0":[0,1],
    #     "del_CIN":[0],
    #     "kinv":[0,1],
    # },
    # "AverageInitialFinalCIN3-Out.nc": {
    #     "qtflx_out":[0,1],
    #     "slflx_out":[0,1],
    #     "uflx_out":[0,1],
    #     "vflx_out":[0,1],
    #     "umf_out":[0,1],
    # },
    # "DefinePrelCbmf-In.nc": {
    #     "condensation": [0],
    #     "cin": [0],
    #     "cinlcl": [0],
    #     "dp0": [0,1],
    #     "exnifc0": [0,1],
    #     "kinv": [0,1],
    #     "klcl": [0,1],
    #     "pifc0": [0,1],
    #     "plcl": [0,1],
    #     "thv0bot": [0,1],
    #     "thv0lcl": [0,1],
    #     "thv0top": [0,1],
    #     "tkeavg": [0,1],
    #     "rkfre": [0],
    # },
    # "DefinePrelCbmf-Out.nc": {
    #     "cbmf": [0,1],
    #     "krel": [0,1],
    #     "prel": [0,1],
    #     "thv0rel": [0,1],
    #     "ufrcinv": [0,1],
    #     "wcrit": [0],
    #     "winv": [0,1],
    # },
    # "DefineUpdraftProperties-In.nc": {
    #     "cbmf": [0,1],
    #     "cinlcl": [0],
    #     "condensation": [0],
    #     "kinv": [0,1],
    #     "krel": [0,1],
    #     "prel": [0,1],
    #     "qtsrc": [0,1],
    #     "rho0inv": [0,1],
    #     "thlsrc": [0,1],
    #     "winv": [0,1],
    #     "thvu":[0,1],
    # },
    # "DefineUpdraftProperties-Out.nc": {
    #     "condensation_out": [0],
    #     "thvu_out": [0,1],
    #     "ufrc": [0,1],
    #     "ufrclcl": [0,1],
    # },
    # "DefineEnvProperties-In.nc": {
    #     "condensation": [0],
    #     "kinv": [0,1],
    #     "krel": [0,1],
    #     "pifc0": [0,1],
    #     "pmid0": [0,1],
    #     "prel": [0,1],
    #     "qt0": [0,1],
    #     "ssqt0": [0,1],
    #     "ssthl0": [0,1],
    #     "sstr0": [0,1,2],
    #     "ssu0": [0,1],
    #     "ssv0": [0,1],
    #     "thl0": [0,1],
    #     "trsrc": [0,2],
    #     "u0": [0,1],
    #     "v0": [0,1],
    #     "usrc": [0,1],
    #     "vsrc":[0,1],
    #     "thv0rel":[0,1],
    #     "tr0_DefineEnvProperties":[0,1,2],
    #     "uu":[0,1],
    #     "vu":[0,1],
    # },
    # "DefineEnvProperties-Out.nc": {
    #     "dpe": [0],
    #     "exne": [0],
    #     "pe": [0],
    #     "qsat_pe": [0,1],
    #     "qte": [0],
    #     "thle":[0],
    #     "thvebot": [0],
    #     "tre":[0,2],
    #     "tru":[0,1,2],
    #     "ue":[0],
    #     "uplus":[0],
    #     "uu":[0,1],
    #     "ve":[0],
    #     "vplus":[0],
    #     "vu":[0,1],
    # },
    # "BuoyancySorting-In.nc": {
    #     "condensation": [0],
    #     "dp0":[0,1],
    #     "exnifc0":[0,1],
    #     "exnmid0":[0,1],
    #     "krel":[0,1],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "qt0":[0,1],
    #     "qtu":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "thl0":[0,1],
    #     "thlu":[0,1],
    #     "thv0bot":[0,1],
    #     "thv0rel":[0,1],
    #     "thv0top":[0,1],
    #     "tscaleh":[0],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "umf":[0,1],
    #     "uu":[0,1],
    #     "vu":[0,1],
    #     "thvu":[0,1],
    #     "wlcl":[0,1],
    #     "wu":[0,1],
    #     "zifc0":[0,1],
    #     "tr0_BuoySort":[0,1,2],
    #     "prel":[0,1],
    #     "zmid0":[0,1],
    #     "qsat_pe":[0,1],
    # },
    # "BuoyancySorting-Out.nc": {
    #     "testvar3D_1":[0,1],
    #     "testvar3D_2":[0,1],
    #     "testvar3D_3":[0,1],
    #     "testvar3D_4":[0,1],
    #     "testvar3D_5":[0,1],
    # },
    # "CalcPpen-In.nc": {
    #     "bogbot":[0],
    #     "bogtop":[0],
    #     "condensation": [0],
    #     "dp0":[0,1],
    #     "drage":[0],
    #     "kpen":[0],
    #     "pifc0": [0,1],
    #     "rhomid0j":[0],
    #     "wtwb":[0],
    #     "wu":[0,1],
    # },
    # "CalcPpen-Out.nc": {
    #     "ppen":[0],
    # },
    # "RecalcCondensate-In.nc": {
    #     "condensation": [0],
    #     "exnifc0":[0,1],
    #     "fer":[0,1],
    #     "kbup":[0],
    #     "kpen": [0,1],
    #     "pifc0":[0,1],
    #     "ppen":[0],
    #     "qt0":[0,1],
    #     "qtu":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "thl0":[0,1],
    #     "thlu":[0,1],
    #     "thv0bot":[0,1],
    #     "thv0top":[0,1],
    #     "zifc0":[0,1],
    #     "krel": [0,1],
    # },
    # "RecalcCondensate-Out.nc": {
    #     "cldhgt": [0],
    #     "diten":[0,1],
    #     "dwten":[0,1],
    #     "emf":[0,1],
    #     "fdr": [0,1],
    #     "qtu_top":[0],
    #     "thlu_top":[0],
    #     "ufrc":[0,1],
    #     "umf":[0,1],
    #     "xco":[0,1],
    # },
    # "CalcEntrainmentMassFlux-In.nc": {
    #     "condensation":[0],
    #     "dp0":[0,1],
    #     "exnifc0":[0,1],
    #     "kbup":[0,1],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "ppen":[0],
    #     "qt0":[0,1],
    #     "qtu":[0,1],
    #     "rei":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "thl0":[0,1],
    #     "thlu":[0,1],
    #     "thv0bot":[0,1],
    #     "thv0top":[0,1],
    #     "tru":[0,1,2],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "uu":[0,1],
    #     "vu":[0,1],
    #     "umf":[0,1],
    #     "kpen":[0,1],
    #     "tr0_CalcEntrain":[0,1,2],
    # },
    # "CalcEntrainmentMassFlux-Out.nc": {
    #   "qtu_emf": [0,1],
    #   "thlu_emf": [0,1],
    #   "tru_emf":[0,1,2],
    #   "uu_emf":[0,1],
    #   "vu_emf":[0,1],
    #   "emf":[0,1],
    # },
    # "CalcPblFluxes-In.nc": {
    #     "cbmf":[0,1],
    #     "condensation":[0],
    #     "kinv":[0,1],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "qt0":[0,1],
    #     "qtsrc":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "thl0":[0,1],
    #     "thlsrc":[0,1],
    #     "tr0_PblFlux":[0,1,2],
    #     "trsrc":[0,2],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "vsrc":[0,1],
    #     "usrc":[0,1],
    #     "exnifc0":[0,1],
    # },
    # "CalcPblFluxes-Out.nc": {
    #     "qtflx": [0,1],
    #     "slflx":[0,1],
    #     "trflx":[0,1,2],
    #     "uflx":[0,1],
    #     "vflx":[0,1],
    #     "xflx":[0,1],
    # },
    # "BuoyancySortingFluxes-In.nc": {
    #     "condensation":[0],
    #     "exnifc0":[0,1],
    #     "kbup":[0,1],
    #     "krel":[0,1],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "qt0":[0,1],
    #     "qtflx":[0,1],
    #     "qtu":[0,1],
    #     "slflx":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "thl0":[0,1],
    #     "thlu":[0,1],
    #     "tr0_BuoySortFlux":[0,1,2],
    #     "trflx":[0,1,2],
    #     "tru":[0,1,2],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "vflx":[0,1],
    #     "vu":[0,1],
    #     "uu":[0,1],
    #     "uflx":[0,1],
    #     "umf":[0,1],
    # },
    # "BuoyancySortingFluxes-Out.nc": {
    #     "qtflx":[0,1],
    #     "slflx":[0,1],
    #     "trflx":[0,1,2],
    #     "vflx":[0,1],
    #     "uflx":[0,1],
    # },
    # "PenetrativeEntrainmentFluxes-In.nc": {
    #     "condensation":[0],
    #     "emf":[0,1],
    #     "kbup":[0,1],
    #     "kpen":[0,1],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "qt0":[0,1],
    #     "qtflx":[0,1],
    #     "qtu_emf":[0,1],
    #     "slflx":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "thl0":[0,1],
    #     "thlu_emf":[0,1],
    #     "tr0_PenEntrainFlux":[0,1,2],
    #     "trflx":[0,1,2],
    #     "tru_emf":[0,1,2],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "vflx":[0,1],
    #     "vu_emf":[0,1],
    #     "uu_emf":[0,1],
    #     "uflx":[0,1],
    #     "umf":[0,1],
    #     "ql0":[0,1],
    #     "qi0":[0,1],
    #     "kinv":[0,1],
    #     "krel":[0,1],
    #     "cbmf":[0,1],
    #     "exnifc0":[0,1],
    # },
    # "PenetrativeEntrainmentFluxes-Out.nc": {
    #     "qtflx":[0,1],
    #     "slflx":[0,1],
    #     "trflx":[0,1,2],
    #     "vflx":[0,1],
    #     "uflx":[0,1],
    #     "qlten_sink":[0,1],
    #     "qiten_sink":[0,1],
    #     "uemf":[0,1],
    # },
    # "MomentumTendency-In.nc": {
    #     "dp0":[0,1],
    #     "condensation":[0],
    #     "kpen":[0,1],
    #     "uflx":[0,1],
    #     "vflx":[0,1],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "uf":[0,1],
    #     "vf":[0,1],
    # },
    # "MomentumTendency-Out.nc": {
    #     "uf":[0,1],
    #     "vf":[0,1],
    #     "uten":[0,1],
    #     "vten":[0,1],
    # },
    # "ThermodynamicTendencies-In.nc": {
    #     "condensation":[0],
    #     "diten":[0,1],
    #     "dwten":[0,1],
    #     "dp0":[0,1],
    #     "kpen":[0,1],
    #     "qtflx":[0,1],
    #     "slflx":[0,1],
    #     "u0":[0,1],
    #     "v0":[0,1],
    #     "uf":[0,1],
    #     "vf":[0,1],
    #     "uflx":[0,1],
    #     "vflx":[0,1],
    #     "umf":[0,1],
    #     "pifc0":[0,1],
    #     "ppen":[0],
    #     "prel":[0,1],
    #     "qtu":[0,1],
    #     "thlu":[0,1],
    #     "thlu_top":[0],
    #     "qtu_top":[0],
    #     "emf":[0,1],
    #     "ql0":[0,1],
    #     "qi0":[0,1],
    #     "pmid0":[0,1],
    #     "kbup":[0,1],
    #     "krel":[0,1],
    #     "thlu_emf":[0,1],
    #     "qtu_emf":[0,1],
    #     "qlten_sink":[0,1],
    #     "qiten_sink":[0,1],
    #     "fdr":[0,1],
    # },
    # "ThermodynamicTendencies-Out.nc": {
    #     "slten":[0,1],
    #     "qc":[0,1],
    #     "sten":[0,1],
    #     "qiten":[0,1],
    #     "qlten":[0,1],
    #     "qvten":[0,1],
    # },
    # "PreventNegativeCondensate-In.nc": {
    #     "condensation":[0],
    #     "dp0":[0,1],
    #     "qi0":[0,1],
    #     "qiten":[0,1],
    #     "ql0":[0,1],
    #     "qlten":[0,1],
    #     "qtten":[0,1],
    #     "qv0":[0,1],
    #     "qvten":[0,1],
    #     "s0":[0,1],
    #     "slten":[0,1],
    #     "sten":[0,1],
    # },
    # "PreventNegativeCondensate-Out.nc": {
    #     "qlten":[0,1],
    #     "qiten":[0,1],
    #     "qvten":[0,1],
    # },
    # "TracerTendencies-In.nc": {
    #     "condensation":[0],
    #     "dp0":[0,1],
    #     "tr0_TracerTendencies":[0,1,2],
    #     "trflx":[0,1,2],
    # },
    # "TracerTendencies-Out.nc": {
    #     "trten":[0,1,2],
    # },
    # "ComputeDiagnosticOutputs-In.nc": {
    #     "condensation":[0],
    #     "prel":[0,1],
    #     "qtu":[0,1],
    #     "thlu":[0,1],
    #     "krel":[0,1],
    # },
    # "ComputeDiagnosticOutputs-Out.nc": {
    #     "qcubelow":[0],
    #     "qiubelow":[0],
    #     "qlubelow":[0],
    #     "rcwp":[0],
    #     "riwp":[0],
    #     "rlwp":[0],
    # },
    #  "CalcCumulusCondensate-In.nc": {
    #     "condensation":[0],
    #     "kpen":[0,1],
    #     "krel":[0,1],
    #     "pifc0":[0,1],
    #     "ppen":[0],
    #     "prel":[0,1],
    #     "qcubelow":[0],
    #     "qiubelow":[0],
    #     "qlubelow":[0],
    #     "qtu":[0,1],
    #     "qtu_top":[0],
    #     "rcwp":[0],
    #     "riwp":[0],
    #     "rlwp":[0],
    #     "thlu":[0,1],
    #     "thlu_top":[0],
    #     "ufrc":[0,1],
    #     "ufrclcl":[0,1],
    # },
    # "CalcCumulusCondensate-Out.nc": {
    #     "qcubelow":[0],
    #     "qiubelow":[0],
    #     "qlubelow":[0],
    #     "rcwp":[0],
    #     "riwp":[0],
    #     "rlwp":[0],
    # },
    # "AdjustImplicitCINInputs1-In.nc": {
    #     "condensation":[0],
    #     "qi0":[0,1],
    #     "qiten":[0,1],
    #     "ql0":[0,1],
    #     "qlten":[0,1],
    #     "qv0":[0,1],
    #     "qvten":[0,1],
    #     "s0":[0,1],
    #     "sten":[0,1],
    #     "t0":[0,1],
    #     "tr0_AdjustCIN":[0,1,2],
    #     "trten":[0,1,2],
    #     "u0":[0,1],
    #     "uten":[0,1],
    #     "v0":[0,1],
    #     "vten":[0,1],
    # },
    # "AdjustImplicitCINInputs1-Out.nc": {
    #    "qi0_s":[0,1],
    #    "qiten_s":[0,1],
    #    "ql0_s":[0,1],
    #    "qlten_s":[0,1],
    #    "qv0_s":[0,1],
    #    "qvten_s":[0,1],
    #    "s0_s":[0,1],
    #    "sten_s":[0,1],
    #    "t0_s":[0,1],
    #    "tr0_s":[0,1,2],
    #    "u0_s":[0,1],
    #    "uten_s":[0,1],
    #    "v0_s":[0,1],
    #    "vten_s":[0,1],
    # },
    # "AdjustImplicitCINInputs2-In.nc": {
    #     "condensation":[0],
    #     "qtflx":[0,1],
    #     "slflx":[0,1],
    #     "uflx":[0,1],
    #     "ufrc":[0,1],
    #     "umf":[0,1],
    #     "vflx":[0,1],
    # },
    # "AdjustImplicitCINInputs2-Out.nc": {
    #     "qtflx_s":[0,1],
    #     "slflx_s":[0,1],
    #     "uflx_s":[0,1],
    #     "ufrc_s":[0,1],
    #     "umf_s":[0,1],
    #     "vflx_s":[0,1],
    # },
    # "RecalcEnvVariables-In.nc": {
    #     "condensation":[0],
    #     "pifc0":[0,1],
    #     "pmid0":[0,1],
    #     "exnmid0":[0,1],
    #     "qi0_s":[0,1],
    #     "ql0_s":[0,1],
    #     "qv0_s":[0,1],
    #     "s0_s":[0,1],
    #     "t0_s":[0,1],
    #     "tr0_RecalcEnv":[0,1,2],
    #     "u0":[0,1],
    #     "v0":[0,1],
        
    # },
    # "RecalcEnvVariables-Out.nc": {
    #     "qi0":[0,1],
    #     "ql0":[0,1],
    #     "ql0":[0,1],
    #     "qt0":[0,1],
    #     "qv0":[0,1],
    #     "s0":[0,1],
    #     "ssqt0":[0,1],
    #     "ssthl0":[0,1],
    #     "sstr0":[0,1,2],
    #     "ssu0":[0,1],
    #     "ssv0":[0,1],
    #     "t0":[0,1],
    #     "thl0":[0,1],
    #     "thv0bot":[0,1],
    #     "thv0top":[0,1],
    #     "thvl0top":[0,1],
    #     "thvl0bot":[0,1],
    # },
    # "UpdateOutputVars1-In.nc": {
    #     "condensation":[0],
    #     "cufrc":[0,1],
    #     "cush":[0],
    #     "dcm":[0,1],
    #     "kinv":[0,1],
    #     "qiten":[0,1],
    #     "qlten":[0,1],
    #     "qrten":[0,1],
    #     "qsten":[0,1],
    #     "qvten":[0,1],
    #     "sten":[0,1],
    #     "umf":[0,1],
    #     "uten":[0,1],
    #     "vten":[0,1],
    #     "zifc0":[0,1],
        
    # },
    # "UpdateOutputVars1-Out.nc": {
    #     "cufrc_out":[0,1],
    #     "cush_inout":[0,1],
    #     "dcm_out":[0,1],
    #     "qiten_out":[0,1],
    #     "qlten_out":[0,1],
    #     "qrten_out":[0,1],
    #     "qsten_out":[0,1],
    #     "qvten_out":[0,1],
    #     "sten_out":[0,1],
    #     "umf_out":[0,1],
    #     "uten_out":[0,1],
    #     "vten_out":[0,1],
    # },
    # "UpdateOutputVars2-In.nc": {
    #     "condensation":[0],
    #     "fdr":[0,1],
    #     "fer":[0,1],
    #     "kpen":[0,1],
    #     "qiten_det":[0,1],
    #     "qiten_sink":[0,1],
    #     "qlten_det":[0,1],
    #     "qlten_sink":[0,1],
    #     "qtflx":[0,1],
    #     "slflx":[0,1],
    #     "tr0_inout":[0,1,2],
    #     "trten":[0,1,2],
    #     "uflx":[0,1],
    #     "vflx":[0,1],
    # },
    # "UpdateOutputVars2-Out.nc": {
    #     "fdr_out":[0,1],
    #     "fer_out":[0,1],
    #     "qidet_out":[0,1],
    #     "qisub_out":[0,1],
    #     "qldet_out":[0,1],
    #     "qlsub_out":[0,1],
    #     "qtflx_out":[0,1],
    #     "slflx_out":[0,1],
    #     "tr0_inout":[0,1,2],
    #     "uflx_out":[0,1],
    #     "vflx_out":[0,1],
    #     "ndrop_out":[0,1],
    #     "nice_out":[0,1],
    # },
    # "ComputeUwshcuInvertAfter-In.nc": {
    #     "cnvtr":[0],
    #     "cufrc_out":[0,1],
    #     "cush_inout":[0],
    #     "dcm_out":[0,1],
    #     "evap":[0],
    #     "fdr_out":[0,1],
    #     "fer_out":[0,1],
    #     "ndrop_out":[0,1],
    #     "nice_out":[0,1],
    #     "qidet_out":[0,1],
    #     "qisub_out":[0,1],
    #     "qiten_out":[0,1],
    #     "qlten_out":[0,1],
    #     "qldet_out":[0,1],
    #     "qlsub_out":[0,1],
    #     "qtflx_out":[0,1],
    #     "slflx_out":[0,1],
    #     "tr0_inout":[0,1,2],
    #     "uflx_out":[0,1],
    #     "vflx_out":[0,1],
    #     "qpert_out":[0],
    #     "qrten_out":[0,1],
    #     "qsten_out":[0,1],
    #     "qvten_out":[0,1],
    #     "shfx":[0],
    #     "sten_out":[0,1],
    #     "tpert_out":[0],
    #     "uten_out":[0,1],
    #     "vten_out":[0,1],
    #     "umf_out":[0,1],
    # },
}

# Make a list of all interface variables
interface_vars = ["trflx", "tru", "tru_emf","thlu","qtu","pifc0", "zifc0","exnifc0","thvu","ufrc","thvu_out","uu","vu",\
    "wu","umf","emf","qtu_emf","uu_emf","vu_emf","thlu_emf","xflx","vflx","uflx","slflx","qtflx","trflx","uemf",\
    "qtflx_s","slflx_s","uflx_s","vflx_s","umf_s","ufrc_s","umf_out","qtflx_out","slflx_out","uflx_out","vflx_out"]
ntracer_2d = ["trsrc", "trsrc_o","tre"]

def determine_keep_dims(da: xr.DataArray, dim_indices: list[int]) -> list[str]:
    """Convert numeric indices (0,1,2) to actual dimension names for a DataArray."""
    dims = list(da.dims)
    result = []

    # always keep savepoint and rank
    result.append(dims[0])
    result.append(dims[1])

    # also keep whichever of the 0/1/2 dims are specified
    for i in dim_indices:
        if i < len(dims):
            result.append(dims[i + 2])
    return result


def flatten_redundant_dims(da: xr.DataArray, keep_dims: list[str]) -> xr.DataArray:
    """
    Collapse any dimensions of the DataArray not included in keep_dims.
    Preserves all dimensions listed in keep_dims.
    """
    # Identify dims to collapse
    dims_to_collapse = [d for d in da.dims if d not in keep_dims]

    # Collapse each unwanted dimension at index 0
    for d in dims_to_collapse:
        da = da.isel({d: 0})

    return da

def post_processing(da: xr.DataArray) -> xr.DataArray:
    """
    Unroll the IJ dimension (576 -> 24,24) while preserving all other dim sizes.
    Slice 72 or 73 levels, depending on size of z_dim.
    Slice ntracer vars.
    """
    # Fill nans with zeros
    da = da.fillna(0)

    if da.shape[2] != 576:
        raise ValueError(f"Expected first dimension length 576, got {da.shape[2]}")

    # Grab first 72 levels, except for interface vars
    if da.name not in interface_vars and da.name not in ntracer_2d:
        if da.ndim > 3:
            da = da.isel({da.dims[3]: slice(0, 72)})

    if da.name in ntracer_2d:
        da = da.isel({da.dims[3]: slice(0, 23)})
    
    # Slice ntracer vars
    if da.ndim == 5:
        da = da.isel({da.dims[4]: slice(0, 23)})

    new_shape = (da.shape[0], da.shape[1], 24, 24) + da.shape[3:]
    
    reshaped = da.values.reshape(new_shape, order="F")

    if reshaped.ndim == 4:
        new_dims=['savepoint','rank',f'dim_{da.name}_0', f"dim_{da.name}_1"]

    if reshaped.ndim == 5:
        new_dims=['savepoint','rank',f'dim_{da.name}_0', f"dim_{da.name}_1", f"dim_{da.name}_2"]
    
    if reshaped.ndim == 6:
        new_dims=['savepoint','rank',f'dim_{da.name}_0', f"dim_{da.name}_1", f"dim_{da.name}_2", f"dim_{da.name}_3"]

    return xr.DataArray(
        reshaped,
        dims=new_dims,
    )


# ------------------------------------------------------------------
for f in keep_dim_maps.keys():
    full_path = moist_path + f
    fname = os.path.basename(full_path)
    print(f"Processing {fname}")
    ds = xr.open_dataset(full_path)

    file_map = keep_dim_maps.get(fname, {})

    subset = xr.Dataset()

    # All other variables
    for v in ds.data_vars:
        da = ds[v]
        if v in file_map:
            keep_dims = determine_keep_dims(da, file_map[v])
            da_flat = flatten_redundant_dims(da, keep_dims)
            da = post_processing(da_flat)
        subset[v] = da

    # Save the new netcdf
    base = os.path.splitext(fname)[0]
    new_base = re.sub(r"-(In|Out)$", rf"-\1", base)
    if new_base == base:
        new_base = f"{base}"
    os.remove(output_dir+f"{base}.nc")
    out_path = os.path.join(output_dir, f"{new_base}.nc")

    subset.to_netcdf(out_path)
    print(f"  Saved {out_path}")

print("✅ Done :)")
