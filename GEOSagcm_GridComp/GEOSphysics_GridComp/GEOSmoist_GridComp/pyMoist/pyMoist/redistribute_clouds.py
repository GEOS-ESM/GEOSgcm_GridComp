import gt4py
from gt4py.cartesian.gtscript import PARALLEL, computation, interval, stencil
from ndsl.dsl.typing import FloatField
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl import QuantityFactory, StencilFactory
import numpy as np
import pyMoist.radiation_coupling_constants as radconstants

def redist_clouds(CF: FloatField, QL: FloatField, QI: FloatField,
                        CLCN: FloatField, CLLS: FloatField,
                        QLCN: FloatField, QLLS: FloatField,
                        QICN: FloatField, QILS: FloatField,
                        QV: FloatField, TE: FloatField):
    
    with computation(PARALLEL), interval(...):
        # Constants from MAPL.h
        alhlbcp = radconstants.ALHLBCP
        alhsbcp = radconstants.ALHSBCP

        # Define FCN as a 3-d array
        FCN = CF

        # Fix cloud quants if too small
        if QL + QI < 1e-8:
            QV = QV + QL + QI
            TE = TE - alhlbcp * QL - alhsbcp * QI
            CF = 0.0
            QL = 0.0
            QI = 0.0
            
        if CF < 1e-5:
            QV = QV + QL + QI
            TE = TE - (alhlbcp*QL) - (alhsbcp*QI)
            CF = 0.0
            QL = 0.0
            QI = 0.0
        
        # Redistribute liquid CN/LS portions based on prior fractions
        FCN = 0.0
        if QLCN + QLLS > 0.0:
            FCN = min(max(QLCN / (QLCN + QLLS), 0.0), 1.0)
        
        # Put all new condensate into LS
        DQC = QL - (QLCN + QLLS)
        if DQC > 0.0:
            QLLS = QLLS + DQC
            DQC = 0.0
        
        # Any loss of condensate uses the FCN ratio
        QLCN = QLCN + DQC * FCN
        QLLS = QLLS + DQC * (1.0 - FCN)
        
        # Redistribute ice CN/LS portions based on prior fractions
        FCN = 0.0
        if QICN + QILS > 0.0:
            FCN = min(max(QICN / (QICN + QILS), 0.0), 1.0)
        
        # Put all new condensate into LS
        DQC = QI - (QICN + QILS)
        if DQC > 0.0:
            QILS = QILS + DQC
            DQC = 0.0
        
        # Any loss of condensate uses the FCN ratio
        QICN = QICN + DQC * FCN
        QILS = QILS + DQC * (1.0 - FCN)
        
        # Redistribute cloud-fraction CN/LS portions based on prior fractions
        FCN = 0.0
        if CLCN + CLLS > 0.0:
            FCN = min(max(CLCN / (CLCN + CLLS), 0.0), 1.0)

        # Put all new condensate into LS
        DQC = CF - (CLCN + CLLS)
        if DQC > 0.0:
            CLLS = CLLS + DQC
            DQC = 0.0
        
        # Any loss of condensate uses the FCN ratio
        CLCN = CLCN + DQC * FCN
        CLLS = CLLS + DQC * (1.0 - FCN)


class RedistributeClouds:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,

    )-> None:
      
        self._redist_clouds = stencil_factory.from_dims_halo(
            func=redist_clouds, 
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
       
    def __call__(
        self,
        RAD_CF: FloatField,
        RAD_QL: FloatField,
        RAD_QI: FloatField,
        CLCN: FloatField,
        CLLS: FloatField,
        QLCN: FloatField,
        QLLS: FloatField,
        QICN: FloatField,
        QILS: FloatField,
        RAD_QV: FloatField,
        T: FloatField,
    ):
        '''
        Perform the redistribute clouds process.

        Parameters:
        RAD_CF (FloatField): Radiation cloud fraction.
        RAD_QL (FloatField): Radiation liquid cloud mixing ratio.
        RAD_QI (FloatField): Radiation ice cloud mixing ratio.
        CLCN (FloatField): Cloud fraction (anvil).
        CLLS (FloatField): Cloud fraction (large-scale).
        QLCN (FloatField): Liquid cloud mixing ratio (anvil).
        QLLS (FloatField): Liquid cloud mixing ratio (large-scale).
        QICN (FloatField): Ice cloud mixing ratio (anvil).
        QILS (FloatField): Ice cloud mixing ratio (large-scale).
        RAD_QV (FloatField): Radiation water vapor mixing ratio.
        T (FloatField): Temperature.
        '''
        self._redist_clouds(CF = RAD_CF, 
                            QL = RAD_QL, 
                            QI = RAD_QI, 
                            CLCN = CLCN,
                            CLLS = CLLS, 
                            QLCN = QLCN,  
                            QLLS = QLLS, 
                            QICN = QICN,
                            QILS = QILS,
                            QV = RAD_QV,
                            TE = T,
                            )
