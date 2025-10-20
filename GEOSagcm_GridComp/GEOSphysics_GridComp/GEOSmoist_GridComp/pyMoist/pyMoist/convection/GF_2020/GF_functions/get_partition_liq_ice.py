from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
    IntFieldIJ,
    IntField,
)
import gt4py.cartesian.gtscript as gtscript
from pyMoist.UW.uwshcu_functions import ice_fraction
import pyMoist.constants as constants

T_0 = 273.16 # K

@gtscript.function
def fract_liq_f(
    temp2: Float,
    cnvfrc: Float,
    srftype: Float,
    FRAC_MODIS: Int,
):

   if FRAC_MODIS == 1:
       fract_liq_f = 1.0 - ice_fraction(temp2,cnvfrc,srftype)
  
   return fract_liq_f

def get_partition_liq_ice(
     # In
    MELT_GLAC: IntField,
    cumulus: IntField,
    ierr: IntField,
    cnvfrc: FloatField,
    srftype: FloatField,
    tn: FloatField,
    z1: FloatField,
    zo_cup: FloatField,
    po_cup: FloatField,
    FRAC_MODIS: IntField,
    norm: FloatFieldIJ,
    # Out
    melting_layer: FloatField,
    p_liq_ice: FloatField,
):
    from __externals__ import k_end

    with computation(FORWARD), interval(...):
       p_liq_ice = 1.0
       melting_layer = 0.0
       delT=3.0
       norm=0.
       ktf = k_end-1

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if K <= ktf:
                if ierr == 0:
                    p_liq_ice = fract_liq_f(tn,cnvfrc,srftype,FRAC_MODIS)
        
    with computation(PARALLEL), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if K <= ktf:
                if ierr == 0:
                    if tn <= (T_0-delT): 
                        melting_layer = 0.

                    elif tn < (T_0+delT) and tn > (T_0-delT):
                        melting_layer =  ((tn-(T_0-delT))/(2.*delT))**2

                    else:
                        melting_layer = 1.
                
                    melting_layer = melting_layer*(1.-melting_layer)

    with computation(FORWARD), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if K <= ktf-1:
                if ierr == 0:
                    dp = 100.*(po_cup-po_cup[0,0,1])
                    norm = norm + melting_layer*dp/constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if ierr == 0:
                melting_layer=melting_layer/(norm+1.e-6)*(100*(po_cup.at(K=0)-po_cup.at(K=ktf))/constants.MAPL_GRAV)
        



class GetPartitionLiqIce:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

     
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._get_partition_liq_ice = self.stencil_factory.from_dims_halo(
            func=get_partition_liq_ice,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.norm = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )

    def __call__(
        self,
        # In
        MELT_GLAC: IntField,
        cumulus: IntField,
        ierr: IntField,
        cnvfrc: FloatField,
        srftype: FloatField,
        tn: FloatField,
        z1: FloatField,
        zo_cup: FloatField,
        po_cup: FloatField,
        FRAC_MODIS: IntField,
        # Out
        melting_layer: FloatField,
        p_liq_ice: FloatField,
    ):
        if FRAC_MODIS.view[:].all() != Int(1):
            raise NotImplementedError(f"Warning: This code has not been ported!! Expecting FRAC_MODIS = 1, got FRAC_MODIS != 1")

        self._get_partition_liq_ice(
            # In
            MELT_GLAC=MELT_GLAC,
            cumulus=cumulus,
            ierr=ierr,
            cnvfrc=cnvfrc,
            srftype=srftype,
            tn=tn,
            z1=z1,
            zo_cup=zo_cup,
            po_cup=po_cup,
            FRAC_MODIS=FRAC_MODIS,
            norm=self.norm,
            # Out
            melting_layer=melting_layer,
            p_liq_ice=p_liq_ice,
        )

           
