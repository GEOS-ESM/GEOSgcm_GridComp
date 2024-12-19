import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, sqrt, exp
import pyMoist.constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    Int,
    IntField,
    IntFieldIJ,
    FloatFieldIJ,
    FloatField
)
from ndsl import StencilFactory, QuantityFactory


def fluxbelowinv(
    cbmf: FloatFieldIJ,
    ps0: FloatField,
    kinv: IntFieldIJ,
    dt: FloatFieldIJ,
    xsrc: FloatFieldIJ,
    xmean: FloatFieldIJ,
    xtopin: FloatFieldIJ,
    xbotin: FloatFieldIJ,
    xflx: FloatField,
    k_idx: FloatField,
    ):
    '''
    Stencil to calculate turbulent fluxes at and below 'kinv-1' interfaces.
    Check in the main program such that input 'cbmf' should not be zero.        
    If the reconstructed inversion height does not go down below the 'kinv-1' 
    interface, then turbulent flux at 'kinv-1' interface  is simply a product 
    of 'cmbf' and 'qtsrc-xbot' where 'xbot' is the value at the top interface 
    of 'kinv-1' layer. This flux is linearly interpolated down to the surface 
    assuming turbulent fluxes at surface are zero. If reconstructed inversion 
    height goes down below the 'kinv-1' interface, subsidence warming &drying 
    measured by 'xtop-xbot', where  'xtop' is the value at the base interface 
    of 'kinv+1' layer, is added ONLY to the 'kinv-1' layer, using appropriate 
    mass weighting ( rpinv and rcbmf, or rr = rpinv / rcbmf ) between current 
    and next provisional time step. Also impose a limiter to enforce outliers 
    of thermodynamic variables in 'kinv' layer  to come back to normal values 
    at the next step.  

    '''
    with computation(FORWARD), interval(...):
        xflx[0, 0, 1] = 0.0 
    
    with computation(PARALLEL), interval(...):
        xflx = 0.0
        k_below=kinv-1
        dp = ps0.at(K=k_below) - ps0.at(K=kinv)    

        xbot = xbotin
        xtop = xtopin

        # Compute reconstructed inversion height
        xtop_ori = xtop
        xbot_ori = xbot
        rcbmf = ( cbmf * constants.MAPL_GRAV * dt ) / dp                  # Can be larger than 1 : 'OK'      

        if xbot >= xtop:
            rpeff = ( xmean - xtop ) / max(  1.e-20, xbot - xtop ) 
        else:
            rpeff = ( xmean - xtop ) / min( -1.e-20, xbot - xtop ) 

        rpeff = min( max(0.0,rpeff), 1.0 )                # As of this, 0<= rpeff <= 1   
        if rpeff == 0.0 or rpeff == 1.0:
            xbot = xmean
            xtop = xmean
            
        #Below two commented-out lines are the old code replacing the above 'if' block.   
        # if(rpeff.eq.1) xbot = xmean
        # if(rpeff.eq.0) xtop = xmean 
        
        rr       = rpeff / rcbmf
        pinv     = ps0.at(K=k_below) - rpeff * dp             # "pinv" before detraining mass
        pinv_eff = ps0.at(K=k_below) + ( rcbmf - rpeff ) * dp # Effective "pinv" after detraining mass

        # Compute turbulent fluxes.                                               
        # Below two cases exactly converges at 'kinv-1' interface when rr = 1. 
        if k_idx <= k_below:
            xflx = cbmf * ( xsrc - xbot ) * ( ps0.at(K=0) - ps0 ) / ( ps0.at(K=0) - pinv )
        if k_idx == k_below:
            if rr <= 1:
                xflx = xflx - ( 1.0 - rr ) * cbmf * ( xtop_ori - xbot_ori )
        

class FluxBelowInv:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._FluxBelowInv = self.stencil_factory.from_dims_halo(
            func=fluxbelowinv,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._k_idx = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for k in range(0, self._k_idx.view[:].shape[2]):
            self._k_idx.view[:, :, k] = k
            
    def __call__(
        self,
        cbmf: FloatFieldIJ,
        ps0: FloatField,
        kinv: IntFieldIJ,
        dt: FloatFieldIJ,
        xsrc: FloatFieldIJ,
        xmean: FloatFieldIJ,
        xtopin: FloatFieldIJ,
        xbotin: FloatFieldIJ,
        xflx: FloatField,
    ):

        self._FluxBelowInv(
            cbmf=cbmf,
            ps0=ps0,
            kinv=kinv,
            dt=dt,
            xsrc=xsrc,
            xmean=xmean,
            xtopin=xtopin,
            xbotin=xbotin,
            xflx=xflx,
            k_idx=self._k_idx,
        )

  
    

