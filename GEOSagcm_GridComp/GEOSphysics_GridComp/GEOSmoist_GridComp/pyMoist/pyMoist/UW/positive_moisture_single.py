import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, BACKWARD
import pyMoist.constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatFieldIJ,
    FloatField,
    Int
)
from ndsl import StencilFactory, QuantityFactory

def positive_moisture_single(
    xlv: FloatFieldIJ,
    xls: FloatFieldIJ,
    mkx: Int,
    dt: FloatFieldIJ,
    qvmin: FloatFieldIJ,
    qlmin: FloatFieldIJ,
    qimin: FloatFieldIJ,
    dp: FloatField,
    qv: FloatField,
    ql: FloatField,
    qi: FloatField,
    s: FloatField,
    qvten: FloatField,
    qlten: FloatField,
    qiten: FloatField,
    sten: FloatField,
    k_idx: FloatField
    ):
    '''
    If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         
    force them to be larger than minimum value by (1) condensating water vapor      
    into liquid or ice, and (2) by transporting water vapor from the very lower     
    layer. '2._r8' is multiplied to the minimum values for safety.

    Update final state variables and tendencies associated with this correction.    
    If any condensation happens, update (s,t) too.      
                                
    Note that (qv,ql,qi,s) are final state variables after applying corresponding   
    input tendencies and corrective tendencies   
    '''
    with computation(BACKWARD), interval(...):
       dql:f64 = max(f64(0.0),f64(1.0)*qlmin-ql)
       dqi:f64 = max(f64(0.0),f64(1.0)*qimin-qi)
       qlten = qlten +  dql/dt
       qiten = qiten +  dqi/dt
       qvten = qvten - (dql+dqi)/dt
       sten  = sten  + xlv * (dql/dt) + xls * (dqi/dt)
       ql    = ql +  dql
       qi    = qi +  dqi
       qv    = qv -  dql - dqi
       s     = s  +  xlv * dql + xls * dqi
       dqv   = max(0.,1.*qvmin-qv)
       qvten = qvten + dqv/dt
       qv    = qv   + dqv

    with computation(BACKWARD), interval(1,None):
       qv[0,0,-1]    = qv[0,0,-1]    - dqv*dp/dp[0,0,-1]
       qvten[0,0,-1] = qvten[0,0,-1] - dqv*dp/dp[0,0,-1]/dt

    with computation(BACKWARD), interval(...):
       qv = max(qv,qvmin)
       ql = max(ql,qlmin)
       qi = max(qi,qimin)

    with computation(PARALLEL), interval(...):
       # Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
       # extracted from all the layers that has 'qv > 2*qvmin'. This fully
       # preserves column moisture. 
       if dqv > f64(1.0e-20):
            sum:f64 = 0.0
            if k_idx <= mkx:
                if qv > f64(2.0)*qvmin:
                    sum = sum + qv*dp
            aa:f64 = dqv*dp.at(K=1)/max(f64(1.0e-20),sum)
            if aa < f64(0.5):
                if k_idx <= mkx:
                    if qv > f64(2.0)*qvmin:
                        dum:f64  = aa*qv
                        qv       = qv - dum
                        qvten    = qvten - dum/dt
            # else: 
            # if scverbose == True:
                # print('Full positive_moisture is impossible in uwshcu')


class PositiveMoistureSingle:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._PositiveMoistureSingle = self.stencil_factory.from_dims_halo(
            func=positive_moisture_single,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._k_idx = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for k in range(0, self._k_idx.view[:].shape[2]):
            self._k_idx.view[:, :, k] = k
            
    def __call__(
        self,
        xlv: FloatFieldIJ,
        xls: FloatFieldIJ,
        mkx: Int,
        dt: FloatFieldIJ,
        qvmin: FloatFieldIJ,
        qlmin: FloatFieldIJ,
        qimin: FloatFieldIJ,
        dp: FloatField,
        qv: FloatField,
        ql: FloatField,
        qi: FloatField,
        s: FloatField,
        qvten: FloatField,
        qlten: FloatField,
        qiten: FloatField,
        sten: FloatField,
    ):

        self._PositiveMoistureSingle(
            xlv=xlv,
            xls=xls,
            mkx=mkx,
            dt=dt,
            qvmin=qvmin,
            qlmin=qlmin,
            qimin=qimin,
            dp=dp,
            qv=qv,
            ql=ql,
            qi=qi,
            s=s,
            qvten=qvten,
            qlten=qlten,
            qiten=qiten,
            sten=sten,
            k_idx=self._k_idx,
        )

  
    

