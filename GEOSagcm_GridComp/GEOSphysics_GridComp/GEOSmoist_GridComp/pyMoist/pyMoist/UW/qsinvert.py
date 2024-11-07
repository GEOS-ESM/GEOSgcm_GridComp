import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, log
import pyMoist.pyMoist_constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatFieldIJ,
    Float,
)
from ndsl import StencilFactory, QuantityFactory
from pyMoist.saturation.qsat import QSat, QSat_Float, FloatField_Extra_Dim
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.UW.compute_uwshcu import ice_fraction

@gtscript.function
def qsinvert(
    qt: Float,
    thl: Float, 
    ps_in: Float,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim
    ):

    '''
    Function calculating saturation pressure ps (or pLCL) from qt and 
    thl ( liquid potential temperature,  NOT liquid virtual potential 
    temperature) by inverting Bolton formula. I should check later if 
    current use of 'leff' instead of 'xlv' here is reasonable or not.
    '''

    psmin: f64  = f64(10000.0) # Default saturation pressure [Pa] if iteration does not converge
    dpsmax: f64 = f64(1.0)     # Tolerance [Pa] for convergence of iteration
    p00 = 1E5
    rovcp = constants.r/constants.cp

    # Calculate best initial guess of pLCL

    Ti: f64       =  thl*(ps_in/p00)**rovcp
    Tgeos: f32    = Ti
    Pgeos: f32    = f32(ps_in)
    qs, dqsdT     = QSat_Float(ese, esx, Tgeos, Pgeos/100.0)
    es: f64       = ps_in * qs  / (constants.ep2 + (f64(1.0)-constants.ep2)*f64(qs) )
    rhi: f64      = qt/f64(qs)


    if rhi <= f64(0.01):
        #print('Source air is too dry and pLCL is set to psmin in uwshcu.F90')
        qsinvert = psmin

    else:

        TLCL: f64     =  f64(55.0) + f64(1.0)/(f64(1.0)/(Ti-f64(55.0))-log(rhi)/f64(2840.0)) # Bolton's formula. MWR.1980.Eq.(22)
        PiLCL: f64    =  TLCL/thl
        ps            =  p00*(PiLCL)**(f64(1.0)/rovcp)

        iteration = 0
        while iteration < 10:
            Pis: f64        =  (ps/p00)**rovcp   # Exner function
            Ts: f64         =  thl*Pis
            Tgeos           = Ts
            Pgeos           = ps
            qs, dqsdT       = QSat_Float(ese, esx, Tgeos, Pgeos/100.0, DQSAT_trigger=True)
            gam: f64        = (constants.xlv/constants.cp)*f64(dqsdT)
            err: f64        =  qt - qs
            nu: f64         =  ice_fraction(f32(Ts),0.0,0.0)       
            leff: f64       =  (f64(1.0) - nu)*constants.xlv + nu*constants.xls                  
            dlnqsdT: f64    =  gam*(constants.cp/leff)/qs
            dTdPis: f64     =  thl
            dPisdps: f64    =  rovcp*Pis/ps 
            dlnqsdps: f64   = f64(-1.0)/(ps - (1. - constants.ep2)*es)
            derrdps: f64    = -qs*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
            #if derrdps = 0:
                #print("QSINVERT: derrdps=0 !!!")
            dps:f64         = -err/derrdps
            ps: f64         =  ps + dps
            
            if ps < f64(0.0):
                qsinvert: f32 = psmin
                iteration=10
                    
            elif abs(dps) <= dpsmax:
                qsinvert: f32 = ps
                iteration=10

            else:
                qsinvert: f32 = psmin
                
            iteration += 1

    return f32(qsinvert)

def Qs_Invert(
    qtsrc: FloatFieldIJ, 
    thlsrc: FloatFieldIJ, 
    pifc0:FloatFieldIJ, 
    plcl: FloatFieldIJ, 
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim
):

    with computation(FORWARD), interval(...):
        plcl = qsinvert(qtsrc,thlsrc,pifc0,ese,esx) 


class QsInvert:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._Qs_Invert = self.stencil_factory.from_dims_halo(
            func=Qs_Invert,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        qtsrc: FloatFieldIJ,
        thlsrc: FloatFieldIJ,
        pifc0: FloatFieldIJ,
        plcl: FloatFieldIJ,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):
        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        self._Qs_Invert(
            qtsrc=qtsrc,
            thlsrc=thlsrc,
            pifc0=pifc0,
            plcl=plcl,
            ese=self.qsat.ese,
            esx=self.qsat.esx
        )

  
    

