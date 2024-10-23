import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log
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

    psmin  = 10000.0 # Default saturation pressure [Pa] if iteration does not converge
    dpsmax = 1.0     # Tolerance [Pa] for convergence of iteration
    p00 = 1E5
    rovcp = constants.r/constants.cp

    # Calculate best initial guess of pLCL

    Ti       =  thl*(ps_in/p00)**rovcp
    Tgeos    = Ti
    Pgeos    = ps_in
    qs, _ = QSat_Float(ese, esx, Tgeos, Pgeos/100.0)
    es       = ps_in * qs  / (constants.ep2 + (1.0-constants.ep2)*qs )
    rhi      = qt/qs


    if rhi <= 0.01:
        #print('Source air is too dry and pLCL is set to psmin in uwshcu.F90')
        qsinvert = psmin

    else:

        TLCL     =  55.0 + 1.0/(1.0/(Ti-55.0)-log(rhi)/2840.0) # Bolton's formula. MWR.1980.Eq.(22)
        PiLCL    =  TLCL/thl
        ps       =  p00*(PiLCL)**(1.0/rovcp)

        iteration = 0
        while iteration < 10:
            Pis      =  (ps/p00)**rovcp   # Exner function
            Ts       =  thl*Pis
            Tgeos    = Ts
            Pgeos    = ps
            #qs, _       = QSat_Float(ese, esx, Tgeos, Pgeos/100.0)
            #Qgeos    = qs
            #dqsdT, _    = QSat_Float(ese, esx, Tgeos, Pgeos/100.0, QSAT=Qgeos)
            qs, dqsdT    = QSat_Float(ese, esx, Tgeos, Pgeos/100.0, DQSAT_trigger=True)
            gam      = (constants.xlv/constants.cp)*dqsdT
            err      =  qt - qs
            nu       =  ice_fraction(Ts,0.0,0.0)       
            leff     =  (1.0 - nu)*constants.xlv + nu*constants.xls                  
            dlnqsdT  =  gam*(constants.cp/leff)/qs
            dTdPis   =  thl
            dPisdps  =  rovcp*Pis/ps 
            dlnqsdps = -1.0/(ps - (1. - constants.ep2)*es)
            derrdps  = -qs*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
            #if derrdps = 0:
                    #print("QSINVERT: derrdps=0 !!!")
            dps      = -err/derrdps
            ps       =  ps + dps
            
            if ps < 0.0:
                qsinvert = psmin
                iteration=10
                    
            elif abs(dps) <= dpsmax:
                qsinvert = ps
                iteration=10

            else:
                qsinvert=psmin
                
            iteration += 1

    return qsinvert

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

  
    

