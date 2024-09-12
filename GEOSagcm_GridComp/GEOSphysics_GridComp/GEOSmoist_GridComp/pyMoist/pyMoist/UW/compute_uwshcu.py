import copy

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD
import pyMoist.pyMoist_constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, Int 
from ndsl import StencilFactory, QuantityFactory
from pyMoist.types import FloatField_NTracers
#from pyMoist.slope import Slope, calc_slope
from ndsl.quantity import Quantity
import numpy as np


'''
def slope(field: FloatField, p0: FloatField, slope: FloatField):
    with computation(PARALLEL), interval(0, 1):
        value = (field[0, 0, 1] - field[0, 0, 0]) / (p0[0, 0, 1] - p0[0, 0, 0])
        if value > 0.0:
            slope = max(0.0, value)
        else:
            slope = min(0.0, value)
    with computation(PARALLEL), interval(1, -1):
        above_value = (field[0, 0, 1] - field[0, 0, 0]) / (p0[0, 0, 1] - p0[0, 0, 0])
        below_value = (field[0, 0, 0] - field[0, 0, -1]) / (p0[0, 0, 0] - p0[0, 0, -1])
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))
    with computation(PARALLEL), interval(-1, None):
        slope = slope[0, 0, -1]

    return slope
'''

@gtscript.function
def compute_slope(field: FloatField, p0: FloatField)-> Float:
    value = (field[0, 0, 1] - field[0, 0, 0]) / (p0[0, 0, 1] - p0[0, 0, 0])
    if value > 0.0:
        slope = max(0.0, value)
    else:
        slope =  min(0.0, value)

    return slope

@gtscript.function
def compute_above_below_slope(field: FloatField, p0: FloatField) -> Float:
    above_value = (field[0, 0, 1] - field[0, 0, 0]) / (p0[0, 0, 1] - p0[0, 0, 0])
    below_value = (field[0, 0, 0] - field[0, 0, -1]) / (p0[0, 0, 0] - p0[0, 0, -1])
    if above_value > 0.0:
        slope = max(0.0, min(above_value, below_value))
    else:
        slope = min(0.0, max(above_value, below_value))

    return slope




def compute_uwshcu(
    dotransport: Float,           
    exnifc0_in: FloatField,     
    pmid0_in: FloatField,       
    zmid0_in: FloatField,       
    exnmid0_in: FloatField,     
    u0_in: FloatField,           
    v0_in: FloatField,          
    qv0_in: FloatField,         
    ql0_in: FloatField,         
    qi0_in: FloatField,         
    th0_in: FloatField,         
    tr0_inout: FloatField_NTracers,      
    tr0: FloatField_NTracers,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    sstr0: FloatField_NTracers,
):
    '''
    University of Washington Shallow Convection Scheme          
                                                                
    Described in Park and Bretherton. 2008. J. Climate :        
                                                                 
    'The University of Washington shallow convection and         
    moist turbulent schemes and their impact on climate         
    simulations with the Community Atmosphere Model'            
                                                                
    Coded in CESM by Sungsu Park. Oct.2005. May.2008.                     
                                                                 
    Coded in GEOS by Nathan Arnold. July 2016.                  
                                                                  
    For general questions, email sungsup@ucar.edu or             
    sungsu@atmos.washington.edu    
                                                                  
    For GEOS-specific questions, email nathan.arnold@nasa.gov                                                         
    

    
    Add description of variables
    '''

    # Start Main Calculation

    with computation(PARALLEL), interval(...):
         #id_exit = False

         zvir = 0.609 # r_H2O/r_air-1
         pmid0 = pmid0_in
         zmid0 = zmid0_in
         u0 = u0_in
         v0 = v0_in
         qv0 = qv0_in
         ql0 = ql0_in
         qi0 = qi0_in

         if dotransport == 1.0:
            n=0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0,0,0][n] = tr0_inout[0,0,0][n]
                n+=1

         #Compute basic thermodynamic variables directly from  
         #input variables for each column                      
         
         # Compute internal environmental variables
         exnmid0 = exnmid0_in
         exnifc0 = exnifc0_in
         t0 = th0_in * exnmid0
         s0 = constants.gravity*zmid0 + constants.cpdry*t0
         qt0 = qv0 + ql0 + qi0
         thl0 = (t0 - constants.latent_heat_vaporization*ql0/constants.cpdry - constants.latent_heat_sublimation*qi0/constants.cpdry) / exnmid0
         thvl0 = (1.0 + zvir*qt0) * thl0

         sstr0 = tr0

    # Compute slopes of environmental variables in each layer
    with computation(PARALLEL), interval(0,1):
         ssthl0 = compute_slope(thl0,pmid0)
         ssqt0 = compute_slope(qt0,pmid0)
         ssu0 = compute_slope(u0,pmid0)
         ssv0 = compute_slope(v0,pmid0)
         #if dotransport == 1.0:
         #   n=0
         #   while n < constants.ncnst:
         #       sstr0 = compute_slope(tr0[n],pmid0)
         #       n+=1

    with computation(PARALLEL), interval(1,-1):
         ssthl0 = compute_above_below_slope(thl0,pmid0)
         ssqt0 = compute_above_below_slope(qt0,pmid0)
         ssu0 = compute_above_below_slope(u0,pmid0)
         ssv0 = compute_above_below_slope(v0,pmid0)
         #if dotransport == 1.0:
         #   n=0
         #   while n < constants.ncnst:
         #       sstr0 = compute_above_below_slope(tr0[n],pmid0)
         #       n+=1

    with computation(PARALLEL), interval(-1,None):
         ssthl0 = ssthl0[0,0,-1]
         ssqt0 = ssqt0[0,0,-1]
         ssu0 = ssu0[0,0,-1]
         ssv0 = ssv0[0,0,-1]
         #if dotransport == 1.0:
         #   n=0
         #   while n < constants.ncnst:
         #       sstr0 = sstr0[0,0,-1][n]
         #       n+=1


class ComputeUwshcu:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:
        
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        #self._slope = stencil_factory.from_dims_halo(
        #    func=slope,
        #    compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #)
        self._compute_uwshcu = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    @staticmethod
    def make_ntracers_quantity_factory(ijk_quantity_factory: QuantityFactory):
        ntracers_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        ntracers_quantity_factory.set_extra_dim_lengths(
            **{
                "ntracers": constants.ncnst,
            }
        )
        return ntracers_quantity_factory

    def __call__(
        self,
        dotransport: Float,          
        exnifc0_in: FloatField,      
        pmid0_in: FloatField,       
        zmid0_in: FloatField,      
        exnmid0_in: FloatField,     
        u0_in: FloatField,          
        v0_in: FloatField,          
        qv0_in: FloatField,       
        ql0_in: FloatField,         
        qi0_in: FloatField,       
        th0_in: FloatField,         
        tr0_inout: FloatField_NTracers,      
        tr0_test: FloatField_NTracers,
        ssthl0_test: FloatField,
        ssqt0_test: FloatField,
        ssu0_test: FloatField,
        ssv0_test: FloatField,
        sstr0_test: FloatField_NTracers,
    ):


        self._compute_uwshcu(
            dotransport=dotransport, 
            exnifc0_in=exnifc0_in, 
            pmid0_in=pmid0_in, 
            zmid0_in=zmid0_in, 
            exnmid0_in=exnmid0_in,
            u0_in=u0_in, 
            v0_in=v0_in, 
            qv0_in=qv0_in, 
            ql0_in=ql0_in, 
            qi0_in=qi0_in, 
            th0_in=th0_in, 
            tr0_inout=tr0_inout,
            tr0=tr0_test,
            ssthl0=ssthl0_test,
            ssqt0=ssqt0_test,
            ssu0=ssu0_test,
            ssv0=ssv0_test,
            sstr0=sstr0_test,
        )
        


      
      
        




            


