import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL
import pyMoist.pyMoist_constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, IntField, Int 
from ndsl import StencilFactory, QuantityFactory 

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


def compute_uwshcu(
    dotransport: Int,           # Transport tracers [1 true]
    exnifc0_in: FloatField,     # Exner function at interfaces 
    pmid0_in: FloatField,       # Environmental pressure at midpoints [Pa] 
    zmid0_in: FloatField,       # Environmental height at midpoints [m] 
    exnmid0_in: FloatField,     # Exner function at midpoints
    u0_in: FloatField,          # Environmental zonal wind [m/s] 
    v0_in: FloatField,          # Environmental meridional wind [m/s] 
    qv0_in: FloatField,         # Environmental specific humidity
    ql0_in: FloatField,         # Environmental liquid water specific humidity 
    qi0_in: FloatField,         # Environmental ice specific humidity
    th0_in: FloatField,         # Environmental potential temperature [K]
    tr0_inout: FloatField,      # Environmental tracers [ #, kg/kg ]
    tr0: FloatField,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    sstr0: FloatField
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
    '''

    # Start Main Calculation

    with computation(PARALLEL), interval(...):
         #id_exit = False

         zvir = Float(0.609) # r_H2O/r_air-1
         pmid0 = pmid0_in
         zmid0 = zmid0_in
         u0 = u0_in
         v0 = v0_in
         qv0 = qv0_in
         ql0 = ql0_in
         qi0 = qi0_in

         #if dotransport == 1:
         #   tr0 = tr0_inout

         '''
         Compute basic thermodynamic variables directly from  
         input variables for each column     
         '''                 
         
         # Compute internal environmental variables
         exnmid0 = exnmid0_in
         exnifc0 = exnifc0_in
         t0 = th0_in * exnmid0
         s0 = constants.gravity*zmid0 + constants.cpdry*t0
         qt0 = qv0 + ql0 + qi0
         thl0 = (t0 - constants.latent_heat_vaporization*ql0/constants.cpdry - constants.latent_heat_sublimation*qi0/constants.cpdry) / exnmid0
         thvl0 = (1.0 + zvir*qt0) * thl0


         # Compute slopes of environmental variables in each layer
         ssthl0 = slope(thl0, pmid0)
         ssqt0 = slope(qt0 , pmid0)
         ssu0 = slope(u0  , pmid0)
         ssv0 = slope(v0  , pmid0)

         #if dotransport == 1:
         #   sstr0 = slope(tr0, pmid0)

      
      
        




            


