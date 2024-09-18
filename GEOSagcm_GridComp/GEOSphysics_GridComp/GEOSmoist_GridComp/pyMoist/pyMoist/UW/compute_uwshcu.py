import copy

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD
import pyMoist.pyMoist_constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Float, Int
from ndsl import StencilFactory, QuantityFactory
from pyMoist.types import FloatField_NTracers


@gtscript.function
def slope(kmask: Float, field: Float, field_above: Float, field_below: Float, p0: Float, p0_above: Float, p0_below: Float):
    """
    Calculates slope of a given field.

    Parameters:
    kmask (Float): K-level mask (e.g., 1 for k=0, 2 for k=1,71, 3 for k=72).
    field (Float): Field of interest.
    field_above (Float): 1 k-level above field (e.g., field[0,0,1]).
    field_below (Float): 1 k-level below field (e.g., field[0,0,-1]).
    p0 (Float): Pressure 
    p0_above (Float): 1 k-level above p0 (e.g., p0[0,0,1]).
    p0_below (Float): 1 k-level below p0 (e.g., p0[0,0,-1]).

    Returns:
    Slope: Slope of the field of interest.
    """
    if kmask == 1:
        value = (field_above - field) / (p0_above - p0)
        if value > 0.0:
            slope = max(0.0, value)
        else:
            slope =  min(0.0, value)
    elif kmask == 2:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))
    else:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
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
    kmask: FloatField,
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

    '''
    Add description of variables
    '''
    
    # Start Main Calculation
    with computation(PARALLEL), interval(0,1):
         pmid0 = pmid0_in
         pmid0_above = pmid0_in[0,0,1]
         u0 = u0_in
         u0_above = u0_in[0,0,1]
         v0 = v0_in
         v0_above = v0_in[0,0,1]
         qv0 = qv0_in
         qv0_above = qv0_in[0,0,1]
         ql0 = ql0_in
         ql0_above = ql0_in[0,0,1]
         qi0 = qi0_in
         qi0_above = qi0_in[0,0,1]
         qt0 = qv0 + ql0 + qi0
         qt0_above = qv0_above + ql0_above + qi0_above
         exnmid0 = exnmid0_in
         exnmid0_above = exnmid0_in[0,0,1]
         t0 = th0_in * exnmid0
         t0_above = th0_in[0,0,1] * exnmid0_above
         thl0 = (t0 - constants.latent_heat_vaporization*ql0/constants.cpdry - constants.latent_heat_sublimation*qi0/constants.cpdry) / exnmid0
         thl0_above = (t0_above - constants.latent_heat_vaporization*ql0_above/constants.cpdry - constants.latent_heat_sublimation*qi0_above/constants.cpdry) / exnmid0_above

         if dotransport == 1.0:
            n=0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0,0,0][n] = tr0_inout[0,0,0][n]
                n+=1

         ssthl0 = slope(kmask,thl0,thl0_above,thl0_above,pmid0,pmid0_above,pmid0_above)
         ssqt0 = slope(kmask,qt0,qt0_above,qt0_above,pmid0,pmid0_above,pmid0_above)
         ssu0 = slope(kmask,u0,u0_above,u0_above,pmid0,pmid0_above,pmid0_above)
         ssv0 = slope(kmask,v0,v0_above,v0_above,pmid0,pmid0_above,pmid0_above)

         if dotransport == 1.0:
             # Raise error if constants.ncnst != ncnst
             sstr0[0,0,0][0] = slope(kmask,tr0[0,0,0][0],tr0[0,0,1][0],tr0[0,0,1][0],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][1] = slope(kmask,tr0[0,0,0][1],tr0[0,0,1][1],tr0[0,0,1][1],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][2] = slope(kmask,tr0[0,0,0][2],tr0[0,0,1][2],tr0[0,0,1][2],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][3] = slope(kmask,tr0[0,0,0][3],tr0[0,0,1][3],tr0[0,0,1][3],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][4] = slope(kmask,tr0[0,0,0][4],tr0[0,0,1][4],tr0[0,0,1][4],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][5] = slope(kmask,tr0[0,0,0][5],tr0[0,0,1][5],tr0[0,0,1][5],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][6] = slope(kmask,tr0[0,0,0][6],tr0[0,0,1][6],tr0[0,0,1][6],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][7] = slope(kmask,tr0[0,0,0][7],tr0[0,0,1][7],tr0[0,0,1][7],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][8] = slope(kmask,tr0[0,0,0][8],tr0[0,0,1][8],tr0[0,0,1][8],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][9] = slope(kmask,tr0[0,0,0][9],tr0[0,0,1][9],tr0[0,0,1][9],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][10] = slope(kmask,tr0[0,0,0][10],tr0[0,0,1][10],tr0[0,0,1][10],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][11] = slope(kmask,tr0[0,0,0][11],tr0[0,0,1][11],tr0[0,0,1][11],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][12] = slope(kmask,tr0[0,0,0][12],tr0[0,0,1][12],tr0[0,0,1][12],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][13] = slope(kmask,tr0[0,0,0][13],tr0[0,0,1][13],tr0[0,0,1][13],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][14] = slope(kmask,tr0[0,0,0][14],tr0[0,0,1][14],tr0[0,0,1][14],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][15] = slope(kmask,tr0[0,0,0][15],tr0[0,0,1][15],tr0[0,0,1][15],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][16] = slope(kmask,tr0[0,0,0][16],tr0[0,0,1][16],tr0[0,0,1][16],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][17] = slope(kmask,tr0[0,0,0][17],tr0[0,0,1][17],tr0[0,0,1][17],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][18] = slope(kmask,tr0[0,0,0][18],tr0[0,0,1][18],tr0[0,0,1][18],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][19] = slope(kmask,tr0[0,0,0][19],tr0[0,0,1][19],tr0[0,0,1][19],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][20] = slope(kmask,tr0[0,0,0][20],tr0[0,0,1][20],tr0[0,0,1][20],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][21] = slope(kmask,tr0[0,0,0][21],tr0[0,0,1][21],tr0[0,0,1][21],pmid0,pmid0_above,pmid0_above)
             sstr0[0,0,0][22] = slope(kmask,tr0[0,0,0][22],tr0[0,0,1][22],tr0[0,0,1][22],pmid0,pmid0_above,pmid0_above)


    with computation(FORWARD), interval(1,-1):
         pmid0 = pmid0_in
         pmid0_above = pmid0_in[0,0,1]
         pmid0_below = pmid0_in[0,0,-1]
         zmid0 = zmid0_in
         u0 = u0_in
         u0_above = u0_in[0,0,1]
         u0_below = u0_in[0,0,-1]
         v0 = v0_in
         v0_above = v0_in[0,0,1]
         v0_below = v0_in[0,0,-1]
         qv0 = qv0_in
         qv0_above = qv0_in[0,0,1]
         qv0_below = qv0_in[0,0,-1]
         ql0 = ql0_in
         ql0_above = ql0_in[0,0,1]
         ql0_below = ql0_in[0,0,-1]
         qi0 = qi0_in
         qi0_above = qi0_in[0,0,1]
         qi0_below = qi0_in[0,0,-1]

         if dotransport == 1.0:
            n=0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0,0,0][n] = tr0_inout[0,0,0][n]
                n+=1

         exnmid0 = exnmid0_in
         exnmid0_above = exnmid0_in[0,0,1]
         exnmid0_below = exnmid0_in[0,0,-1]
         t0 = th0_in * exnmid0
         t0_above = th0_in[0,0,1] * exnmid0_above
         t0_below = th0_in[0,0,-1] * exnmid0_below
         qt0 = qv0 + ql0 + qi0
         qt0_above = qv0_above + ql0_above + qi0_above
         qt0_below = qv0_below + ql0_below + qi0_below
         thl0 = (t0 - constants.latent_heat_vaporization*ql0/constants.cpdry - constants.latent_heat_sublimation*qi0/constants.cpdry) / exnmid0
         thl0_above = (t0_above - constants.latent_heat_vaporization*ql0_above/constants.cpdry - constants.latent_heat_sublimation*qi0_above/constants.cpdry) / exnmid0_above
         thl0_below = (t0_below - constants.latent_heat_vaporization*ql0_below/constants.cpdry - constants.latent_heat_sublimation*qi0_below/constants.cpdry) / exnmid0_below
        
         ssthl0 = slope(kmask,thl0,thl0_above,thl0_below,pmid0,pmid0_above,pmid0_below)
         ssqt0 = slope(kmask,qt0,qt0_above,qt0_below,pmid0,pmid0_above,pmid0_below)
         ssu0 = slope(kmask,u0,u0_above,u0_below,pmid0,pmid0_above,pmid0_below)
         ssv0 = slope(kmask,v0,v0_above,v0_below,pmid0,pmid0_above,pmid0_below)

      
         if dotransport == 1.0:
             # Raise error if constants.ncnst != ncnst
             sstr0[0,0,0][0] = slope(kmask,tr0[0,0,0][0],tr0[0,0,1][0],tr0[0,0,-1][0],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][1] = slope(kmask,tr0[0,0,0][1],tr0[0,0,1][1],tr0[0,0,-1][1],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][2] = slope(kmask,tr0[0,0,0][2],tr0[0,0,1][2],tr0[0,0,-1][2],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][3] = slope(kmask,tr0[0,0,0][3],tr0[0,0,1][3],tr0[0,0,-1][3],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][4] = slope(kmask,tr0[0,0,0][4],tr0[0,0,1][4],tr0[0,0,-1][4],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][5] = slope(kmask,tr0[0,0,0][5],tr0[0,0,1][5],tr0[0,0,-1][5],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][6] = slope(kmask,tr0[0,0,0][6],tr0[0,0,1][6],tr0[0,0,-1][6],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][7] = slope(kmask,tr0[0,0,0][7],tr0[0,0,1][7],tr0[0,0,-1][7],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][8] = slope(kmask,tr0[0,0,0][8],tr0[0,0,1][8],tr0[0,0,-1][8],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][9] = slope(kmask,tr0[0,0,0][9],tr0[0,0,1][9],tr0[0,0,-1][9],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][10] = slope(kmask,tr0[0,0,0][10],tr0[0,0,1][10],tr0[0,0,-1][10],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][11] = slope(kmask,tr0[0,0,0][11],tr0[0,0,1][11],tr0[0,0,-1][11],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][12] = slope(kmask,tr0[0,0,0][12],tr0[0,0,1][12],tr0[0,0,-1][12],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][13] = slope(kmask,tr0[0,0,0][13],tr0[0,0,1][13],tr0[0,0,-1][13],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][14] = slope(kmask,tr0[0,0,0][14],tr0[0,0,1][14],tr0[0,0,-1][14],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][15] = slope(kmask,tr0[0,0,0][15],tr0[0,0,1][15],tr0[0,0,-1][15],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][16] = slope(kmask,tr0[0,0,0][16],tr0[0,0,1][16],tr0[0,0,-1][16],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][17] = slope(kmask,tr0[0,0,0][17],tr0[0,0,1][17],tr0[0,0,-1][17],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][18] = slope(kmask,tr0[0,0,0][18],tr0[0,0,1][18],tr0[0,0,-1][18],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][19] = slope(kmask,tr0[0,0,0][19],tr0[0,0,1][19],tr0[0,0,-1][19],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][20] = slope(kmask,tr0[0,0,0][20],tr0[0,0,1][20],tr0[0,0,-1][20],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][21] = slope(kmask,tr0[0,0,0][21],tr0[0,0,1][21],tr0[0,0,-1][21],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][22] = slope(kmask,tr0[0,0,0][22],tr0[0,0,1][22],tr0[0,0,-1][22],pmid0,pmid0_above,pmid0_below)

    with computation(PARALLEL), interval(-1,None):
         pmid0 = pmid0_in[0,0,-1]
         pmid0_above = pmid0_in
         pmid0_below = pmid0_in[0,0,-2]
         u0 = u0_in[0,0,-1]
         u0_above = u0_in
         u0_below = u0_in[0,0,-2]
         v0 = v0_in[0,0,-1]
         v0_above = v0_in
         v0_below = v0_in[0,0,-2]
         qv0 = qv0_in[0,0,-1]
         qv0_above = qv0_in
         qv0_below = qv0_in[0,0,-2]
         ql0 = ql0_in[0,0,-1]
         ql0_above = ql0_in
         ql0_below = ql0_in[0,0,-2]
         qi0 = qi0_in[0,0,-1]
         qi0_above = qi0_in
         qi0_below = qi0_in[0,0,-2]

         exnmid0 = exnmid0_in[0,0,-1]
         exnmid0_above = exnmid0_in
         exnmid0_below = exnmid0_in[0,0,-2]
         t0 = th0_in[0,0,-1] * exnmid0
         t0_above = th0_in * exnmid0_above
         t0_below = th0_in[0,0,-2] * exnmid0_below
         qt0 = qv0 + ql0 + qi0
         qt0_above = qv0_above + ql0_above + qi0_above
         qt0_below = qv0_below + ql0_below + qi0_below
         thl0 = (t0 - constants.latent_heat_vaporization*ql0/constants.cpdry - constants.latent_heat_sublimation*qi0/constants.cpdry) / exnmid0
         thl0_above = (t0_above - constants.latent_heat_vaporization*ql0_above/constants.cpdry - constants.latent_heat_sublimation*qi0_above/constants.cpdry) / exnmid0_above
         thl0_below = (t0_below - constants.latent_heat_vaporization*ql0_below/constants.cpdry - constants.latent_heat_sublimation*qi0_below/constants.cpdry) / exnmid0_below

         if dotransport == 1.0:
            n=0
            # Loop over tracers
            while n < constants.ncnst:
                tr0[0,0,0][n] = tr0_inout[0,0,0][n]
                n+=1

         ssthl0 = slope(kmask,thl0,thl0_above,thl0_below,pmid0,pmid0_above,pmid0_below)
         ssqt0 = slope(kmask,qt0,qt0_above,qt0_below,pmid0,pmid0_above,pmid0_below)
         ssu0 = slope(kmask,u0,u0_above,u0_below,pmid0,pmid0_above,pmid0_below)
         ssv0 = slope(kmask,v0,v0_above,v0_below,pmid0,pmid0_above,pmid0_below)

         if dotransport == 1.0:
             # Raise error if constants.ncnst != ncnst
             sstr0[0,0,0][0] = slope(kmask,tr0[0,0,-1][0],tr0[0,0,0][0],tr0[0,0,-2][0],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][1] = slope(kmask,tr0[0,0,-1][1],tr0[0,0,0][1],tr0[0,0,-2][1],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][2] = slope(kmask,tr0[0,0,-1][2],tr0[0,0,0][2],tr0[0,0,-2][2],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][3] = slope(kmask,tr0[0,0,-1][3],tr0[0,0,0][3],tr0[0,0,-2][3],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][4] = slope(kmask,tr0[0,0,-1][4],tr0[0,0,0][4],tr0[0,0,-2][4],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][5] = slope(kmask,tr0[0,0,-1][5],tr0[0,0,0][5],tr0[0,0,-2][5],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][6] = slope(kmask,tr0[0,0,-1][6],tr0[0,0,0][6],tr0[0,0,-2][6],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][7] = slope(kmask,tr0[0,0,-1][7],tr0[0,0,0][7],tr0[0,0,-2][7],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][8] = slope(kmask,tr0[0,0,-1][8],tr0[0,0,0][8],tr0[0,0,-2][8],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][9] = slope(kmask,tr0[0,0,-1][9],tr0[0,0,0][9],tr0[0,0,-2][9],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][10] = slope(kmask,tr0[0,0,-1][10],tr0[0,0,0][10],tr0[0,0,-2][10],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][11] = slope(kmask,tr0[0,0,-1][11],tr0[0,0,0][11],tr0[0,0,-2][11],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][12] = slope(kmask,tr0[0,0,-1][12],tr0[0,0,0][12],tr0[0,0,-2][12],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][13] = slope(kmask,tr0[0,0,-1][13],tr0[0,0,0][13],tr0[0,0,-2][13],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][14] = slope(kmask,tr0[0,0,-1][14],tr0[0,0,0][14],tr0[0,0,-2][14],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][15] = slope(kmask,tr0[0,0,-1][15],tr0[0,0,0][15],tr0[0,0,-2][15],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][16] = slope(kmask,tr0[0,0,-1][16],tr0[0,0,0][16],tr0[0,0,-2][16],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][17] = slope(kmask,tr0[0,0,-1][17],tr0[0,0,0][17],tr0[0,0,-2][17],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][18] = slope(kmask,tr0[0,0,-1][18],tr0[0,0,0][18],tr0[0,0,-2][18],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][19] = slope(kmask,tr0[0,0,-1][19],tr0[0,0,0][19],tr0[0,0,-2][19],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][20] = slope(kmask,tr0[0,0,-1][20],tr0[0,0,0][20],tr0[0,0,-2][20],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][21] = slope(kmask,tr0[0,0,-1][21],tr0[0,0,0][21],tr0[0,0,-2][21],pmid0,pmid0_above,pmid0_below)
             sstr0[0,0,0][22] = slope(kmask,tr0[0,0,-1][22],tr0[0,0,0][22],tr0[0,0,-2][22],pmid0,pmid0_above,pmid0_below)




class ComputeUwshcu:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        ncnst: Int,
    ) -> None:
        """
        Initialize the ComputeUwshcu class.

        Parameters:
        stencil_factory (StencilFactory): Factory for creating stencil computations.
        quantity_factory (QuantityFactory): Factory for creating quantities.
        ncnst (Int): Number of tracers.

        Raises:
        NotImplementedError: If the number of tracers is not equal to the expected number.
        """
        if constants.ncnst != ncnst:
            raise NotImplementedError(
                f"Coding limitation: 23 tracers are expected, getting {ncnst}"
            )
        
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._compute_uwshcu = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._k_mask = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    if k == 0:
                        self._k_mask.view[i, j, k] = 1
                    elif k == 71:
                        self._k_mask.view[i, j, k] = 3
                    else:
                        self._k_mask.view[i, j, k] = 2


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
            kmask=self._k_mask,
        )
        


      
      
        




            


