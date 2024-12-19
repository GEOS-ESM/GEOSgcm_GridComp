## UW Convection Port
## Date: 9/26/24

### GEOS_UW_InterfaceMod.F90
### uwshcu.F90
- The code contains 3 larger subroutines: UW_Run, compute_uwshcu_inv, and compute_uwshcu.

- 'compute_uwshcu' contains several smaller functions and subroutines, including slope, conden, roots, fluxbelowinv, etc. 

- However, many of the smaller functions in compute_uwshcu are not always triggered in the Fortran. Errors should be raised in the event that this code is triggered under a different configuration. 

- Currently, slope, conden, roots, compute_alpha, compute_mumin2, compute_ppen, getbuoy, single_cin, fluxbelowinv, positive_moisture_single have been ported and verified. To do still: qsinvert

- 10/7/24 compute_uwshcu has been ported and passes for all but 3 variables linked to the conden function. These variables are defined as real*8 (64-bit) in uwshcu.F90. This causes mixed precision and will lead to problems. For now, we plan to write something into the translate test to 'PASS' and we will move on.

## Other notes
- The main calculations in uwshcu.F90 are very long and involve loops along the i and k dimensions. To serialize the data within loops, data_buffered and data_append can be used.

- In the Fortran, 'slope' and 'conden' use 2D input variables. These variables can be serialized as 2D arrays, however, they need to be reshaped to shape(i,j,k) prior to the stencil computation. This workaround has been implemented in translate_slope.py and translate_conden.py

- Slope and conden stencils must be rewritten as gt4py functions in order to be used within the larger stencils. Currently, Slope has been rewritten as a function and verified. Conden is in progress. These stencils have been deleted (conden.py and slope.py).

- 'compute_uwshcu' creates tracer bundles, which have an extra dimension e.g., tr0_inout[576,72,23]. These variables can be serialized as 3D variables, but need to be reshaped to 4D and converted to 4D quantities.

-'pifc0_in' and some other unused variables in 'compute_uwshcu' have 73 k-levels. These variables should be initialized as quantities using Z_INTERFACE_DIM instead of Z_DIM.

- 9/26/24 Ran into some translate test issues that stem from the gt4py branch. 'compute_uwshcu' was not verifying on gt4py cast_to_int branch. Switching to unstable/develop led to much better numerics: 
Better numerics: 60294592b47910816915adc1f51e7c715481405d
Probable change: e0fb2a2410909a9ce63d0c15f33325dc3577eb70
UW original port : cfc8a721829fdaad245715def3aa8721b89b1d09

- 10/4/24 Fortran constant MAPL_UNDEF is used in compute_uwshcu, which is currently assigned the value 1E15. We plan to keep this constant as is for now to match 11.5.2 GEOS. However, we plan to change this at some point.


## Next steps
- Finish porting smaller functions/stencils in uwshcu
- Port and test compute_uwshcu
- Port and test compute_uwshcu_inv
- Port and test UW_Run

