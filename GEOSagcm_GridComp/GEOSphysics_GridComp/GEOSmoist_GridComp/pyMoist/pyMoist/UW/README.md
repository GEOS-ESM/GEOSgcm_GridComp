## UW Convection Port
## Date: 9/18/24

### GEOS_UW_InterfaceMod.F90
### uwshcu.F90
- The code contains 3 larger subroutines: UW_Run, compute_uwshcu_inv, and compute_uwshcu.

- 'compute_uwshcu' contains several smaller functions and subroutines. including slope and conden. 

- However, many of the smaller functions in compute_uwshcu are not actually triggered in the Fortran. Errors should be raised in the event that this code is triggered under a different configuration. 

- Currently, 'slope' and 'conden' have been ported and verified. Small errors still exist when testing 'conden'. These errors likely stem from errors in the QSat function.

## General notes
- The main calculations in uwshcu.F90 are very long and involve loops along the i and k dimensions. To serialize the data within loops, a workaround has been implemented in the .SER file. This workaround involves defining test variables outside of the loop, creating a savepoint before the loop, and saving the data after the loop has finished.

- In the Fortran, 'slope' and 'conden' uses 2D input variables. These variables can be serialized as 2D arrays, however, they need to be reshaped to shape(i,j,k) prior to the stencil computation. This workaround has been implemented in translate_slope.py and translate_conden.py

- Slope and conden stencils must be rewritten as gt4py functions in order to be used within the larger stencils. Currently, Slope has been rewritten as a function and verified. Conden is in progress.

- 'compute_uwshcu' creates tracer bundles, which have an extra dimension e.g., tr0_inout[576,72,23]. These variables can be serialized as 3D variables, but need to be reshaped to 4D and converted to 4D quantities.


## Next steps
- Finish porting and testing compute_uwshcu
- Port and test compute_uwshcu_inv
- Port and test UW_Run

