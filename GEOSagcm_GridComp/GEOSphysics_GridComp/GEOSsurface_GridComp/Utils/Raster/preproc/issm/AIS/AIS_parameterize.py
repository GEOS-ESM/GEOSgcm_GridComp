import numpy as np
from paterson import paterson
from netCDF4 import Dataset
from xy2ll import xy2ll
from InterpFromGridToMesh import InterpFromGridToMesh
from SetMarineIceSheetBC import SetMarineIceSheetBC
from m1qn3inversion import m1qn3inversion 
from setflowequation import setflowequation
from pathlib import Path

#Name and Coordinate system
md.miscellaneous.name="AIS"
md.mesh.epsg=3031

nsidc_vel = Dataset('/discover/nobackup/agstubbl/ISSM/data/AIS/Antarctica_ice_velocity.nc')
xmin = nsidc_vel.xmin
xmin = float(xmin.lstrip()[0:10])
ymax = nsidc_vel.ymax
ymax = float(ymax.lstrip()[0:10])
spacing = nsidc_vel.spacing
spacing = float((spacing.lstrip())[0:4])
nx = nsidc_vel.nx
ny = nsidc_vel.ny
velx = nsidc_vel['vx'][:].data
vely = nsidc_vel['vy'][:].data
# Build coordinates
x2 = xmin + np.arange(nx + 1) * spacing
y2 = (ymax - ny * spacing) + np.arange(ny + 1) * spacing

# print('   Set observed velocities')
md.initialization.vx = InterpFromGridToMesh(x2, y2, np.flipud(velx), md.mesh.x, md.mesh.y, 0)
md.initialization.vy = InterpFromGridToMesh(x2, y2, np.flipud(vely), md.mesh.x, md.mesh.y, 0)
md.initialization.vz = np.zeros(md.mesh.numberofvertices)
md.initialization.vel = np.sqrt(md.initialization.vx**2 + md.initialization.vy**2)
del velx, vely

# Parameters to change/Try
friction_coefficient = 10  # default [10]
Temp_change = 0           # default [0 K]

# NetCDF Loading
print('   Loading SeaRISE data from NetCDF')
ncdata = Dataset('/discover/nobackup/agstubbl/ISSM/data/AIS/Antarctica_5km_withshelves_v0.75.nc')
x1 = ncdata['x1'][:].data
y1 = ncdata['y1'][:].data
usrf = ncdata['usrf'][:].data[0]
topg = ncdata['topg'][:].data[0]
temp = ncdata['presartm'][:].data[0]
smb = ncdata['presprcp'][:].data[0]
gflux = ncdata['bheatflx_fox'][:].data[0]

# Geometry
print('   Interpolating surface and ice base')
md.geometry.base = InterpFromGridToMesh(x1, y1, topg, md.mesh.x, md.mesh.y, 0)
md.geometry.surface = InterpFromGridToMesh(x1, y1, usrf, md.mesh.x, md.mesh.y, 0)
del usrf, topg

thkmask=ncdata['thkmask'][:].data[0]

##interpolate onto our mesh vertices
groundedice=   InterpFromGridToMesh(x1,y1,thkmask,md.mesh.x,md.mesh.y,0)
groundedice[groundedice<=0]=-1
del thkmask

#fill in the md.mask structure
md.mask.ocean_levelset = groundedice  #ice is grounded for mask equal one
md.mask.ice_levelset = -1*np.ones(np.shape(md.mesh.x))  #ice is present when negatvie

print('   Constructing thickness')
md.geometry.thickness = md.geometry.surface - md.geometry.base

# Ensure hydrostatic equilibrium on ice shelf
di = md.materials.rho_ice / md.materials.rho_water

# Get the node numbers of floating nodes
pos = np.where(md.mask.ocean_levelset < 0)

# Apply flotation criterion
md.geometry.thickness[pos] = 1 / (1 - di) * md.geometry.surface[pos]
md.geometry.base[pos] = md.geometry.surface[pos] - md.geometry.thickness[pos]
md.geometry.hydrostatic_ratio = np.ones(md.mesh.numberofvertices)  

# Set min thickness to 1 meter
pos0 = np.where(md.geometry.thickness <= 1)
md.geometry.thickness[pos0] = 1
md.geometry.surface = md.geometry.thickness + md.geometry.base
md.geometry.bed = md.geometry.base.copy()
md.geometry.bed[pos] = md.geometry.base[pos] - 1000


# Initialization parameters
print('   Interpolating temperatures')
md.initialization.temperature = InterpFromGridToMesh(
    x1, y1, temp, md.mesh.x, md.mesh.y, 0
) + 273.15 + Temp_change

[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, -1)


print('   Set Pressure')
md.initialization.pressure = md.materials.rho_ice * md.constants.g * md.geometry.thickness

print('   Construct ice rheological properties')
md.materials.rheology_n = 3 * np.ones(md.mesh.numberofelements)
md.materials.rheology_B = paterson(md.initialization.temperature)

# Forcings
print('   Interpolating surface mass balance')
mass_balance = InterpFromGridToMesh(x1, y1, smb, md.mesh.x, md.mesh.y, 0)
md.smb.mass_balance = mass_balance * md.materials.rho_water / md.materials.rho_ice

print('   Set geothermal heat flux')
md.basalforcings.geothermalflux = InterpFromGridToMesh(x1, y1, gflux, md.mesh.x, md.mesh.y, 0)

# Friction and inversion set up
print('   Construct basal friction parameters')
md.friction.coefficient = friction_coefficient * np.ones(md.mesh.numberofvertices)
md.friction.p = np.ones(md.mesh.numberofelements)
md.friction.q = np.ones(md.mesh.numberofelements)

# No friction applied on floating ice
pos = np.where(md.mask.ocean_levelset < 0)[0]
md.friction.coefficient[pos] = 0
md.groundingline.migration = 'SubelementMigration'

md.inversion = m1qn3inversion()
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy
md.inversion.vel_obs = md.initialization.vel

print('   Set flow equations')
md = setflowequation(md,'SSA','all')

print('   Set boundary conditions')
md = SetMarineIceSheetBC(md)
md.basalforcings.floatingice_melting_rate = np.zeros(md.mesh.numberofvertices)
md.basalforcings.groundedice_melting_rate = np.zeros(md.mesh.numberofvertices)
md.thermal.spctemperature = md.initialization.temperature
md.masstransport.spcthickness = np.full(md.mesh.numberofvertices, np.nan)




