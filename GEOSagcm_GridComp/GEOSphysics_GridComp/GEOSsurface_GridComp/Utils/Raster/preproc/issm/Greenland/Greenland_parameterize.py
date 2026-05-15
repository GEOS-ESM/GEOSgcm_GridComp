import numpy as np
from paterson import paterson
from netCDF4 import Dataset
from ll2xy import ll2xy
from xy2ll import xy2ll
from InterpFromGridToMesh import InterpFromGridToMesh
from SetIceSheetBC import SetIceSheetBC
import os
ISSM_DIR = os.getenv('ISSM_DIR')


#Name and Coordinate system
md.miscellaneous.name = 'GreenlandGEOS'
md.mesh.epsg = 3413

# interpolate velocities
ncdata_x = Dataset('/discover/nobackup/agstubbl/ISSM/data/greenland_vel_mosaic200_2017_2018_vx_v02.1.nc', mode='r')
ncdata_y = Dataset('/discover/nobackup/agstubbl/ISSM/data/greenland_vel_mosaic200_2017_2018_vy_v02.1.nc', mode='r')

x1 = np.squeeze(ncdata_x.variables['x'][:].data)
y1 = np.squeeze(ncdata_x.variables['y'][:].data)
velx = np.squeeze(ncdata_x.variables['Band1'][:].data)
vely = np.squeeze(ncdata_y.variables['Band1'][:].data)
ncdata_x.close()
ncdata_y.close()

# set missing data points to zero???
velx[np.abs(velx)>1e9] = 0
vely[np.abs(vely)>1e9] = 0

vx = InterpFromGridToMesh(x1, y1, velx, md.mesh.x, md.mesh.y, 0)
vy = InterpFromGridToMesh(x1, y1, vely, md.mesh.x, md.mesh.y, 0)
speed = np.sqrt(vx**2 + vy**2)

md.initialization.vx = vx
md.initialization.vy = vy
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.vel = speed

md.inversion.vx_obs = vx
md.inversion.vy_obs = vy
md.inversion.vel_obs = speed

# initialize ice thickness
ncdata_H = Dataset('/discover/nobackup/agstubbl/ISSM/data/BedMachineGreenland-v5.nc', mode='r')
H = np.flipud(ncdata_H['thickness'][:].data.astype(np.float64))
surf = np.flipud(ncdata_H['surface'][:].data.astype(np.float64))
bed = np.flipud(ncdata_H['bed'][:].data.astype(np.float64))
x1 = ncdata_H['x'][:].data.astype(np.float64)
y1 = np.flipud(ncdata_H['y'][:].data.astype(np.float64))
ncdata_H.close()

md.geometry.base = InterpFromGridToMesh(x1, y1, bed, md.mesh.x, md.mesh.y, 0)
md.geometry.surface = InterpFromGridToMesh(x1, y1, surf, md.mesh.x, md.mesh.y, 0)

md.geometry.thickness = md.geometry.surface - md.geometry.base

#Set min thickness to 1 meter
pos0 = np.nonzero(md.geometry.thickness <= 0)
md.geometry.thickness[pos0] = 1
md.geometry.surface = md.geometry.thickness + md.geometry.base
 
 #--------------------------------------------------
 ## EDITS
## Project mesh onto old coordinate system temporarily
[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, + 1, 45, 70)
[xi, yi] = ll2xy(md.mesh.lat, md.mesh.long, + 1, 39, 71)
md.mesh.x = xi
md.mesh.y = yi

ncdata = Dataset('/discover/nobackup/agstubbl/ISSM/data/Greenland_5km_dev1.2.nc', mode='r')
x1 = np.squeeze(ncdata.variables['x1'][:].data)
y1 = np.squeeze(ncdata.variables['y1'][:].data)
usrf = np.squeeze(ncdata.variables['usrf'][:].data)
topg = np.squeeze(ncdata.variables['topg'][:].data)
velx = np.squeeze(ncdata.variables['surfvelx'][:].data)
vely = np.squeeze(ncdata.variables['surfvely'][:].data)
temp = np.squeeze(ncdata.variables['airtemp2m'][:].data)
smb = np.squeeze(ncdata.variables['smb'][:].data)
gflux = np.squeeze(ncdata.variables['bheatflx'][:].data)
ncdata.close()

# initialize temperature
md.initialization.temperature = InterpFromGridToMesh(x1, y1, temp, md.mesh.x, md.mesh.y, 0) + 273.15

# impose observed temperature on surface
md.thermal.spctemperature = md.initialization.temperature
md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))

# initialize surface mass balance (zero smb example)
md.smb.mass_balance = InterpFromGridToMesh(x1, y1, smb, md.mesh.x, md.mesh.y, 0)
md.smb.mass_balance = md.smb.mass_balance * md.materials.rho_water / md.materials.rho_ice

# initialize basal friction
md.friction.coefficient = 30 * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(md.mask.ocean_levelset < 0)
md.friction.coefficient[pos] = 0  #no friction applied on floating ice
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# initialize ice rheology
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))
md.materials.rheology_B = paterson(md.initialization.temperature)

# set geothermal heat flux
md.basalforcings.geothermalflux = InterpFromGridToMesh(x1, y1, gflux, md.mesh.x, md.mesh.y, 0)

# set other boundary conditions
md.mask.ice_levelset[np.nonzero(md.mesh.vertexonboundary == 1)] = 0
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))

# initialize pressure
md.initialization.pressure = md.materials.rho_ice * md.constants.g * md.geometry.thickness

# initialize single point constraints
md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

## Re-project mesh onto coordinate system
[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, + 1, 39, 71)
[xi, yi] = ll2xy(md.mesh.lat, md.mesh.long, + 1, 45, 70)
md.mesh.x = xi
md.mesh.y = yi

md = SetIceSheetBC(md) 
