import devpath
import numpy as np
from interpMIPkitGreenland import interpMIPkitGreenland
from reinitializelevelset import reinitializelevelset
import sethydrostaticmask
from effectivepressure import effectivepressure
from averaging import averaging
from drivingstress import drivingstress
from removeponds import removeponds
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from friction import friction
from setflowequation import setflowequation
from cuffey import cuffey
from loadmodel import loadmodel
from clusters.discover_geos import export_discover

md = loadmodel('./netcdfs/GRIS_mesh.nc')

#Name and Coordinate system
md.miscellaneous.name = "GRIS"
md.mesh.epsg = 3413

# path to ISMIP7 data 
mipkit_path = '/discover/nobackup/agstubbl/ISSM/data/GRIS/GreenlandObsISMIP7-v1.3.nc'

# Setting mask:
mask = interpMIPkitGreenland(md.mesh.x,md.mesh.y,'mask',ncfile=mipkit_path,default_value=0)
md.mask.ice_levelset= -1*np.ones(md.mesh.numberofvertices)
pos = np.where(mask<1)[0]
md.mask.ice_levelset[pos]=1

#Tricky part here: we want to offset the mask by one element so that we don't end up with a cliff at the transition
elements = md.mesh.elements.astype(int) - 1                  # 0-based, int
levelset = np.asarray(md.mask.ice_levelset, dtype=float).ravel().copy()

# offset the mask by one element so we don't end up with a cliff at the transition
pos = np.flatnonzero(levelset[elements].max(axis=1) > 0)
levelset[elements[pos, :]] = 1

md.mask.ice_levelset = levelset.reshape(np.shape(md.mask.ice_levelset))
md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset)

# Reading geometry:
md.geometry.surface = interpMIPkitGreenland(md.mesh.x,md.mesh.y,'surface_gimp',ncfile=mipkit_path,default_value=0)
md.geometry.bed     = interpMIPkitGreenland(md.mesh.x,md.mesh.y,'bed',ncfile=mipkit_path,default_value=0)

# Determining ice thickness based on floatation:
min_thickness = 1.0

md.geometry.base = md.geometry.bed.copy()  

# using the corrected surface
floatation_base = (md.geometry.surface * md.materials.rho_ice
                   / (md.materials.rho_ice - md.materials.rho_water))

# grounded ice
pos = floatation_base <= md.geometry.bed
md.geometry.base[pos] = md.geometry.bed[pos]

# floating ice
pos = floatation_base > md.geometry.bed
md.geometry.base[pos] = floatation_base[pos]

md.geometry.thickness = md.geometry.surface - md.geometry.base

# dealing with rocks or ocean part
pos = md.geometry.thickness < min_thickness
md.geometry.thickness[pos] = min_thickness
md.geometry.surface[pos] = md.geometry.thickness[pos] + md.geometry.base[pos]

# Add some correction for "outliers"
# ice-free areas (levelset > 0): force surface to base + min_thickness
pos = (md.mask.ice_levelset > 0)
md.geometry.surface[pos] = md.geometry.base[pos] + min_thickness

md.geometry.thickness = md.geometry.surface - md.geometry.base

# enforce minimum thickness everywhere
pos = md.geometry.thickness <= min_thickness
md.geometry.surface[pos] = md.geometry.base[pos] + min_thickness

md.geometry.thickness = md.geometry.surface - md.geometry.base

md.masstransport.min_thickness = min_thickness

# Set ocean levelset now that geometry is set
md = sethydrostaticmask.setmask(md)

# remove any little isolated 'ocean' ponds
md = removeponds(md, min_area=1e3,min_vertices=10,raise_bed=1)

# Reading velocities:
md.inversion.vx_obs =  interpMIPkitGreenland(md.mesh.x, md.mesh.y, 'vx_mosaic',mipkit_path,default_value=0)
md.inversion.vy_obs =  interpMIPkitGreenland(md.mesh.x, md.mesh.y, 'vy_mosaic',mipkit_path,default_value=0)
pos = np.isnan(md.inversion.vx_obs) | np.isnan(md.inversion.vy_obs)
md.inversion.vx_obs[pos] = 0
md.inversion.vy_obs[pos] = 0
md.inversion.vel_obs = np.sqrt(md.inversion.vx_obs**2 + md.inversion.vy_obs**2)
md.initialization.vx  = md.inversion.vx_obs.copy()
md.initialization.vy  = md.inversion.vy_obs.copy()
md.initialization.vz  = np.zeros(md.mesh.numberofvertices)
md.initialization.vel = md.inversion.vel_obs.copy()

# Basal friction initialization (based on smoothed driving stress)
Neff = np.fmax(0, np.asarray(effectivepressure(md)).flatten())
_, _, taub = drivingstress(md)
taub = np.asarray(averaging(md, taub, 3)).flatten()

vb = md.initialization.vel / md.constants.yts

with np.errstate(divide='ignore', invalid='ignore'):
    C = np.fmin(400, np.sqrt(taub / (Neff * vb)))

# extrapolate Cs where undefined
flags = (vb == 0) | (Neff < 1e5) | np.isnan(C)
good = ~flags

x = np.asarray(md.mesh.x).flatten()
y = np.asarray(md.mesh.y).flatten()

pts = np.column_stack((x[good], y[good]))
F = LinearNDInterpolator(pts, C[good])
C_new = F(x[flags], y[flags])

# linear interpolation returns NaN outside the convex hull, so fall back to nearest
bad = np.isnan(C_new)
if bad.any():
    C_new[bad] = NearestNDInterpolator(pts, C[good])(x[flags][bad], y[flags][bad])
C[flags] = C_new

print(f"      → Interpolated {flags.sum()} friction coefficients")

# Set Budd's friction law
md.friction = friction()
md.friction.coefficient = C # try constant for test: 50*np.ones(md.mesh.numberofvertices) 
md.friction.p = np.ones(md.mesh.numberofelements)
md.friction.q = np.ones(md.mesh.numberofelements)

# Flow law (assume ice is at -10°C for now)
md=setflowequation(md,'SSA','all')
md.materials.rheology_n = 3*np.ones(md.mesh.numberofelements)
md.materials.rheology_B = cuffey(273.15 - 10)*np.ones(md.mesh.numberofvertices)

# Set boundary conditions
md.stressbalance.spcvx=np.nan*np.ones(md.mesh.numberofvertices)
md.stressbalance.spcvy=np.nan*np.ones(md.mesh.numberofvertices)
md.stressbalance.spcvz=np.nan*np.ones(md.mesh.numberofvertices)
md.stressbalance.referential=np.nan*np.ones((md.mesh.numberofvertices,6))
md.stressbalance.loadingforce=np.zeros((md.mesh.numberofvertices,3))
pos=np.where((md.mask.ice_levelset<0) & (md.mesh.vertexonboundary))[0]
md.stressbalance.spcvx[pos]=md.initialization.vx[pos]
md.stressbalance.spcvy[pos]=md.initialization.vy[pos]
md.stressbalance.spcvz[pos]=0

# transient options
md.transient.ismovingfront = 1
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.isgroundingline = 1
md.transient.amr_frequency = 0

# BC for mass transport
onboundary = np.asarray(md.mesh.vertexonboundary).flatten().astype(bool)
spcthickness = np.full((md.mesh.numberofvertices, 1), np.nan)
spcthickness[onboundary, 0] = np.asarray(md.geometry.thickness).flatten()[onboundary]
md.masstransport.spcthickness = spcthickness
md.masstransport.stabilization = 5  # SUPG

# Levelset options (https://gmd.copernicus.org/articles/17/6227/2024/)
spclevelset = np.full((md.mesh.numberofvertices, 1), np.nan)
spclevelset[onboundary, 0] = np.asarray(md.mask.ice_levelset).flatten()[onboundary]
md.levelset.spclevelset = spclevelset
md.levelset.reinit_frequency = 20  # Should be more than 10
md.levelset.stabilization = 5      # SUPG

# basal forcings
md.basalforcings.groundedice_melting_rate = np.zeros(md.mesh.numberofvertices)
md.basalforcings.floatingice_melting_rate = np.zeros(md.mesh.numberofvertices)

# zero frontal melt rate and calving rate, placeholder zero smb
md.frontalforcings.meltingrate = np.zeros(md.mesh.numberofvertices)
md.calving.calvingrate = np.zeros(md.mesh.numberofvertices)
md.smb.mass_balance = np.zeros(md.mesh.numberofvertices) 
md.masstransport.spcthickness = np.nan * np.ones(md.mesh.numberofvertices)

# export parameterization
export_discover(md, './netcdfs/GRIS_parameterization.nc')