import devpath
import sys,os
import numpy as np
from triangle import triangle
from model import *
from interpMIPkitAntarctica import interpMIPkitAntarctica
from bamg import bamg
from clusters.discover_geos import export_discover
from xy2ll import xy2ll

h_max = sys.argv[1] if len(sys.argv) > 1 else 24000
h_min = sys.argv[2] if len(sys.argv) > 2 else 2000
h_max = float(h_max)
h_min = float(h_min)

print(f'h_max: {h_max}')
print(f'h_min: {h_min}')

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
# Generate initial ~uniform mesh (resolution = h_max)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/landice/issm/v1/AIS/AntarcticaOutline.exp', h_max)

mipkit_path = '/discover/nobackup/agstubbl/ISSM/data/AIS/AntarcticaObsISMIP7-v1.2.nc'

vx = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'vx_mosaic',ncfile=mipkit_path,default_value=0)
vy = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'vy_mosaic',ncfile=mipkit_path,default_value=0)
v = np.sqrt(vx**2 + vy**2)

# # # # Refine mesh according to flow speed
md = bamg(md, 'hmax', h_max, 'hmin', h_min, 'gradation', 1.4, 'field', v, 'err', 8)

# update lat, long
[md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, -1)

# export mesh
export_discover(md, './netcdfs/AIS_mesh.nc')
