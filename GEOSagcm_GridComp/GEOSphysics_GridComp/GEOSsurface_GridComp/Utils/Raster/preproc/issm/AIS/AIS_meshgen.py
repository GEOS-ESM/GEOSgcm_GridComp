import devpath
import sys,os
import numpy as np
from triangle import triangle
from model import *
from interpMIPkitAntarctica import interpMIPkitAntarctica
from bamg import bamg
from clusters.discover_geos import export_discover

h_max = sys.argv[1] if len(sys.argv) > 1 else 24000
h_min = sys.argv[2] if len(sys.argv) > 2 else 2000
h_max = float(h_max)
h_min = float(h_min)

print(f'h_max: {h_max}')
print(f'h_min: {h_min}')

if not os.path.exists('./netcdfs'):
    os.mkdir('./netcdfs')

# Step 1: Mesh generation
# Generate initial ~uniform mesh (resolution = 20000 m)
# project mesh onto new coordinate system
md = triangle(model(), '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/landice/issm/v1/AIS/AntarcticaOutline.exp', 20000)

vx = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'vx_mosaic',ncfile=mipkit_path,default_value=0)
vy = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'vy_mosaic',ncfile=mipkit_path,default_value=0)
v = np.sqrt(vx**2 + vy**2)

# # # # Refine mesh according to flow speed
md = bamg(md, 'hmax', 24000, 'hmin', 2000, 'gradation', 1.4, 'field', v, 'err', 8)

# export mesh
export_discover(md, './netcdfs/AIS_mesh.nc')
