import devpath
import numpy as np
import os,sys
from pathlib import Path
from contextlib import redirect_stdout
from netCDF4 import Dataset
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

from model import *
from loadmodel import loadmodel

num_nodes = np.array([ ])
num_edges = np.array([ ])
mean_edges = np.array([ ])

# node coordinates
nodeCoords_lon = np.array([])
nodeCoords_lat = np.array([])

# element connectivity
elementConn_n1 = np.array([])
elementConn_n2 = np.array([])
elementConn_n3 = np.array([])


# find all subdirectories with python code
ROOT = Path(".").resolve()
EXCLUDE = (ROOT / "utils_issm").resolve()

models = sorted({
    str(p.parent)
    for p in ROOT.rglob("*.py")
    if EXCLUDE not in p.resolve().parents
    and p.parent != ROOT
})

# top-level glacier names
model_names = sorted({Path(p).relative_to(ROOT).parts[0] for p in models})

i=0
for model in models:
   with open(os.devnull, "w") as f: 
      with redirect_stdout(f): 
         md = loadmodel(f'{model}/netcdfs/{model_names[i]}_initialization.nc')

   v1_idx = md.mesh.edges[:,0]-1
   v2_idx = md.mesh.edges[:,1]-1

   v1_x = md.mesh.x[v1_idx]
   v1_y = md.mesh.y[v1_idx]

   v2_x = md.mesh.x[v2_idx]
   v2_y = md.mesh.y[v2_idx]

   edge_lengths = np.sqrt((v1_x-v2_x)**2 + (v1_y-v2_y)**2)
   mean_edge_length = np.mean(edge_lengths)
   # print(f'mean edge length: {int(np.ceil(mean_edge_length))} m')
   # print(f'number of edges: {v1_x.size}')
   # print(f'number of nodes: {md.mesh.x.size}')
   # print('\n')

   nodeCoords_lon = np.append(nodeCoords_lon,md.mesh.long)
   nodeCoords_lat = np.append(nodeCoords_lat,md.mesh.lat)

   # shift nodeIds by total number of nodes from previous models	
   elementConn_n1 = np.append(elementConn_n1,md.mesh.elements[:,0] - 1 + np.sum(num_nodes)  )
   elementConn_n2 = np.append(elementConn_n2,md.mesh.elements[:,1] - 1 + np.sum(num_nodes)  )
   elementConn_n3 = np.append(elementConn_n3,md.mesh.elements[:,2] - 1 + np.sum(num_nodes)  )
	
   num_edges = np.append(num_edges,[v1_x.size])
   num_nodes = np.append(num_nodes,[md.mesh.x.size])
   mean_edges = np.append(mean_edges,[mean_edge_length])
   i += 1

nodeCoords = np.column_stack((nodeCoords_lon, nodeCoords_lat))
elementConn = np.column_stack((elementConn_n1, elementConn_n2,elementConn_n3))

#print(nodeCoords.shape)

global_mean = 0
total_nodes = int(np.sum(num_nodes))
total_edges = int(np.sum(num_edges))
# print(f'total nodes: {total_nodes}')
for j in range(np.size(num_nodes)):
   global_mean += num_edges[j]*mean_edges[j]/total_edges

global_mean = int(np.round(global_mean,0))

#print(f'mean edge length method: {global_mean} m')

#print('\n')
#print('============================================================')
#print(f'Domain name:)
print(f'ISSM_ME{global_mean}_N{total_nodes}_{"_".join(model_names)}')
#print('============================================================')


N = nodeCoords.shape[0]
M = elementConn.shape[0]

with Dataset("ISSM_MESH.nc", "w", format="NETCDF4") as nc:

	# --- Dimensions ---
	nc.createDimension("nNodes", N)
	nc.createDimension("nElements", M)
	nc.createDimension("nVertices", 3)  # triangles

	# --- Mesh topology variable (UGRID convention) ---
	mesh = nc.createVariable("mesh", "i4")
	mesh.cf_role = "mesh_topology"
	mesh.topology_dimension = 2
	mesh.node_coordinates = "node_lon node_lat"
	mesh.face_node_connectivity = "element_conn"

	# --- Node coordinates ---
	node_lon = nc.createVariable("node_lon", "f8", ("nNodes",))
	node_lon.standard_name = "longitude"
	node_lon.units = "degrees_east"
	node_lon[:] = nodeCoords[:, 0]
	
	node_lat = nc.createVariable("node_lat", "f8", ("nNodes",))
	node_lat.standard_name = "latitude"
	node_lat.units = "degrees_north"
	node_lat[:] = nodeCoords[:, 1]

	# --- Element connectivity ---
	conn = nc.createVariable("element_conn", "i4", ("nElements", "nVertices"))
	conn.cf_role = "face_node_connectivity"
	conn.start_index = 0  # 0-based indexing
	conn[:] = elementConn

	# --- Global attributes ---
	nc.Conventions = "UGRID-1.0"

