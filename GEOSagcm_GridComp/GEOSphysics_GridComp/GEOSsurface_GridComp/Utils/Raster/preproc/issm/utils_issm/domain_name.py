import devpath
import numpy as np
import os
from pathlib import Path
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

from model import *
from loadmodel import loadmodel

num_nodes = np.array([ ])
mean_edges = np.array([ ])

models = [
    p.name
    for p in Path("..").iterdir()
    if p.is_dir() and p.name != "utils_issm"
]

print(f'available glaciers: {models}')

for model in models:
   md = loadmodel(f'./{model}/netcdfs/{model}_mesh.nc')

   v1_idx = md.mesh.edges[:,0]-1
   v2_idx = md.mesh.edges[:,1]-1

   v1_x = md.mesh.x[v1_idx]
   v1_y = md.mesh.y[v1_idx]

   v2_x = md.mesh.x[v2_idx]
   v2_y = md.mesh.y[v2_idx]

   edge_lengths = np.sqrt((v1_x-v2_x)**2 + (v1_y-v2_y)**2)
   mean_edge_length = np.mean(edge_lengths)
   print(f'mean edge length: {int(np.ceil(mean_edge_length))} m')
   print(f'number of edges: {v1_x.size}')
   print(f'number of nodes: {md.mesh.x.size}')
   print('\n')
   
   num_nodes = np.append(num_nodes,[md.mesh.x.size])
   mean_edges = np.append(mean_edges,[mean_edge_length])

global_mean = 0
total_nodes = int(np.sum(num_nodes))
print(f'total nodes: {total_nodes}')
for j in range(np.size(num_nodes)):
   global_mean += num_nodes[j]*mean_edges[j]/total_nodes

global_mean = int(np.round(global_mean/1000,0)*1000)

print(f'mean edge length: {global_mean} m')

print('\n')
print('============================================================')
print(f'Domain name: ISSM_ME{global_mean}_N{total_nodes}_AIS_GRIS')
print('============================================================')
