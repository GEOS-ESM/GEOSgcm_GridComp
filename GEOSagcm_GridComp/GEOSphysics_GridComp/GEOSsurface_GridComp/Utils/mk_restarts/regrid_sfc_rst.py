#!/usr/bin/env python3
'''
This program needs to run on compute nodes with at leat 56 cpus.
This program regrids surface restarts.

Usage
Change to a working dir with write permissions and create the following symlinks:
intilefile.data: input tile file
outtilefile:     output tile file
clsm:            output clsm directory (usually in same dir as output tile file)
inrst:           input restart directory
bindir:          directory with mk_LakeLaniceSaltRestarts and mk_CatchRestarts programs

Run on compute nodes as usual. Command line args interface may be added later.
'''

import os, sys

intilefile=os.path.realpath('intile.data')
outtilefile=os.path.realpath('outtile.data')
clsm=os.path.realpath('clsm')
inrstdir=os.path.realpath('inrst')
outdata='OutData'
outrstdir=outdata
bindir=os.path.realpath('bin')

zoom=8     # zoom value to send to land regridding codes [8]
surflay=50 # number of surface layers (catch & catchcn)
wemin=26   # minimum snow water equivalent threshold for input catch/cn [$weminDFLT]
wemout=13  # weminOUT minimum snow water equivalent threshold for output catch/cn [$weminDFLT]

# Create output directory (hardcoded in mk_* utils)
os.makedirs(outdata, exist_ok=True)

maskdict={'openwater_internal_rst': 0,
          'seaicethermo_internal_rst': 0,
          'lake_internal_rst': 19,
          'landice_internal_rst': 20}

rstlist=os.listdir('inrst')
cmdlist=[]
for rst, mask in maskdict.items():
    if rst in rstlist: 
        cmdlist.append(f'{bindir}/mk_LakeLandiceSaltRestarts {outtilefile} {intilefile} {inrstdir}/{rst} {mask} {zoom}')
    
# Regrid lake, landice, saltwater restarts
for cmd in cmdlist:
    print(cmd)
    os.system(cmd)

# Regrid land restarts
if 'catch_internal_rst' in rstlist:
    # Step 1
    # remove clsm if exists
    if os.path.exists(f'{outdata}/clsm'): os.remove(f'{outdata}/clsm')
    cmd=f'mpirun -np 56 {bindir}/mk_CatchRestarts {outtilefile} {intilefile} {inrstdir}/catch_internal_rst {surflay}'
    print(cmd)
    os.system(cmd)
    
    # Step 2
    # need to move catch restart to other dir
    os.makedirs('tmp',exist_ok=True)
    os.replace(f'{outdata}/catch_internal_rst',f'tmp/catch_internal_rst')
    # link clsm
    os.symlink(clsm,f'{outdata}/clsm')
    cmd=f'mpirun -np 56 {bindir}/mk_CatchRestarts {outtilefile} {outtilefile} tmp/catch_internal_rst {surflay}'
    print(cmd)
    os.system(cmd)
    
    # Step 3 scale catch restart
    os.replace(f'{outdata}/catch_internal_rst',f'tmp/catch_internal_rst1')
    cmd=f'{bindir}/Scale_Catch tmp/catch_internal_rst tmp/catch_internal_rst1 {outdata}/catch_internal_rst {surflay} {wemin} {wemout}'
    print(cmd)
    os.system(cmd)
    
    # cleanup
    [os.remove(f'tmp/{file}') for file in os.listdir('tmp')]; os.rmdir('tmp')

print('All done')
