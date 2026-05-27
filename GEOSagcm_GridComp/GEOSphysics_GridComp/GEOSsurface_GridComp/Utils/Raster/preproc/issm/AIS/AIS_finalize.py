import devpath
import os
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

from model import *
from loadmodel import loadmodel
from clusters.discover_geos import export_discover
from loadresultsfromdisk import loadresultsfromdisk
from marshall import marshall
from verbose import verbose 

md = loadmodel('./netcdfs/AIS_inversion.nc')
md = loadresultsfromdisk(md,'ISSM_AIS.outbin')
md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient


# Write the binary input file
# Additional options
md.inversion.iscontrol = 0
md.transient.requested_outputs = ['default']
md.transient.isthermal=0
md.settings.waitonlock = 0
md.private.solution = 'Transient'
md.verbose = verbose('000000000')
md.toolkits = toolkits()
marshall(md,'ISSM_'+md.miscellaneous.name+'.bin') # create .bin file
md.toolkits.ToolkitsFile('ISSM_'+md.miscellaneous.name + '.toolkits')
export_discover(md,'./netcdfs/AIS_initialization.nc',delete_rundir=True)

