generate_issm_bcs.sh creates ISSM input files (ISSM*.bin and ISSM*.toolkits) to be used in GEOS. 
The "ISSM" prefix is used to make sure ISSM doesn't inadvertently try to read other binary files.   
Run via: sbatch generate_issm_bcs.sh 
Data is currently read from: /discover/nobackup/agstubbl/ISSM/data/
That directory is organized into AIS (Antarctica) and GRIS (Greenland) subdirectories.

The script finds all subdirectories containing files of the form:
glaciername/glaciername_meshgen.py
glaciername/glaciername_parameterize.py
glaciername/glaciername_control.py
glaciername/glaciername_finalize.py

where "glaciername" (e.g., AIS, GRIS, etc...) is the name of the subdirectory. 

The subdirectories can be nested or organized in any way (i.e., directories without the required 
python files are passed over), so future development could add a structure like:

iceland/vatnajokull/vatnajokull*.py
iceland/snaefellsjokull/snaefellsjokull*.py
...

Upon running the required python files, two files will be produced for each glacier found: 
ISSM_glaciername.bin and ISSM_glaciername.toolkits

The bin file contains the mesh, boundary conditions, physical parameters, etc., while the toolkits
file contains configuration options for external packages (usually just PETSc solver options).

The script then calls utils_issm/domain_name.py, which calculates the mean length of every triangle
edge across all glaciers (i.e. mean node spacing) and calculates the total number of nodes.

A domain name is then prescribed as:

ISSM_ME{ mean edge length }_N{ total nodes }_{ top-level glacier names separated by _ }

*The mean edge length (in meters) is rounded to the nearest 1000m.

*Top-level glacier names are those that exist at the issm directory level. In the iceland example
 above, iceland would appear in domain_name while vatnajokull and snaefellsjokull would not. 

Finally, the script creates a directory named domain_name and copies all ISSM*.bin and 
ISSM*toolkits files there.



