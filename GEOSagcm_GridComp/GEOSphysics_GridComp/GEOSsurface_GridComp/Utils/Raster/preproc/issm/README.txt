generate_issm_bcs.sh creates ISSM input files (extension *.bin and *toolkits) to be used in GEOS.   

The script finds all subdirectories containing files of the form:
glaciername/glaciername_meshgen.py
glaciername/glaciername_parameterize.py
glaciername/glaciername_control.py
glaciername/glaciername_finalize.py

where "glaciername" (e.g., AIS, GRIS, etc...) is the name of the subdirectory. 

The subdirectories can be nested or organized in any way (i.e., directories without the require 
python files are passed over), so future development could have a structure like:

iceland/vatnajokull/vatnajokull*.py
iceland/snaefellsjokull/snaefellsjokull*.py
...

Upon running the required python files, two files will be produced for each glacier found: 
glaciername.bin and glaciername.toolkits.

The script then calls utils_issm/domain_name.py, which calculates the mean length of every triangle
edge across all glaciers and calculates the total number of nodes across all glaciers. 

A domain name is then prescribed as:

ISSM_ME{ mean edge length }_N{ total nodes }_{ top-level glacier names separated by _ }

Top-level glacier names are those that exist at the issm directory level. So, in the iceland
example above, iceland would appear in domain_name while vatnajokull and snaefellsjökull would not. 

Finally, the script creates a directory named domain_name and copies all *.bin and *toolkits there.



