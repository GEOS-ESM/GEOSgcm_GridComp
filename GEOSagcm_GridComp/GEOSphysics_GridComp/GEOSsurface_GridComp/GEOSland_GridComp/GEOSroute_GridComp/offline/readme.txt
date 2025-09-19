README - River Routing Model Offline Version
Last Updated: 10/28/2024
Contact: yujin.zeng@nasa.gov
Overview

This directory contains the code required to run the offline version of the river routing model. Note that not all files in this directory pertain to the offline model. Key files include:

    run: Script for building and running the model.
    ncdioMod.f90: Local NetCDF library.
    rwncMod.f90: Local NetCDF I/O library.
    interp_M36toPfaf.f90: Interpolation module.
    river_io_mod.f90: I/O module.
    res_mod.f90: Reservoir module.
    lake_mod.f90: Lake module.
    river_routing.f90: Main program.

Running the Offline Model

    Set Directory Paths
    In river_io_mod.f90, set:
        input_dir: Path for input data, e.g., /discover/nobackup/yzeng3/work/river_routing_model_offline/input/
        runoff_dir: Path for runoff data (e.g., Catchment model 2D output in M36 or M09 resolutions).
            Example for M36 resolution: /discover/nobackup/yzeng3/GEOldas_output
        output_dir: Path for output data.

    Define Start and End Dates
    In river_routing.f90, set step_start (start date) and step_end (end date) as days since January 1, 1990 (Day 1). Ensure these dates align with the runoff forcing period.

    Build and Run
    Compile and run the model using:
    ./run river_routing.f90

Output Format

The output files are in .txt format, generated daily with date information in each filename. The output variables are as follows:

    Main river discharge (Pfaf_Qr_Kv_*.txt) [m³/s]
    Main river storage (Pfaf_Wr_Kv_*.txt) [kg]
    Local stream storage (Pfaf_Ws_Kv_*.txt) [kg]
    Reservoir outflow (Pfaf_Q_res_Kv_*.txt) [m³/s] (0 for catchments without reservoirs)
    Reservoir water storage (Pfaf_Wr_res_Kv_*.txt) [kg] (0 for catchments without reservoirs)
    Lake outflow (Pfaf_Q_lake_Kv_*.txt) [m³/s] (0 for catchments without lakes)
    Lake water storage (Pfaf_Wr_lake_Kv_*.txt) [kg] (0 for catchments without lakes)

Each .txt file contains a list of 291,284 values corresponding to catchments indexed from 1 to 291,284. To convert these lists into spatial maps, use the catchment distribution map at 1-minute resolution in CatchIndex from SRTM_PfafData.nc:

    Path: /discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc

For further assistance, please contact yujin.zeng@nasa.gov.