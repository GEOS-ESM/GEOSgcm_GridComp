from geos_state_bias.processor import Processor
import xarray as xr
import numpy as np
def main():
    # Load sample data and call the predict function
    # This part can be replaced with any data source
    # For example, you can load data from a database, csv file, etc.
    # You can also pass the data as an argument to the Processor class
    # Each array should have shape of (lev, lat, lon), 
    #  excpet for PS which has shape of (lat, lon)
    ## Sample data load from a netcdf file
    ds = xr.open_dataset("/discover/nobackup/jli30/data/geos_prog/REPLAY_M2-100KM-L181-AMIP-GFDL.geosgcm_prog.20021230_1500z.nc4")
    variables = ['U', 'V', 'T', 'QV', 'QI', 'QL', 'QG', 'QR', 'QS', 'PS']
    print(ds["U"].squeeze().shape, flush=True)
    print(ds["PS"].squeeze().shape, flush=True)
    arrays = [ds[var].to_numpy().squeeze() for var in variables]

    ## Another argument needed is the path to the checkpoint directory
    # ckpt_root_path = "/path/to/checkpoint/directory"
    # If the path is static, you can hardcode it in the Processor class
    #ckpt_root_path = "/discover/nobackup/jli30/data/geos_prog/ckpt" 
    ckpt_root_path = "/discover/nobackup/jli30/geos_state/SmaAt-UNet/geos_state_bias/checkpoints/batch_2"
    ############################################################    
    
    # Processor class is responsible for making predictions
    processor = Processor(ckpt_root_path, *arrays)
    outputs = processor.predict()
    # np.savez("sample_out.npz", *outputs)

if __name__ == "__main__":
    main()
