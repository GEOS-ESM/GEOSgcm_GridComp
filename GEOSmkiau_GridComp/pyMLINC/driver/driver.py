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
    levs= [
        2, 5, 8,
        11, 14, 17,
        20, 23, 26, 29,
        32, 35, 38,
        41, 44, 47,
        50, 53, 56,
        60, 64, 69,
        74,
        80, 85,
        90, 95,
        100, 105,
        110, 115,
        120, 125,
        130, 135,
        140, 145,
        152,
        160,
        170,
        181,
    ]
    levs = [x-1 for x in levs] # Fortran->C indexing
    print(ds["U"].squeeze().to_numpy()[levs, :, :].shape, flush=True)
    arrays = []
    for var in variables[:-1]: # excluding PS
        arrays.append(ds[var].squeeze().to_numpy()[levs, :, :])
    arrays.append(ds["PS"].squeeze().to_numpy())

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
