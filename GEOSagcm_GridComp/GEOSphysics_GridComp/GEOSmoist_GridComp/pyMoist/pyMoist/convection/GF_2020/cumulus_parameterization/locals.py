import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020CumulusParameterizationLocals(State):
    t_new: Quantity = dataclasses.field(
        metadata={
            "name": "t_new",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_excess: Quantity = dataclasses.field(
        metadata={
            "name": "t_excess",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_new: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_new",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_excess: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_excess",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_new_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "t_new_pbl",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_new_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_new_pbl",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    moist_static_energy: Quantity = dataclasses.field(
        metadata={
            "name": "moist_static_energy",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    kbmax: Quantity = dataclasses.field(
        metadata={
            "name": "kbmax",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    kstamb: Quantity = dataclasses.field(
        metadata={
            "name": "kstamb",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ocean_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "ocean_fraction",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cap_max: Quantity = dataclasses.field(
        metadata={
            "name": "cap_max",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ierrc: Quantity = dataclasses.field(
        metadata={
            "name": "ierrc",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    ierr2: Quantity = dataclasses.field(
        metadata={
            "name": "ierr2",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    ierr3: Quantity = dataclasses.field(
        metadata={
            "name": "ierr3",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    max_increment: Quantity = dataclasses.field(
        metadata={
            "name": "max_increment",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_modified: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_0: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_0",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_2: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_2",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_3: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_3",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_0_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_0_pbl",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_pbl",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1_fa: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_fa",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cin1: Quantity = dataclasses.field(
        metadata={
            "name": "cin1",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    k_x_modified: Quantity = dataclasses.field(
        metadata={
            "name": "k_x_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    epsilon: Quantity = dataclasses.field(
        metadata={
            "name": "epsilon",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pbl_time_scale: Quantity = dataclasses.field(
        metadata={
            "name": "pbl_time_scale",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_wetbulb: Quantity = dataclasses.field(
        metadata={
            "name": "t_wetbulb",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_wetbulb: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_wetbulb",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    tau_ecmwf: Quantity = dataclasses.field(
        metadata={
            "name": "tau_ecmwf",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    f_dicycle_modified: Quantity = dataclasses.field(
        metadata={
            "name": "f_dicycle_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    add_buoy_modified: Quantity = dataclasses.field(
        metadata={
            "name": "add_buoy_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hcdo: Quantity = dataclasses.field(
        metadata={
            "name": "hcdo",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cupclw: Quantity = dataclasses.field(
        metadata={
            "name": "cupclw",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    qrcdo: Quantity = dataclasses.field(
        metadata={
            "name": "qrcdo",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hcot: Quantity = dataclasses.field(
        metadata={
            "name": "hcot",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    c1d: Quantity = dataclasses.field(
        metadata={
            "name": "c1d",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evap_bcb: Quantity = dataclasses.field(
        metadata={
            "name": "evap_bcb",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_ensemble: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_ensemble",
            "dims": [X_DIM, Y_DIM, "ensemble_members"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitation_ensemble: Quantity = dataclasses.field(
        metadata={
            "name": "precipitation_ensemble",
            "dims": [X_DIM, Y_DIM, "ensemble_members"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    scale_dependence_factor: Quantity = dataclasses.field(
        metadata={
            "name": "scale_dependence_factor",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    scale_dependence_factor_downdraft: Quantity = dataclasses.field(
        metadata={
            "name": "scale_dependence_factor_downdraft",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    random_number: Quantity = dataclasses.field(
        metadata={
            "name": "random_number",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
