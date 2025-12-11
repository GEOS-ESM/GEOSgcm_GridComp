import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float


@dataclasses.dataclass
class GFDL1MState(State):
    area: Quantity = dataclasses.field(
        metadata={
            "name": "area",
            "dims": [X_DIM, Y_DIM],
            "units": "m2",
            "intent": "?",
            "dtype": Float,
        }
    )
    z_interface: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_interface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_interface: Quantity = dataclasses.field(
        metadata={
            "name": "p_interface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "Pa",
            "intent": "?",
            "dtype": Float,
        }
    )
    t: Quantity = dataclasses.field(
        metadata={
            "name": "t",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "K",
            "intent": "?",
            "dtype": Float,
        }
    )
    u: Quantity = dataclasses.field(
        metadata={
            "name": "u",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    v: Quantity = dataclasses.field(
        metadata={
            "name": "v",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    land_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "land_fraction",
            "dims": [X_DIM, Y_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    scalar_diffusivity_interface: Quantity = dataclasses.field(
        metadata={
            "name": "scalar_diffusivity_interface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "m2 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    pdf_first_plume_fractional_area: Quantity = dataclasses.field(
        metadata={
            "name": "pdf_first_plume_fractional_area",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    covariance_liquid_water_static_energy_and_total_water_specific_humidity: Quantity = dataclasses.field(
        metadata={
            "name": "covariance_liquid_water_static_energy_and_total_water_specific_humudity",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "K",
            "intent": "?",
            "dtype": Float,
        }
    )
    surface_temperature: Quantity = dataclasses.field(
        metadata={
            "name": "surface_temperature",
            "dims": [X_DIM, Y_DIM],
            "units": "K",
            "intent": "?",
            "dtype": Float,
        }
    )
    sensible_heat_flux: Quantity = dataclasses.field(
        metadata={
            "name": "sensible_heat_flux",
            "dims": [X_DIM, Y_DIM],
            "units": "W m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    omega: Quantity = dataclasses.field(
        metadata={
            "name": "omega",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "Pa s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "convection_fraction",
            "dims": [X_DIM, Y_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    surface_type: Quantity = dataclasses.field(
        metadata={
            "name": "surface_type",
            "dims": [X_DIM, Y_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_liquid_evaporation: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_liquid_evaporation",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_ice_sublimation: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_ice_sublimation",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    icefall: Quantity = dataclasses.field(
        metadata={
            "name": "icefall",
            "dims": [X_DIM, Y_DIM],
            "units": "kg m-2 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    freezing_rainfall: Quantity = dataclasses.field(
        metadata={
            "name": "freezing_rainfall",
            "dims": [X_DIM, Y_DIM],
            "units": "kg m-2 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    relative_humidity_after_pdf: Quantity = dataclasses.field(
        metadata={
            "name": "relative_humidity_after_pdf",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    buoyancy_flux: Quantity = dataclasses.field(
        metadata={
            "name": "buoyancy_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_water_flux: Quantity = dataclasses.field(
        metadata={
            "name": "liquid_water_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    hydrostatic_pdf_iterations: Quantity = dataclasses.field(
        metadata={
            "name": "hydrostatic_pdf_iterations",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    lower_tropospheric_stability: Quantity = dataclasses.field(
        metadata={
            "name": "lower_tropospheric_stability",
            "dims": [X_DIM, Y_DIM],
            "units": "K",
            "intent": "?",
            "dtype": Float,
        }
    )
    estimated_inversion_strength: Quantity = dataclasses.field(
        metadata={
            "name": "estimated_inversion_strength",
            "dims": [X_DIM, Y_DIM],
            "units": "K",
            "intent": "?",
            "dtype": Float,
        }
    )
    lcl_height: Quantity = dataclasses.field(
        metadata={
            "name": "lcl_height",
            "dims": [X_DIM, Y_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    shallow_convection_rain: Quantity = dataclasses.field(
        metadata={
            "name": "shallow_convection_rain",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    shallow_convection_snow: Quantity = dataclasses.field(
        metadata={
            "name": "shallow_convective_snow",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    critical_relative_humidity_for_pdf: Quantity = dataclasses.field(
        metadata={
            "name": "critical_relative_humidity_for_pdf",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "1",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_rainwater_source: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_rainwater_source",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg kg-1 s-1",
            "intent": "?",
            "dtype": Float,
        }
    )

    @dataclasses.dataclass
    class VerticalMotion:
        velocity: Quantity = dataclasses.field(
            metadata={
                "name": "velocity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        variance: Quantity = dataclasses.field(
            metadata={
                "name": "variance",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m2 s-2",
                "intent": "?",
                "dtype": Float,
            }
        )
        third_moment: Quantity = dataclasses.field(
            metadata={
                "name": "third_moment",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m3 s-3",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class MixingRatio:
        vapor: Quantity = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        rain: Quantity = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Quantity = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_ice: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "convective_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_ice: Quantity = dataclasses.field(
            metadata={
                "name": "convective_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class CloudFraction:
        large_scale: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective: Quantity = dataclasses.field(
            metadata={
                "name": "convective_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Concentration:
        liquid: Quantity = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m-3",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m-3",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class LiquidWaterStaticEnergy:
        """
        Units:
            flux: K m s-1
            variance: K+2
            third_moment: K+3
        """

        flux: Quantity = dataclasses.field(
            metadata={
                "name": "flux",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "K m s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        variance: Quantity = dataclasses.field(
            metadata={
                "name": "variance",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "K+2",
                "intent": "?",
                "dtype": Float,
            }
        )
        third_moment: Quantity = dataclasses.field(
            metadata={
                "name": "third_moment",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "K+3",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class TotalWater:
        """
        Optional outputs of Hydrostatic PDF for PDF_Shape=5

        Units:
            flux: kg kg-1 m s-1
            variance: 1
            third_moment: 1
        """

        flux: Quantity = dataclasses.field(
            metadata={
                "name": "flux",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 m s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        variance: Quantity = dataclasses.field(
            metadata={
                "name": "variance",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )
        third_moment: Quantity = dataclasses.field(
            metadata={
                "name": "third_moment",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class RadiationField:
        cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor: Quantity = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        liquid: Quantity = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        rain: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Quantity = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class CloudParticleEffectiveRadius:
        liquid: Quantity = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "m",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class PrecipitationAtSurface:
        rain: Quantity = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Quantity = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        shallow_convective_precipitation: Quantity = dataclasses.field(
            metadata={
                "name": "shallow_convective_precipitation",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        deep_convective_precipitation: Quantity = dataclasses.field(
            metadata={
                "name": "deep_convective_precipitation",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        anvil_precipitation: Quantity = dataclasses.field(
            metadata={
                "name": "anvil_precipitation",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        shallow_convective_snow: Quantity = dataclasses.field(
            metadata={
                "name": "shallow_convective_snow",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        deep_convective_snow: Quantity = dataclasses.field(
            metadata={
                "name": "deep_convective_snow",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        anvil_snow: Quantity = dataclasses.field(
            metadata={
                "name": "anvil_snow",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class NonAnvilLargeScale:
        precip: Quantity = dataclasses.field(
            metadata={
                "name": "precip",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        evaporation: Quantity = dataclasses.field(
            metadata={
                "name": "evaporation",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        sublimation: Quantity = dataclasses.field(
            metadata={
                "name": "sublimation",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        liquid_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "liquid_precip_flux",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "ice_precip_flux",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Anvil:
        liquid_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "liquid_precip_flux",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "ice_precip_flux",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "kg m-2 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Tendencies:
        dcloud_fractiondt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dsurface_specific_humuditydt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dicedt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dicedt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dliquiddt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dliquiddt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        draindt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "draindt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dgraupeldt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dgraupeldt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dsnowdt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dsnowdt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dudt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt_macro: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_macro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloud_fractiondt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dsurface_specific_humuditydt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dicedt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dicedt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dliquiddt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dliquiddt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        draindt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "draindt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dgraupeldt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dgraupeldt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dsnowdt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dsnowdt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dudt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt_micro: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_micro",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1 s-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt_friction_pressure_weighted: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_friction_pressure_weighted",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "Pa K s-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Radar:
        simulated_reflectivity: Quantity = dataclasses.field(
            metadata={
                "name": "simulated_reflectivity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "dBZ",
                "intent": "?",
                "dtype": Float,
            }
        )
        maximum_composite_reflectivity: Quantity = dataclasses.field(
            metadata={
                "name": "maximum_composite_reflectivity",
                "dims": [X_DIM, Y_DIM],
                "units": "dBZ",
                "intent": "?",
                "dtype": Float,
            }
        )
        base_1km_agl_reflectivity: Quantity = dataclasses.field(
            metadata={
                "name": "base_1km_agl_reflectivity",
                "dims": [X_DIM, Y_DIM],
                "units": "dBZ",
                "intent": "?",
                "dtype": Float,
            }
        )
        echo_top_reflectivity: Quantity = dataclasses.field(
            metadata={
                "name": "echo_top_reflectivity",
                "dims": [X_DIM, Y_DIM],
                "units": "dBZ",
                "intent": "?",
                "dtype": Float,
            }
        )
        minus_10c_reflectivity: Quantity = dataclasses.field(
            metadata={
                "name": "minus_10c_reflectivity",
                "dims": [X_DIM, Y_DIM],
                "units": "dBZ",
                "intent": "?",
                "dtype": Float,
            }
        )

    vertical_motion: VerticalMotion
    mixing_ratio: MixingRatio
    cloud_fraction: CloudFraction
    concentration: Concentration
    liquid_water_static_energy: LiquidWaterStaticEnergy
    total_water: TotalWater
    radiation_field: RadiationField
    cloud_particle_effective_radius: CloudParticleEffectiveRadius
    precipitation_at_surface: PrecipitationAtSurface
    non_anvil_large_scale: NonAnvilLargeScale
    anvil: Anvil
    tendencies: Tendencies
    radar: Radar
