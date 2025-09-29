import dataclasses
from typing import Any, Self

import dacite
from numpy.typing import ArrayLike

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM, Float


@dataclasses.dataclass
class NDSLState:
    """[! Experimental !] Collection of Quantities.
    Can do zero copy with checks."""

    @classmethod
    def zeros(cls, quantity_factory: QuantityFactory) -> Self:
        """Init all quantities to zeros - included nested ones"""

        def _zeros(cls):
            initial_quantities = {}
            for _field in dataclasses.fields(cls):
                if dataclasses.is_dataclass(_field.type):
                    initial_quantities[_field.name] = _zeros(_field.type)
                else:
                    if "dims" not in _field.metadata.keys():
                        raise ValueError(
                            "Malformed state - no dims to init "
                            f"Quantity in  {_field.name} of type {_field.type}"
                        )

                    initial_quantities[_field.name] = quantity_factory.zeros(
                        _field.metadata["dims"],
                        _field.metadata["units"],
                        dtype=_field.metadata["dtype"],
                        allow_mismatch_float_precision=True,
                    )

            return initial_quantities

        dict_of_qty = _zeros(cls)
        return dacite.from_dict(data_class=cls, data=dict_of_qty)

    def init_from_memory(self, memory_map: dict[str, Any]):
        """Copy data from the input by following dataclass naming convention"""

        def _init_from_memory(dataclss, memory_map: dict[str, Any]):
            for name, array in memory_map.items():
                if isinstance(array, dict):
                    _init_from_memory(dataclss.__getattribute__(name), array)
                else:
                    try:
                        dataclss.__getattribute__(name).field[:] = array
                    except ValueError as e:
                        e.add_note(f"Error when initializing field {name} on state {type(self)}")
                        raise e

        _init_from_memory(self, memory_map)

    def init_zero_copy(self, memory_map: dict[str, Any], check: bool = True):
        """Swap buffers given into the Quantities carried by the state
        by following dataclass naming convention"""

        def _init_from_memory(dataclss, memory_map: dict[str, Any | ArrayLike]):
            for name, array in memory_map.items():
                if isinstance(array, dict):
                    _init_from_memory(dataclss.__getattribute__(name), array)
                else:
                    qty = dataclss.__getattribute__(name)
                    if check:
                        if array.shape != qty.field.shape:
                            e = ValueError("Shape mismatch on zero copy for")
                            e.add_note(f"  Error on {name} for {type(dataclss)}")
                            e.add_note(f"  Shapes: {array.shape} != {qty.field.shape}")
                            raise e
                        if array.strides != qty.data.strides:
                            e = ValueError("Stride mismatch on zero copy for")
                            e.add_note(f"  Error on {name} for {type(dataclss)}")
                            e.add_note(f"  Strides: {array.strides} != {qty.data.strides}")
                            raise e

                    qty.data = array

        _init_from_memory(self, memory_map)


@dataclasses.dataclass
class MicrophysicState(NDSLState):
    @dataclasses.dataclass
    class MixingRatios:
        """
        Mixing ratios of water species
        """

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
        convective_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "convective liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_ice: Quantity = dataclasses.field(
            metadata={
                "name": "convective ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "large scale liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_ice: Quantity = dataclasses.field(
            metadata={
                "name": "large scale ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "kg kg-1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class CloudFractions:
        """
        ?
        """

        convective: Quantity = dataclasses.field(
            metadata={
                "name": "cloud fractions convective",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale: Quantity = dataclasses.field(
            metadata={
                "name": "cloud fractions large scale",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class VerticalMotion:
        """
        ?
        """

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
                "units": " m3 s-3",
                "intent": "?",
                "dtype": Float,
            }
        )

    mixing_ratios: MixingRatios
    cloud_fractions: CloudFractions
    vertical_motion: VerticalMotion

    # @classmethod
    # def zeros(cls, quantity_factory: QuantityFactory) -> "MicrophysicState":
    #     """Init all quantities to zeros - included nested ones"""

    #     def _zeros(cls):
    #         initial_quantities = {}
    #         for _field in dataclasses.fields(cls):
    #             if dataclasses.is_dataclass(_field.type):
    #                 initial_quantities[_field.name] = _zeros(_field.type)
    #             else:
    #                 if "dims" not in _field.metadata.keys():
    #                     raise ValueError(
    #                         "Malformed state - no dims to init "
    #                         f"Quantity in  {_field.name} of type {_field.type}"
    #                     )

    #                 initial_quantities[_field.name] = quantity_factory.zeros(
    #                     _field.metadata["dims"],
    #                     _field.metadata["units"],
    #                     dtype=_field.metadata["dtype"],
    #                     allow_mismatch_float_precision=True,
    #                 )

    #         return initial_quantities

    #     dict_of_qty = _zeros(cls)
    #     return dacite.from_dict(data_class=cls, data=dict_of_qty)

    # def init_from_memory(self, memory_map: dict[str, Any]):
    #     """Will copy data from the memory map if it follows the nested
    #     naming convention of the dataclass"""

    #     def _init_from_memory(dataclss, memory_map: dict[str, Any]):
    #         for name, array in memory_map.items():
    #             if isinstance(array, dict):
    #                 _init_from_memory(dataclss.__getattribute__(name), array)
    #             else:
    #                 dataclss.__getattribute__(name).field[:] = array

    #     _init_from_memory(self, memory_map)


@dataclasses.dataclass
class Outputs(NDSLState):
    """
    Collection of all fields computed locally and returned to the rest of the model
    """

    lower_tropospheric_stability: Quantity = dataclasses.field(
        metadata={
            "name": "lower_tropospheric_stability",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    estimated_inversion_strength: Quantity = dataclasses.field(
        metadata={
            "name": "estimated_inversion_strength",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    du_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "du_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dv_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dv_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dt_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dt_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapor_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dvapor_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dliquid_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dliquid_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dice_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dice_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dcloud_fraction_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dcloud_fraction_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    drain_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "drain_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dsnow_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dsnow_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dgraupel_dt_macro: Quantity = dataclasses.field(
        metadata={
            "name": "dgraupel_dt_macro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    du_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "du_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dv_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dv_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dt_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dt_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapor_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dvapor_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dliquid_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dliquid_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dice_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dice_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dcloud_fraction_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dcloud_fraction_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    drain_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "drain_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dsnow_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dsnow_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dgraupel_dt_micro: Quantity = dataclasses.field(
        metadata={
            "name": "dgraupel_dt_micro",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_cloud_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_cloud_fraction",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_ice: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_ice",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_liquid: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_liquid",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_vapor",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_rain: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_rain",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_snow: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_snow",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    radiation_graupel: Quantity = dataclasses.field(
        metadata={
            "name": "radiation_graupel",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_radius: Quantity = dataclasses.field(
        metadata={
            "name": "ice_radius",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_radius: Quantity = dataclasses.field(
        metadata={
            "name": "liquid_radius",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_nonanvil_precipitation_evaporation: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_nonanvil_precipitation_evaporation",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_nonanvil_precipitation_sublimation: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_nonanvil_precipitation_sublimation",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    relative_humidity_after_pdf: Quantity = dataclasses.field(
        metadata={
            "name": "relative_humidity_after_pdf",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitated_rain: Quantity = dataclasses.field(
        metadata={
            "name": "precipitated_rain",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitated_snow: Quantity = dataclasses.field(
        metadata={
            "name": "precipitated_snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitated_ice: Quantity = dataclasses.field(
        metadata={
            "name": "precipitated_ice",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitated_graupel: Quantity = dataclasses.field(
        metadata={
            "name": "precipitated_graupel",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_precip: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_precip",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_snow: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    icefall: Quantity = dataclasses.field(
        metadata={
            "name": "icefall",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    freezing_rainfall: Quantity = dataclasses.field(
        metadata={
            "name": "freezing_rainfall",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_nonanvil_ice_flux: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_nonanvil_ice_flux",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_nonanvil_liquid_flux: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_nonanvil_liquid_flux",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    anvil_liquid_flux: Quantity = dataclasses.field(
        metadata={
            "name": "anvil_liquid_flux",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    anvil_ice_flux: Quantity = dataclasses.field(
        metadata={
            "name": "anvil_ice_flux",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    deep_convective_precipitation: Quantity = dataclasses.field(
        metadata={
            "name": "deep_convective_precipitation",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    anvil_precipitation: Quantity = dataclasses.field(
        metadata={
            "name": "anvil_precipitation",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    shallow_convective_precipitation: Quantity = dataclasses.field(
        metadata={
            "name": "shallow_convective_precipitation",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    deep_convective_snow: Quantity = dataclasses.field(
        metadata={
            "name": "deep_convective_snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    anvil_snow: Quantity = dataclasses.field(
        metadata={
            "name": "anvil_snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    shallow_convective_snow: Quantity = dataclasses.field(
        metadata={
            "name": "shallow_convective_snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    # Optional outputs
    z_lcl: Quantity = dataclasses.field(
        metadata={
            "name": "z_lcl",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_rainwater_source: Quantity | None = dataclasses.field(
        metadata={
            "name": "large_scale_rainwater_source",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    moist_friction_temperature_tendency: Quantity | None = dataclasses.field(
        metadata={
            "name": "moist_friction_temperature_tendency",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    simulated_reflectivity: Quantity | None = dataclasses.field(
        metadata={
            "name": "simulated_reflectivity",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    maximum_reflectivity: Quantity | None = dataclasses.field(
        metadata={
            "name": "maximum_reflectivity",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    one_km_agl_reflectivity: Quantity | None = dataclasses.field(
        metadata={
            "name": "one_km_agl_reflectivity",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    echo_top_reflectivity: Quantity | None = dataclasses.field(
        metadata={
            "name": "echo_top_reflectivity",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    minus_10c_reflectivity: Quantity | None = dataclasses.field(
        metadata={
            "name": "minus_10c_reflectivity",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )


@dataclasses.dataclass
class MicrophysicsDiagnostics:
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
                "name": "third moment",
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
                "name": "third moment",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "1",
                "intent": "?",
                "dtype": Float,
            }
        )

    liquid_water_static_energy: LiquidWaterStaticEnergy
    total_water: TotalWater
