from dataclasses import dataclass

import numpy.typing as npt

from ndsl import Quantity


@dataclass
class LiquidWaterStaticEnergy:
    """
    Units:
        flux: K m s-1
        variance: K+2
        third_moment: K+3
    """

    flux: Quantity | npt.NDArray
    variance: Quantity | npt.NDArray
    third_moment: Quantity | npt.NDArray


@dataclass
class MixingRatios:
    """
    Units:
        vapor: kg kg-1
        rain: kg kg-1
        snow: kg kg-1
        graupel: kg kg-1
        convective_liquid: kg kg-1
        convective_ice: kg kg-1
        large_scale_liquid: kg kg-1
        large_scale_ice: kg kg-1
    """

    vapor: Quantity | npt.NDArray
    rain: Quantity | npt.NDArray
    snow: Quantity | npt.NDArray
    graupel: Quantity | npt.NDArray
    convective_liquid: Quantity | npt.NDArray
    convective_ice: Quantity | npt.NDArray
    large_scale_liquid: Quantity | npt.NDArray
    large_scale_ice: Quantity | npt.NDArray


@dataclass
class CloudFractions:
    """
    Units:
        convective: 1
        large_scale: 1
    """

    convective: Quantity | npt.NDArray
    large_scale: Quantity | npt.NDArray


@dataclass
class TotalWater:
    """
    Optional outputs of Hydrostatic PDF for PDF_Shape=5

    Units:
        flux: kg kg-1 m s-1
        variance: 1
        third_moment: 1
    """

    flux: Quantity | npt.NDArray
    variance: Quantity | npt.NDArray
    third_moment: Quantity | npt.NDArray


@dataclass
class VerticalMotion:
    """
    Units:
        velocity: m s-1
        variance: m2 s-2
        third_moment: m3 s-3
    """

    velocity: Quantity | npt.NDArray
    variance: Quantity | npt.NDArray
    third_moment: Quantity | npt.NDArray
