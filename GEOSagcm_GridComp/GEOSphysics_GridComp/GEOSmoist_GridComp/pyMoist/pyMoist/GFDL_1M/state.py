from ndsl import Quantity
from dataclasses import dataclass


@dataclass
class LiquidWaterStaticEnergy:
    """
    Units:
        flux: K m s-1
        variance: K+2
        third_moment: K+3
    """

    flux: Quantity
    variance: Quantity
    third_moment: Quantity


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

    vapor: Quantity
    rain: Quantity
    snow: Quantity
    graupel: Quantity
    convective_liquid: Quantity
    convective_ice: Quantity
    large_scale_liquid: Quantity
    large_scale_ice: Quantity


@dataclass
class CloudFractions:
    """
    Units:
        convective: 1
        large_scale: 1
    """

    convective: Quantity
    large_scale: Quantity


@dataclass
class TotalWater:
    """
    Units:
        flux: kg kg-1 m s-1
        variance: 1
        third_moment: 1
    """

    flux: Quantity
    variance: Quantity
    third_moment: Quantity


@dataclass
class VericalMotion:
    """
    Units:
        velocity: m s-1
        variance: m2 s-2
        third_moment: m3 s-3
    """

    velocity: Quantity
    variance: Quantity
    third_moment: Quantity
