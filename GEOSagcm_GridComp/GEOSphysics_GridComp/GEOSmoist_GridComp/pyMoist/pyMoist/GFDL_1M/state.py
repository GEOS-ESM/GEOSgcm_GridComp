from ndsl.dsl.typing import FloatField
from dataclasses import dataclass


@dataclass
class LiquidWaterStaticEnergy:
    """
    Units:
        flux: K m s-1
        variance: K+2
        third_moment: K+3
    """

    flux: FloatField
    variance: FloatField
    third_moment: FloatField


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

    vapor: FloatField
    rain: FloatField
    snow: FloatField
    graupel: FloatField
    convective_liquid: FloatField
    convective_ice: FloatField
    large_scale_liquid: FloatField
    large_scale_ice: FloatField


@dataclass
class CloudFractions:
    """
    Units:
        convective: 1
        large_scale: 1
    """

    convective: FloatField
    large_scale: FloatField


@dataclass
class TotalWater:
    """
    Units:
        flux: kg kg-1 m s-1
        variance: 1
        third_moment: 1
    """

    flux: FloatField
    variance: FloatField
    third_moment: FloatField


@dataclass
class VericalMotion:
    """
    Units:
        velocity: m s-1
        variance: m2 s-2
        third_moment: m3 s-3
    """

    velocity: FloatField
    variance: FloatField
    third_moment: FloatField
