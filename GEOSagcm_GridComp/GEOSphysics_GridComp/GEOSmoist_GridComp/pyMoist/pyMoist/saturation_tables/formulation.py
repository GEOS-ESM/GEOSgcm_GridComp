import enum


class SaturationFormulation(enum.Enum):
    """
    The choice of saturation vapor pressure formulation is a compile-time
    option. Three choices are currently supported: The CAM formulation,
    Murphy and Koop (2005, QJRMS), and Staars formulation from NSIPP-1.
    """

    Staars = 1
    CAM = 2
    MurphyAndKoop = 3
