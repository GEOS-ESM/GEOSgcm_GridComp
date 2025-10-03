from ndsl import StencilFactory
from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, min, abs, floor
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.convective_parameterization.GF_2020.config import GF2020Config
import pyMoist.constants as constants
from pyMoist.convective_parameterization.GF_2020.temporaries import Temporaries
from pyMoist.convective_parameterization.GF_2020.state import MixingRatios
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.types import GlobalTable_saturation_tables
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


def setup(
    p_interface: FloatField,
    p: FloatField,
    p_kappa: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    geopotential_height_interface: FloatField,
    t: FloatField,
    th: FloatField,
    vapor: FloatField,
    mass: FloatField,
    tpwi: FloatFieldIJ,
    tpwi_star: FloatFieldIJ,
    seed_convection: FloatFieldIJ,
    area: FloatFieldIJ,
    modified_area: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    """
    Performs initial setup for the GF 2020 convection scheme:
     - Compute derived states
     - Initalize stochastic variability for convection
     - Modify area (m^2) here so GF scale dependence has a convection_fraction dependence

    This stencil MUST be built using Z_INTERFACE_DIM to function properly.
    """
    from __externals__ import k_end, STOCHASTIC_CONVECTION, STOCH_TOP, STOCH_BOT, GF_MIN_AREA

    # compute derived states
    with computation(PARALLEL), interval(...):
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)

    with computation(PARALLEL), interval(0, -1):
        p = 0.5 * (p_interface + p_interface[0, 0, 1])
        p_kappa = (p / constants.MAPL_P00) ** (constants.MAPL_KAPPA)
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        th = t / p_kappa
        mass = (p_interface[0, 0, 1] - p_interface) / constants.MAPL_GRAV

    with computation(FORWARD), interval(0, 1):
        tpwi = vapor * mass
        qsat, _ = saturation_specific_humidity(t, p, ese, esx)
        tpwi_star = qsat * mass

    with computation(FORWARD), interval(1, -1):
        tpwi = tpwi + vapor * mass
        qsat, _ = saturation_specific_humidity(t, p, ese, esx)
        tpwi_star = tpwi_star + qsat * mass

        # Initalize stochastic variability for convection
        if STOCHASTIC_CONVECTION == True:  # noqa
            # Create bit-processor-reproducible random white noise for convection [0:1]
            seedini = 1000000 * (100 * t.at(K=k_end) - floor(100 * t.at(K=k_end)))
            seed_convection = sqrt(max(min(seedini / 1000000.0, 1.0), 0.0))
            # Create stochastic variability to GF sigma
            seed_convection = sqrt(1.0 - (1.0 - seed_convection)) * (STOCH_TOP - STOCH_BOT) + STOCH_BOT
        else:
            seed_convection = 1.0

        # Modify area (m^2) here so GF scale dependence has a convection_fraction dependence
        if GF_MIN_AREA > 0:
            if area > GF_MIN_AREA:
                modified_area = GF_MIN_AREA * convection_fraction + area * (1.0 - convection_fraction)
            else:
                modified_area = area
        elif GF_MIN_AREA < 0:
            if area > abs(GF_MIN_AREA):
                modified_area = area * convection_fraction + abs(GF_MIN_AREA) * (1.0 - convection_fraction)
            else:
                modified_area = area
        else:
            modified_area = area


class Setup:
    """
    Performs initial setup for the GF 2020 convection scheme:
     - Compute derived states
     - Initalize stochastic variability for convection
     - Modify area (m^2) here so GF scale dependence has a convection_fraction dependence

    """

    def __init__(self, stencil_factory: StencilFactory, GF_2020_config: GF2020Config):
        self.setup = stencil_factory.from_dims_halo(
            func=setup,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            externals={
                "STOCHASTIC_CONVECTION": GF_2020_config.STOCHASTIC_CONVECTION,
                "STOCH_TOP": GF_2020_config.STOCH_TOP,
                "STOCH_BOT": GF_2020_config.STOCH_BOT,
                "GF_MIN_AREA": GF_2020_config.GF_MIN_AREA,
            },
        )

    def __call__(
        self,
        p_interface,
        geopotential_height_interface,
        t,
        mixing_ratios: MixingRatios,
        area,
        convection_fraction,
        temporaries: Temporaries,
        saturation_tables: SaturationVaporPressureTable,
    ):
        self.setup(
            p_interface,
            temporaries.p,
            temporaries.p_kappa,
            temporaries.edge_height_above_surface,
            temporaries.layer_height_above_surface,
            geopotential_height_interface,
            t,
            temporaries.th,
            mixing_ratios.vapor,
            temporaries.mass,
            temporaries.tpwi,
            temporaries.tpwi_star,
            temporaries.seed_convection,
            area,
            temporaries.modified_area,
            convection_fraction,
            saturation_tables.ese,
            saturation_tables.esx,
        )
