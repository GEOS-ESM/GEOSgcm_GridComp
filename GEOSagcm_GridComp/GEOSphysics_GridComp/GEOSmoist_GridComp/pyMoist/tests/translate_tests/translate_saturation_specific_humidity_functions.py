from f90nml import Namelist
from ndsl import Backend, StencilFactory
from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, computation, interval
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    get_saturation_vapor_pressure_table,
    saturation_specific_humidity,
    saturation_specific_humidity_frozen_surface,
    saturation_specific_humidity_liquid_surface,
)


meshgrid_domain = [100, 100, 1]
nhalo = 0

meshgrid_stencil_factory, meshgrid_quantity_factory = get_factories_single_tile(
    meshgrid_domain[0], meshgrid_domain[1], meshgrid_domain[2], 0, backend=Backend("st:dace:cpu:IJK")
)


def test_saturation_specific_humidity_functions(
    t: FloatField,
    p: FloatField,
    sat_over_ice: FloatField,
    dqsat_over_ice: FloatField,
    sat_over_liquid: FloatField,
    dqsat_over_liquid: FloatField,
    sat: FloatField,
    dqsat: FloatField,
    ese: GlobalTable_saturation_tables,
    esw: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
    frz: Float,
    lqu: Float,
):
    with computation(FORWARD), interval(...):
        sat_over_ice, dqsat_over_ice = saturation_specific_humidity_frozen_surface(ese=ese, frz=frz, t=t, p=p * 100.0)
        sat_over_liquid, dqsat_over_liquid = saturation_specific_humidity_liquid_surface(esw=esw, lqu=lqu, t=t, p=p * 100.0)
        sat, dqsat = saturation_specific_humidity(t=t, p=p * 100.0, ese=ese, esx=esx)


def test_saturation_specific_humidity_functions_2d(
    t: FloatFieldIJ,
    p: FloatFieldIJ,
    sat_over_ice: FloatFieldIJ,
    dqsat_over_ice: FloatFieldIJ,
    sat_over_liquid: FloatFieldIJ,
    dqsat_over_liquid: FloatFieldIJ,
    sat: FloatFieldIJ,
    dqsat: FloatFieldIJ,
    ese: GlobalTable_saturation_tables,
    esw: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
    frz: Float,
    lqu: Float,
):
    with computation(FORWARD), interval(0, 1):
        sat_over_ice, dqsat_over_ice = saturation_specific_humidity_frozen_surface(ese=ese, frz=frz, t=t, p=p * 100.0)
        sat_over_liquid, dqsat_over_liquid = saturation_specific_humidity_liquid_surface(esw=esw, lqu=lqu, t=t, p=p * 100.0)
        sat, dqsat = saturation_specific_humidity(t=t, p=p * 100.0, ese=ese, esx=esx)


class Translatesaturation_specific_humidity_functions(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        self.data_saturation_specific_humidity_functions = self.stencil_factory.from_dims_halo(
            func=test_saturation_specific_humidity_functions,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self.meshgrid_saturation_specific_humidity_functions = meshgrid_stencil_factory.from_dims_halo(
            func=test_saturation_specific_humidity_functions_2d,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PRES_ARRAY": {},
            "TEMP_ARRAY": {},
            "PLmb": {},
            "T": {},
        }

        # FloatField Outputs
        self.out_vars: dict = {
            # regular data fields
            "SER_QSATICE": {},
            "SER_DQSI": {},
            "SER_QSATLQU": {},
            "SER_DQSL": {},
            "SER_QSAT_QS": {},
            "SER_QSAT_DQ": {},
            "SER_DQSAT_DQ": {},
            "SER_DQSAT_QS": {},
            # meshgrid size fields
            "SER_QSATICE_2D": {},  # meshgrid_stencil_factory.grid_indexing.domain_compute(),
            "SER_DQSI_2D": {},
            "SER_QSATLQU_2D": {},  # meshgrid_stencil_factory.grid_indexing.domain_compute(),
            "SER_DQSL_2D": {},
            "SER_QSAT0_QS_2D": {},  # meshgrid_stencil_factory.grid_indexing.domain_compute(),
            "SER_QSAT0_DQ_2D": {},
            "SER_DQSAT0_DQ_2D": {},  # meshgrid_stencil_factory.grid_indexing.domain_compute(),
            "SER_DQSAT0_QS_2D": {},
        }

    def compute(self, inputs):

        # Inputs
        t = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        p_mb = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        safe_assign_array(t.field[:, :, :], inputs["T"])
        safe_assign_array(p_mb.field[:, :, :], inputs["PLmb"])

        temp_grid = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        pres_grid = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        safe_assign_array(temp_grid.field[:, :], inputs["TEMP_ARRAY"])
        safe_assign_array(pres_grid.field[:, :], inputs["PRES_ARRAY"])

        # Outputs
        data_sat_over_ice = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        data_dqsat_over_ice = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        data_sat_over_liquid = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        data_dqsat_over_liquid = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        data_sat = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")
        data_dqsat = self.quantity_factory.zeros([I_DIM, J_DIM, K_DIM], "n/a")

        meshgrid_sat_over_ice = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        meshgrid_dqsat_over_ice = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        meshgrid_sat_over_liquid = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        meshgrid_dqsat_over_liquid = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        meshgrid_sat = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")
        meshgrid_dqsat = meshgrid_quantity_factory.zeros([I_DIM, J_DIM], "n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)

        self.data_saturation_specific_humidity_functions(
            t=t,
            p=p_mb,
            sat_over_ice=data_sat_over_ice,
            dqsat_over_ice=data_dqsat_over_ice,
            sat_over_liquid=data_sat_over_liquid,
            dqsat_over_liquid=data_dqsat_over_liquid,
            sat=data_sat,
            dqsat=data_dqsat,
            ese=saturation_vapor_pressure_table.ese,
            esw=saturation_vapor_pressure_table.esw,
            esx=saturation_vapor_pressure_table.esx,
            frz=saturation_vapor_pressure_table.frz,
            lqu=saturation_vapor_pressure_table.lqu,
        )

        self.meshgrid_saturation_specific_humidity_functions(
            t=temp_grid,
            p=pres_grid,
            sat_over_ice=meshgrid_sat_over_ice,
            dqsat_over_ice=meshgrid_dqsat_over_ice,
            sat_over_liquid=meshgrid_sat_over_liquid,
            dqsat_over_liquid=meshgrid_dqsat_over_liquid,
            sat=meshgrid_sat,
            dqsat=meshgrid_dqsat,
            ese=saturation_vapor_pressure_table.ese,
            esw=saturation_vapor_pressure_table.esw,
            esx=saturation_vapor_pressure_table.esx,
            frz=saturation_vapor_pressure_table.frz,
            lqu=saturation_vapor_pressure_table.lqu,
        )

        return {
            # regular data fields
            "SER_QSATICE": data_sat_over_ice.field[:],
            "SER_DQSI": data_dqsat_over_ice.field[:],
            "SER_QSATLQU": data_sat_over_liquid.field[:],
            "SER_DQSL": data_dqsat_over_liquid.field[:],
            "SER_QSAT_QS": data_sat.field[:],
            "SER_QSAT_DQ": data_dqsat.field[:],
            "SER_DQSAT_DQ": data_dqsat.field[:],
            "SER_DQSAT_QS": data_sat.field[:],
            # meshgrid fields
            "SER_QSATICE_2D": meshgrid_sat_over_ice.field[:],
            "SER_DQSI_2D": meshgrid_dqsat_over_ice.field[:],
            "SER_QSATLQU_2D": meshgrid_sat_over_liquid.field[:],
            "SER_DQSL_2D": meshgrid_dqsat_over_liquid.field[:],
            "SER_QSAT0_QS_2D": meshgrid_sat.field[:],
            "SER_QSAT0_DQ_2D": meshgrid_dqsat.field[:],
            "SER_DQSAT0_DQ_2D": meshgrid_dqsat.field[:],
            "SER_DQSAT0_QS_2D": meshgrid_sat.field[:],
        }
