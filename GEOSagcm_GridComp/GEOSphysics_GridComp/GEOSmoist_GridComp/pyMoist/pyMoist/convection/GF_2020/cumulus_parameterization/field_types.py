from ndsl.dsl.gt4py import IJ, IJK, Field, K
from ndsl.dsl.typing import Float, Int
from pyMoist.constants import NUMBER_OF_TRACERS
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    MAXENS1,
    MAXENS2,
    MAXENS3,
    NUMBER_OF_PLUMES,
)


# NOTE must cast to int because numpy types are not acceptable for data dimensions
# plume field types
IntFieldIJ_Plume = Field[IJ, (Int, int(NUMBER_OF_PLUMES))]
FloatFieldIJ_Plume = Field[IJ, (Float, int(NUMBER_OF_PLUMES))]
FloatField_Plume = Field[IJK, (Float, int(NUMBER_OF_PLUMES))]

# internal ensemble field types
FloatFieldIJ_ensemble_1 = Field[IJ, (Float, int(MAXENS1))]
FloatFieldIJ_ensemble_2 = Field[IJ, (Float, int(MAXENS2))]
FloatFieldIJ_ensemble_3 = Field[IJ, (Float, int(MAXENS3))]
FloatFieldIJ_Ensemble = Field[IJ, (Float, (int(MAXENS1 * MAXENS2 * MAXENS3)))]

# NOTE THESE SHOULD BE MIGRATED TO PYMOIST GLOBAL FIELD TYPES
# convection tracer field types
FloatFieldIJ_ConvectionTracers = Field[IJ, (Float, int(NUMBER_OF_TRACERS))]
FloatField_ConvectionTracers = Field[IJK, (Float, int(NUMBER_OF_TRACERS))]
FloatField_ConvectionTracers_Plume = Field[IJK, (Float, (int(NUMBER_OF_PLUMES), int(NUMBER_OF_TRACERS)))]
