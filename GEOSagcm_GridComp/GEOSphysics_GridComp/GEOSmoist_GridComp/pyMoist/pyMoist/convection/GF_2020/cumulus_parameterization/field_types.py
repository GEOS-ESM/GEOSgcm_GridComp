from ndsl.dsl.gt4py import IJ, IJK, Field
from ndsl.dsl.typing import Float, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3

IntFieldIJ_Plume = Field[IJ, (Int, (3))]
FloatFieldIJ_Plume = Field[IJ, (Float, (3))]
FloatField_Plume = Field[IJK, (Float, (3))]

FloatFieldIJ_Ensemble = Field[
    IJ, (Float, (int(MAXENS1 * MAXENS2 * MAXENS3)))
]  # cast to int because numpy types are not acceptable for data dimensions
