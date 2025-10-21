from ndsl.dsl.gt4py import IJ, IJK, Field
from ndsl.dsl.typing import Float
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3


FloatFieldIJ_Plumes = Field[IJ, (Float, (3))]
FloatField_Plumes = Field[IJK, (Float, (3))]
FloatField_Ensemble = Field[IJK, (Float, (MAXENS1 * MAXENS2 * MAXENS3))]
