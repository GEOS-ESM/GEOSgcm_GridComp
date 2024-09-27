"""Stencils and functions called by multiple pyMoist modules.
These functions manupilate data to formats more friendly to
gt4py and/or deal with gt4py I/O shortcomings."""

from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    computation,
    interval,
)
import pyMoist.pyMoist_constants as constants
from ndsl.dsl.typing import FloatField, FloatFieldIJ


def get_last(
    in_field: FloatField, temporary_field: FloatFieldIJ, out_field: FloatField
):
    with computation(FORWARD), interval(-1, None):
        temporary_field = in_field

    with computation(PARALLEL), interval(...):
        out_field = temporary_field


def hybrid_index_2dout(
    data_field: FloatField,
    k_mask: FloatField,
    k_index_desired: FloatFieldIJ,
    out_field: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        if k_mask == k_index_desired:
            out_field = data_field
