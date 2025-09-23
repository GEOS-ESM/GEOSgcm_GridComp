"""Stencils and functions called by multiple pyMoist modules.
These functions manupilate data to formats more friendly to
gt4py and/or deal with gt4py I/O shortcomings."""

from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField, FloatFieldIJ


def get_last(in_field: FloatField, temporary_field: FloatFieldIJ, out_field: FloatField):
    """Workaround for getting last value, e.g. `Field[max(K)]`"""
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
    """Workaround for absolute index in K and relative in I/J,
    e.g. `Float[0,0,Absolute(k_index_desired)]`"""
    with computation(FORWARD), interval(...):
        if k_mask == k_index_desired:
            out_field = data_field
