import pytest

from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.warm_rain import get_v3_autoconversion_adjustments


@pytest.mark.parametrize(
    ("eis", "convection_fraction", "expected_cpaut_factor"),
    [
        (0.0, 0.0, 1.0),
        (5.0, 0.25, 0.875),
        (10.0, 1.0, 0.5),
        (20.0, 1.0, 0.5),
    ],
)
def test_v3_autoconversion_adjustments(eis, convection_fraction, expected_cpaut_factor):
    rthreshs = 10.0e-6
    rthreshu = 7.0e-6

    cpaut_factor, fac_rc = get_v3_autoconversion_adjustments(
        eis,
        convection_fraction,
        rthreshs,
        rthreshu,
    )

    expected_fac_rc = constants.RC * (
        rthreshu * convection_fraction + rthreshs * (1.0 - convection_fraction)
    ) ** 3
    assert cpaut_factor == pytest.approx(expected_cpaut_factor)
    assert fac_rc == pytest.approx(expected_fac_rc)
