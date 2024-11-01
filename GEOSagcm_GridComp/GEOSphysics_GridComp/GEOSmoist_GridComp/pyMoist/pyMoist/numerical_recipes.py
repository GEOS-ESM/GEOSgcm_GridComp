import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import exp, log, f64


@gtscript.function
def GammLn(xx: f64) -> f64:
    """
    See numerical recipes, w. press et al., 2nd edition.

    Compute the natural logarithm of the gamma function for a given value xx.

    Parameters:
    xx (Float in): Input value for which the natural logarithm of the
                gamma function is to be computed.

    Returns:
    Float: The natural logarithm of the gamma function value for the input xx.
    """
    stp = f64(2.5066282746310005)

    x = xx
    y = x
    tmp = x + f64(5.5)
    tmp = (x + f64(0.5)) * log(tmp) - tmp
    ser = f64(1.000000000190015)

    ser += f64(76.18009172947146) / (y + 1)
    ser += f64(-86.50532032941677) / (y + 1)
    ser += f64(24.01409824083091) / (y + 1)
    ser += f64(-1.231739572450155) / (y + 1)
    ser += f64(0.001208650973866179) / (y + 1)
    ser += f64(-0.000005395239384953) / (y + 1)

    return tmp + log(stp * ser / x)


@gtscript.function
def gser(a: f64, x: f64, gln: f64) -> f64:
    """
    See numerical recipes, w. press et al., 2nd edition.

    Compute the series representation of the incomplete gamma function.

    Parameters:
    a (Float in): Parameter a for the incomplete gamma function.
    x (Float in): Parameter x for the incomplete gamma function.
    gln (Float in): Natural logarithm of the gamma function.

    Returns:
    Float: The series representation of the incomplete gamma function.
    """
    eps = f64(3.0e-9)  # was eps=3.0d-07 in press et al.
    itmax = 10000  # was itmax=100   in press et al.
    gln = GammLn(a)
    gamser: f64 = f64(0.0)
    if x <= f64(0):
        # Fortran messages here x < 0 in gser
        # TODO: Allow print in GT4Py
        # if x < 0:
        # raise ValueError('aero_actv: subroutine gser: x < 0 in gser')
        gamser = f64(0.0)
    else:
        ap = a
        sum_ = f64(1.0) / a
        del_ = sum_
        n = 0
        while n < itmax:
            ap += f64(1.0)
            del_ *= x / ap
            sum_ += del_
            if abs(del_) < abs(sum_) * eps:
                gamser = sum_ * exp(-x + a * log(x) - gln)
                n = itmax
            n += 1
        gamser = sum_ * exp(-x + a * log(x) - gln)
    return gamser


@gtscript.function
def gcf_matrix(a: f64, x: f64, gln: f64) -> f64:
    """
    See numerical recipes, w. press et al., 2nd edition.

    Compute the continued fraction representation of the incomplete gamma function.

    Parameters:
    a (Float in): Parameter a for the incomplete gamma function.
    x (Float in): Parameter x for the incomplete gamma function.
    gln (Float in): Natural logarithm of the gamma function.

    Returns:
    Float: The continued fraction representation of the incomplete gamma function.
    """
    itmax = 10000
    eps: f64 = f64(3.0e-7)
    fpmin: f64 = f64(1.0e-30)
    gln = GammLn(a)
    b: f64 = x + f64(1.0) - a
    c: f64 = f64(1.0) / fpmin
    d: f64 = f64(1.0) / b
    h: f64 = d

    i = 1
    while i <= itmax:
        an = -i * (i - a)
        b += f64(2.0)
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = f64(1.0) / d
        del_ = d * c
        h *= del_
        if abs(del_ - f64(1.0)) < eps:
            i = itmax + 1
        i += 1
    return exp(-x + a * log(x) - gln) * h


@gtscript.function
def GammP(a: f64, x: f64) -> f64:
    """
    See numerical recipes, w. press et al., 2nd edition.

    Compute the incomplete gamma function for given values a and x.

    Parameters:
    a (Float in): Parameter a for the incomplete gamma function.
    x (Float in): Parameter x for the incomplete gamma function.

    Returns:
    Float: The incomplete gamma function value for the input parameters a and x.
    """
    # Fortran messages here potential bad arguments
    # TODO: Allow print in GT4Py
    # if (x < 0.0) or (a <= 0.0):
    #    raise ValueError("aero_actv: function gammp: bad arguments")
    gln = GammLn(a)
    if x < a + f64(1.0):
        gammp = gser(a, x, gln)
    else:
        gammp = f64(1.0) - gcf_matrix(a, x, gln)
    return gammp


@gtscript.function
def Erf(x: f64) -> f64:
    """
    See numerical recipes, w. press et al., 2nd edition.

    Compute the error function for a given value x.

    Parameters:
    x (Float in): Input value for which the error function is to be computed.

    Returns:
    Float: The error function value for the input x.
    """
    erf: f64 = f64(0.0)
    if x < f64(0.0e00):
        erf = f64(-1.0) * GammP(f64(0.5), x**2)
    else:
        erf = GammP(f64(0.5), x**2)
    return erf
