from curvie.black.normal import cdf, pdf, quantile
from math import log, exp, sqrt, fabs

TOL64 = 2.2204460492503131E-16

def pair(r, v):

    if v > 0.:
        d = log(r) / v + 0.5 * v
        return (r * cdf(d) - cdf(d - v), r * pdf(d))
    elif r > 1:
        return (r - 1, 0.)
    elif r < 1:
        return (0., 0.)
    else:
        return (0., pdf(0.))

def base(r, v):

    if v > 0.:
        d = log(r) / v + 0.5 * v
        return r * cdf(d) - cdf(d - v)
    elif r > 1:
        return r - 1
    else:
        return 0.

def grad(r, v):

    if v > 0.:
        return cdf(log(r) / v + 0.5 * v)
    elif r > 1:
        return 1
    else:
        return 0.

def hess(r, v):

    if v > 0.:
        return pdf(log(r) / v + 0.5 * v) / (r * v)
    else:
        return 0.

def vega(r, v):

    if v > 0:
        d = log(r) / v + 0.5 * v
        return r * pdf(d)
    elif r == 1:
        return pdf(0.)
    else:
        return 0.

def cintro(f, k):

    return (f - k) if f > k else 0
    
def cprice(f, k, v):

    return base(f / k, v) * k

def cdelta(f, k, v):

    return grad(f / k)

def cgamma(f, k, v):

    return hess(f / k) / k

def cvegas(f, k, v):

    return vega(f / k, v) * k

def cdetox(f, d, v):

    return f * exp(-quantile(d) * v + 0.5 * v * v)

def ctodev(f, k, p):

    r = f / k
    q = p / k
    v = 2 * float(quantile((1 + q) / (1 + r)))    
    b, a = pair(r, v)
    while fabs(b - q) > 10 * (1 + q) * TOL64:
        v += (q - b) / a
        b, a = pair(r, v)
    return v

def pintro(f, k):

    return cintro(k, f)

def pprice(f, k, v):

    return base(k / f, v) * f

def pdelta(f, x, v):

    return cdelta(f, k, v) - 1

def pgamma(f, x, v):

    return cgamma(f, k, v)

def pvegas(f, x, v):

    return cvegas(f, k, v)

def pdetox(f, d, v):

    return cdetox(f, 1 - d, v)
