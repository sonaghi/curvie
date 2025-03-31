from math import erf, sqrt, exp
from scipy.special import erfinv

SQRT = 1.4142135623730951  # \sqrt(2)
TRQS = 0.7071067811865475  # 1 / \sqrt(2)
RTTP = 0.3989422804014327  # 1 / \sqrt(2 * \pi)

def cdf(x):

    return 0.5 * (1 + erf(x * TRQS))

def quantile(x):

    return SQRT * erfinv(2 * x - 1)

def pdf(x):

    return RTTP * exp(-0.5 * x * x)

    
