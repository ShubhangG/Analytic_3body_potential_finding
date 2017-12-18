import numpy as np
from scipy.optimize import curve_fit

def func(x, *p): return p[0] + p[1] * x

popt, pcov = curve_fit(func, np.arange(10), np.arange(10), p0=(0, 0))
print popt,pcov