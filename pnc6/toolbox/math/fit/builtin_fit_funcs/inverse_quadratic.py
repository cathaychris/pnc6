import numpy as np

def func(xs, a=1, b=0, c=0, x0=0):
    return a/(b**2 + (xs-x0)**2) + c

def guess(xs, ys):
    p = np.polyfit(xs, 1./ys, 2)
    return dict(a=1./p[0], b=1./p[1], c=1./p[2], x0=(p[2]-p[1]**2/(4*p[0])))

