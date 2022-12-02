import numpy as np

def func(xs, a=1, b=0):
    return np.abs(a * (xs - b))

def guess(xs, ys):
    p = np.polyfit(xs, ys, 2)
    bval = xs[np.argmin(ys)]
    return dict(a=p[1], b=bval)

