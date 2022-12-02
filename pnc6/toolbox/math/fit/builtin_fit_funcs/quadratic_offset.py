import numpy as np

def func(xs, a=1, b=1, c=0):
    return a * (xs-b)**2 + c

def guess(xs, ys):
    p = np.polyfit(xs, ys, 2)
    return dict(a=p[0], b=1.0, c=p[2]+1.0*p[1])

