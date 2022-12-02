import numpy as np
import common

def func(xs, ofs=0.0, A=1.0, x0=0.0, delta=2.0, n=1.0, Tc=5000.0):
    # assuming an on resonant ramsey expt.
    # assuming no cavity decay

    loss = np.exp(-xs/Tc)
    exp_term = np.exp(n*(1-loss**2-0.5*loss*np.cos(2*np.pi*delta*(xs-x0))))
    exp_term = np.exp(-n*(1-np.cos(2*np.pi*delta*(xs-x0))))
    return A * exp_term + ofs

def guess(xs, ys):
    return dict(
        delta = np.abs(4.5 / (xs[-1] - xs[0])),
        x0 = -2./np.abs(4.5 / (xs[-1] - xs[0])),
        n = 2.0,
        A = (max(ys) - min(ys)),
        Tc = 10000.0,
#        ofs = -(max(ys) - min(ys))/2.0,
        ofs = 0.0,
    )