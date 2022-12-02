import numpy as np
import matplotlib.pyplot as plt
import sine

def func(xs, A=1, f=0.05, dphi=np.pi/4, ofs=0, tau=0.5, tau2=2):
    return A * np.sin(2*np.pi*xs*f + dphi) * np.exp(-xs / tau) + b*xs + ofs

def guess(xs, ys):
    d = sine.guess(xs, ys)
    d['tau'] = np.average(xs)
    d['b'] = (np.average(ys[3*int(len(a)/4):])-np.average(ys[:int(len(a)/4)]))/(xs[-1]-xs[0])
    return d
