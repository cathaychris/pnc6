import numpy as np
import common

def func(xs, A=1, tau=1, ofs=0, pseudo = 0.5):
    return A /( np.exp(xs / tau) + pseudo) + ofs

def guess(xs, ys):
    yofs = ys[-1]
    ys = ys - yofs
    pseudo = 0.5
    return dict(
        A=ys[0],
        tau=xs[common.find_index_of(ys, ys[0]/2)],
        ofs=yofs,
        pseudo=pseudo
    )

TEST_RANGE = 0, 100
TEST_PARAMS = dict(A=-5, tau=10)
