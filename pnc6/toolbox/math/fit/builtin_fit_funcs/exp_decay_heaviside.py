import numpy as np
import common

def func(xs, x0=0, A=1, tau=1, ofs=0):
    def heaviside(x):
        return 0.5 * (np.sign(x) + 1)
    return heaviside(xs-x0) * (A * np.exp(-(xs-x0) / tau) + ofs) + heaviside(-xs+x0) * (A + ofs)

def guess(xs, ys):
    yofs = ys[-1]
    ys = ys - yofs
    return dict(
        x0=xs[np.argmax(np.array(ys)<0.9*(np.max(ys)-np.min(ys))+np.min(ys))],
        A=ys[0],
        tau=xs[common.find_index_of(ys, ys[0]/2)],
        ofs=yofs,
    )

TEST_RANGE = 0, 100
TEST_PARAMS = dict(A=-5, tau=10)
