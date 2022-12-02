import numpy as np
import common

def func(xs, ofs=0, x0=0, sigma=2, A=-1, k = 5, sigma_k = 1, A_k = -1):
    '''
    Sum of two lorentzian dips separated by k/2
    '''
    def lorentzian(x, cntr, width):
        return width**2/((x-cntr)**2 + width*2)

    return ofs + A*lorentzian(xs, x0, sigma) + A_k*lorentzian(xs, x0+k/2, sigma_k)

def guess(xs, ys):
    yofs = max(ys)
    ys = ys - yofs
    A = min(ys)
    minidx = np.argmin(ys)

    sigma = (xs[-1] - xs[1])/50.

    k = -sigma*5
    sigma_k = sigma/2
    x0=xs[minidx]-k

    return dict(
        ofs=yofs,
        x0=x0,
        k=k,
        sigma=sigma,
        sigma_k = sigma_k,
        A = A,
        A_k = A
    )



TEST_RANGE = 0, 10
TEST_PARAMS = dict(ofs=0, x0=0, sigma=2, A=-1, k = 5, sigma_k = 1, A_k = -1)
