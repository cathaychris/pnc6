# import logging
# import time
import numpy as np

from scipy import polyval, fliplr, conjugate
from scipy.misc import factorial
from scipy.signal import convolve2d
import scipy.linalg as la

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

try:
    import qutip as qp
except ImportError:
    qp = None

from toolbox.plotting import mpltools


def Ws(xs_in, ys_in, zs_in, s_in, eta, normalization='int'):
    """
    Calculate the s-parametrized W-function according to (Leonhardt, 1997).
    Note that the factor of 1/2 in the gaussian is because we choose g=2
    in the the qutip wigner and qfunc functions (gives the vacuum we convolve
    with the correct sigma).

    Normalization is either
    - int : normalize integral to one
    - sum : normalize sum to one
    - any other value : returns just the convolution
    """
    s_out = (s_in - 1. + eta)/eta
    if s_out >= s_in:
        raise Exception("s_in must be larger than s_out")

    xx, yy = np.meshgrid(xs_in, ys_in)
    gaussian = np.exp(-(xx**2 + yy**2)/((s_in-s_out)/2))/(np.pi * (s_in-s_out))
    W = convolve2d(zs_in, gaussian, boundary='symm', mode='same')
    if normalization == 'sum':
        W /= W.sum()
    elif normalization == 'int':
        W /= np.trapz(np.trapz(W, x=xs_in), x=ys_in)
    return W


def qfunc(state, alpha):
    scalar = False
    if type(alpha) is not np.ndarray:
        alpha = np.array([alpha])
        scalar = True
    qmat = np.zeros(alpha.size)

    if qp.isket(state):
        qmat = _qfunc_pure(state, alpha)
    elif qp.isoper(state):
        d, v = la.eig(state.full())
        qmat = np.zeros(np.shape(alpha))
        for k in np.arange(0, len(d)):
            qmat1 = _qfunc_pure(v[:, k], alpha)
            qmat += (d[k] * qmat1).real

    if scalar:
        return float(qmat)
    return qmat


def _qfunc_pure(psi, alpha_mat):
    n = np.prod(psi.shape)
    if isinstance(psi, qp.qobj.Qobj):
        psi = psi.full().flatten()
    else:
        psi = psi.T

    qmat = abs(polyval(fliplr([psi / np.sqrt(factorial(np.arange(n)))])[0],
                       conjugate(alpha_mat))) ** 2

    return qmat.real * np.exp(-np.abs(alpha_mat) ** 2) / np.pi
