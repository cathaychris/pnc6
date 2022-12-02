import time
import numpy as np

import scipy as sp
from scipy.integrate import ode, odeint, quad
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import root
import lmfit
from pnc6.toolbox import h5tools


PAD = 10


def take_f_der(f, t, dt, tl, tr, *args):
    if t - dt/2.0 < tl:
        der = (f(t + dt, *args) - f(t, *args))/dt
    elif t + dt/2.0 > tr:
        der = (f(t, *args) - f(t-dt, *args))/dt
    else:
        der = (f(t + dt/2, *args) - f(t - dt/2, *args))/dt
    return der

def f_der(f, t, *args):
    #return the derivative of a function f at time t
    #will break if f is defined on too narrow a domain
    #best if the function f is analytic rather than interpolated
    if isinstance(t, np.ndarray):
        der = np.array([])
        for tval in t:
            der = np.append(der, take_f_der(f, tval, *args))
    else:
        der = take_f_der(f, t, *args)
    return der

def interp_pad(ts, data, pad = PAD):
    #ode solver tries to reach too far for interpolatd functions
    #pulses will be zero outside shaped region anyway

    tpad_l = ts[0:pad] - ts[pad]
    tpad_r = ts[-pad:] + ts[pad]
    dpad_l = np.ones(pad)*data[0]
    dpad_r = np.ones(pad)*data[-1]
    ts_pad = np.concatenate((tpad_l, ts, tpad_r))
    data_pad = np.concatenate((dpad_l,data, dpad_r))

    interped = interp1d(ts_pad, data_pad)
    return interped

def normalize_b(b_out, ts, n_out):
    '''
    correctly scale a waveform
    '''
    scale = np.trapz(abs(b_out)**2, ts)/n_out
    b_out_scaled = b_out/scale**.5
    return b_out_scaled


def model_2dg(p,x,y,z=None):
    ofs = p['ofs'].value
    a11 = p['a11'].value
    a13 = p['a13'].value
    a31 = p['a31'].value
    # a22 = p['a22'].value

    model = ofs \
            + a11*(x)*(y) \
            + a13*(x)*(y*abs(y)**2) \
            + a31*(x*abs(x)**2)*(y) #\
            # + a22*(x**2)*(y**2)

    if z is None:
        return model
    return model - z

def model_g_prime(p, x, y, x_prime):
    ofs = p['ofs'].value
    a11 = p['a11'].value
    a13 = p['a13'].value
    a31 = p['a31'].value
    # a22 = p['a22'].value

    model = ofs \
            + a11*(x_prime)*(y) \
            + a13*(x_prime)*(y*abs(y)**2) \
            + 2*a31*(x_prime)*(abs(x)**2)*(y) + a31*(x_prime.conjugate())*(x**2)*(y) #\#\
            # + 2*a22*(x_prime)*(y**2)

    return model


def fit_2dg(x,y,z=None, vary_a13=False, vary_a31=False,):
    p0 = lmfit.Parameters()
    p0.add('ofs',   0,      vary=False)
    p0.add('a11',   0,    vary=True)
    p0.add('a13',   0,     vary=vary_a13)#True
    p0.add('a31',   0,     vary=vary_a31)
    # p0.add('a22',   0,      vary=False)

    if z is None:
        return p0

    fit = lmfit.minimize(model_2dg, p0, args=(x,y,z), nan_policy='omit')
    return fit


def fit_g(x, y, z, kappa, verbose=False):
    """
    args:
        - x : amp1 (DAC units)
        - y : amp2 (DAC units)
        - z : kappa effective (1/ns)
              should only have positive values and NAN
        - kappa : output mode decay rate (1/ns)
    returns:
        fit_func, lmfitparams
    """
    g = (z*kappa/4)**.5

    xx, yy = np.meshgrid(x, y, indexing='ij')
    fitres = fit_2dg(xx, yy, g)
    if verbose:
        print(lmfit.fit_report(fitres))
    ff = lambda x, y: model_2dg(fitres.params, x, y)
    return ff, fitres.params


def fit_g_from_file(file_name, msmt_name, date, time, x_name, y_name, z_name, kappa):
    data = h5tools.get_content_from_file(file_name)

    found_date = False
    for name, d in data.items():
        if msmt_name in name and date in name:
            cal_grp = d
            found_date = True
    if not found_date:
        raise KeyError('Data for that datestamp not found in calibration data file')

    found_time = False
    for name, d in cal_grp.items():
        if time == name:
            cal_exp = d
            found_time = True
    if not found_time:
        raise KeyError('Data for that timestamp not found in calibration data file')

    x = cal_exp[x_name]['value']
    y = cal_exp[y_name]['value']
    z = cal_exp[z_name]['value']/1e6 # k_eff in twopi*GHz
    z = np.maximum(z,0*z)

    return fit_g(x, y, z, kappa)

