# -*- coding: utf-8 -*-
"""
Created on Mon Oct 06 14:15:10 2014

@author: rsl
"""

import numpy as np
import matplotlib.pyplot as plt
# import cavity_analysis
from lmfit import minimize, Parameters
from numpy import random
from scipy import optimize
import lmfit
import copy

def quick_lor_fit(freq, db_mag, show=False, bkg=False, grads_are_data=False):
    db2v = from_dB
    fit_v_lorentzian(freq, db_mag-np.min(db_mag), show=True)
    # fit_v_lorentzian(freq, db2v(db_mag)-np.median(db2v(db_mag)), show=True)
    # return quick_hanger_fit(freq, db_mag, show=show, bkg=bkg, grads_are_fitdata=True)

def quick_hanger_fit(freq, db_mag, show=False, bkg=False, grads_are_data=False, grads_are_fitdata=False, extra_info='', returnparams=False):
    x = freq
    y0 = db_mag
    if bkg:
        slope = np.linspace(y0[0],y0[-1],np.size(x))
        y0 = np.subtract(y0,slope)
    params, fitdata = lmfit_asymmetric_db_hanger(x, y0)
    f0 = params['f0'].value
    qi = params['qi'].value
    qcr = params['qcr'].value
    qci = params['qci'].value
    bw = f0/qi
    fit_label = 'center = %0.5f GHz\n' % (f0/1e9)
    fit_label += 'Qi = %0.3e\n' % (qi)
    fit_label += 'Qc = %0.3e (%0.3e + j%0.3e)\n' % (qcr,qcr,qci)
    fit_label += 'BW = %0.2f kHz' % (bw/1e3)
    fit_label += extra_info
    if not grads_are_fitdata:
        if grads_are_data:
            graddata = y0
        else:
            graddata = fitdata
        grad = np.gradient(np.array([x, graddata])) # derivative of fit
        grads = grad[1][1]
        grad = np.gradient(np.array([x, grads]))
        grads = np.abs(grad[1][1]) # second derivative
    else:
        grads = fitdata
    if show:
        from matplotlib import gridspec
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        fig = plt.figure()
        fig.add_subplot(gs[0])
        fig.add_subplot(gs[1])
        fig.axes[0].plot(x, y0, marker='s', ms=3, label='')
        fig.axes[0].plot(x, fitdata, ls='-', ms=3, label=fit_label)
        fig.axes[1].plot(x, grads, ls='-', ms=3, label='')
        fig.axes[0].legend(loc='best')
        fig.axes[0].set_xlabel('Frequency (Hz)')
        fig.axes[0].set_ylabel('Magnitude (dB)')
    if returnparams:
        return params, bw, grads
    else:
        return f0, qi, bw, grads # all in Hz
    
if __name__ == '__main__':
    # an example:
    testfname = r"Z:\_Data\Si_uMachine\201506_InStriplines8-14\InStrips12-14\In13_LT_resonanceOnly-30.00dB_22.41mK.dat"
    x, y0, y1 = np.loadtxt(testfname, unpack=True, delimiter='\t')
    fitout = quick_hanger_fit(x, y0, show=True)
    
def dumb_fit_v_lorentzian(freqs,data,show=False):
    params, fitdata = dumb_fit_v_lorentzian(freqs,np.abs(data)-np.min(np.abs(data)),show=show)
    f0 = params['f0'].value
    bw = params['fwhm'].value
    # qcr = params['qcr'].value
    # qci = params['qci'].value
    return params, bw # all in Hz


def to_dB(x):
    return 20*np.log10(x)


def from_dB(x):
    return 10**(x/20.)


def fit_hfss_transmission(filename, lossy=False, show=False):
    with open(filename) as f:
        headers = f.readline()[1:-1].split('","')
        data = np.transpose([map(float, line.split(',')) for line in f.readlines()])
        freqs = data[0]*1e9
        res = {}
        for name, trace in zip(headers[1:], data[1:]):
            cutoff = trace.max() / 1e2
            region = np.where(trace > cutoff)[0]
            imin, imax = region[0], region[-1]
            freq_clip = freqs[imin:imax]
            trace_clip = trace[imin:imax]
            if 'mag' in name:
                params = fit_v_lorentzian(freq_clip, trace_clip, show=show)
            else:
                params = fit_db_lorentzian(freq_clip, trace_clip, show=show)
            res[name] = {n: params[n].value for n in ("qi", "qc", "f0")}
            res[name]['q'] = 1/(1/res[name]['qc'] + 1/res[name]['qi'])
        return res


def open_text_data(filename, delim=None):
    with open(filename) as f:
        try:
            return np.transpose([map(float, line.split(delim)) for line in f.readlines()])
        except ValueError:
            f.seek(0)
            headers = f.readline()
            print headers
            return np.transpose([map(float, line.split(delim)) for line in f.readlines()])


def send_to_igor(arr, wavenames, overwrite=True):
    arr = np.array(arr)
    if len(arr.shape) == 1:
        arr = np.array([arr])
        wavenames = [wavenames]
    if len(wavenames) != arr.shape[1]:
        arr = arr.transpose()
    if len(wavenames) != arr.shape[1]:
        raise ValueError('wavenames not of length of the data array')
    import tempfile

    f = tempfile.NamedTemporaryFile(delete=False)
    fn = f.name
    fn = fn.replace('\\', ':').replace('::', ':')
    np.savetxt(f, arr, delimiter=',', header=','.join(wavenames))
    f.close()
    args = "/J/D/W/K=0/A"
    if overwrite:
        args += '/O'
    cmd = 'LoadWave%s "%s"'%(args, fn)
    igor_do(cmd)


def plot_in_igor(x, y, xname, yname, overwrite=False):
    xname = igorify_name(xname)
    yname = igorify_name(yname)
    send_to_igor([x, y], [xname, yname], overwrite=overwrite)
    igor_do(
        'Display %s vs %s'%(yname, xname),
        'Label left "%s"'%yname,
        'Label bottom "%s"'%xname
    )


def igor_do(*cmds):
    for cmd in cmds:
        print cmd
        import win32com.client

        igor = win32com.client.Dispatch("IgorPro.Application")
        igor.Execute(cmd)


def igorify_name(n):
    if len(n) > 32:
        print n, 'is too long shortening to 32 chars'
        n = n[:31]
    for c in set(n):
        if not c.isalnum():
            n = n.replace(c, '_')
    while n.endswith('_'):
        n = n[:-1]
    print n
    return n


def plot_complex_data(freqs, mags, phases):
    s21 = mags*np.exp(1j*np.pi*phases/180.)
    s21inv = 1/s21
    plt.figure(2)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.plot(np.real(s21inv), np.imag(s21inv))
    plt.axes().set_aspect('equal')
    plt.axes().grid(True)

def do_hangerfit(p0, model, x, y=None):
    if y is None:
        return p0
    
    fit = lmfit.minimize(model, p0, args=[x,y])
    print lmfit.fit_report(fit)
    fit_yvals = model(fit.params, x)
    return fit.params, fit_yvals

def do_lmfit(xpts, data_pts, fit_func, params, show=False):
    param_names = params.keys()
    initial_params = copy.deepcopy(params)
    fit = lmfit.minimize(fit_func, params, args=[xpts, data_pts])
    print lmfit.fit_report(fit)
 
    def apply_params(params, xs):
        values = {s: params[s].value for s in param_names}
        return fit_func(xs, **values)
    
    fitted_data = fit_func(fit.params, xpts)

    if show:
        for n in param_names:
            print n, 'initial: ', initial_params[n].value, 'fitted: ', params[n].value
        print data_pts
        plt.plot(xpts, data_pts, label='data')
        plt.plot(xpts, apply_params(initial_params, xpts), label='initial')
        plt.plot(xpts, apply_params(params, xpts), label='fitted')

        plt.legend()
        if show is True:
            plt.show()
        else:
            plt.title(show)
            plt.savefig("fit_plot_%s.png"%show)
            plt.clf()

    return fit.params, fitted_data

def do_fit(xpts, data_pts, fit_func, params, show=False, force=False):
    show=True
    """
    params is a dictionary of dictionaries containing optionally value, min, and max fields
    if show is enabled, matplotlib is used to compare the data, initial guess, and fitted result
    """
    param_names = params.keys()
    is_complex = data_pts.imag.any()
    if is_complex:
        print 'Complex Data Found'

    def apply_params(params, xs):
        values = {s: params[s].value for s in param_names}
        return fit_func(xs, **values)

    if is_complex:
        residual = lambda params, xs, ys: np.abs(ys - apply_params(params, xs))**2
    else:
        residual = lambda params, xs, ys: ys - apply_params(params, xs)

    print residual
    print params

    if isinstance(params, Parameters):
        nlm_params = params
    else:
        nlm_params = Parameters()
        for name, kwargs in params.items():
            nlm_params.add(name, **kwargs)
    initial_params = copy.deepcopy(nlm_params)

    if not force and (apply_params(initial_params, xpts) == np.nan).any():
        raise ValueError("Function produced NaNs on initial params")

    # min_result = minimize(residual, nlm_params, args=(xpts, data_pts))
    print xpts.shape
    print data_pts.shape
    from scipy import optimize
    min_result = minimize(residual, nlm_params, args=(xpts, data_pts))
    print min_result.lmdif_message
    print min_result.success, 'success'
    print min_result.message
    print min_result.nfev, 'evaluations'
    fitted_data = apply_params(nlm_params, xpts)

    if show:
        for n in param_names:
            print n, 'initial: ', initial_params[n].value, 'fitted: ', nlm_params[n].value
        if is_complex:
            plt.figure(0)
            plt.plot(xpts, np.abs(data_pts), label='data')
            plt.plot(xpts, np.abs(apply_params(initial_params, xpts)), label='initial')
            plt.plot(xpts, np.abs(apply_params(nlm_params, xpts)), label='fitted')
            plt.legend()
            plt.figure(1)
            plt.plot(xpts, np.angle(data_pts), label='data')
            plt.plot(xpts, np.angle(apply_params(initial_params, xpts)), label='initial')
            plt.plot(xpts, np.angle(apply_params(nlm_params, xpts)), label='fitted')
            plt.legend()
        else:
            print data_pts
            plt.plot(xpts, data_pts, label='data')
            plt.plot(xpts, apply_params(initial_params, xpts), label='initial')
            plt.plot(xpts, apply_params(nlm_params, xpts), label='fitted')

        plt.legend()
        if show is True:
            plt.show()
        else:
            plt.title(show)
            plt.savefig("fit_plot_%s.png"%show)
            plt.clf()

    return nlm_params, fitted_data


def param(value, min, max):
    'A parameter which can be passed as a keyword argument to do_fit'
    return {'value': value, 'min': min, 'max': max}

def complex_v_hanger(f, f0, qi, qc, scale):
    q = 1/(1/qi + 1/qc)
    x = (f - f0)/f0
    return scale*(1 - (q/qc) / (1 + 2j*q*x))

def v_hanger(f, f0, qi, qc, scale):
    return np.abs(complex_v_hanger(f, f0, qi, qc, scale))

def model_asymmetric_complex_v_hanger(p, f, y=None):
    f0 = p['f0'].value
    qi = p['qi'].value
    qcr = p['qcr'].value
    qci = p['qci'].value
    scale = p['scale'].value
    phase=0
    qc = qcr + 1j*qci
    x = (f - f0)/f0
    model = np.exp(1j*phase)*scale*qc*(2*qi*x - 1j)/(2*qi*qc*x - 1j*(qi + qc))
#    return model
    if y is None:
        return model
    else:
        return model - y

def asymmetric_complex_v_hanger(f, f0, qi, qcr, qci, scale, phase=0):
    'See Kurtis: Asymmetric Hanger Equations'
    qc = qcr + 1j*qci
    x = (f - f0)/f0
    return np.exp(1j*phase)*scale*qc*(2*qi*x - 1j)/(2*qi*qc*x - 1j*(qi + qc))

def asymmetric_v_hanger(f, f0, qi, qcr, qci, scale):
    return np.abs(asymmetric_complex_v_hanger(f, f0, qi, qcr, qci, scale))

def model_asymmetric_v_hanger(p, f, y=None):
    return np.abs(model_asymmetric_complex_v_hanger(p, f, y=y))

def angle_asymmetric_v_hanger(f, f0, qi, qcr, qci, scale):
    return np.angle(asymmetric_complex_v_hanger(f, f0, qi, qcr, qci, scale))

def asymmetric_db_hanger(f, f0, qi, qcr, qci, offset):
    return to_dB(asymmetric_v_hanger(f, f0, qi, qcr, qci, from_dB(offset)))

def fit_v_hanger(f, s21, show=False):
    params = asymmetric_hanger_guess(f, s21)
    params['qc'] = params['qcr']
    params.pop('qcr')
    params.pop('qci')
    print params
    return do_fit(f, s21, v_hanger, params, show=show)

def fit_shifted_circle(x, y, x_m=0, y_m=0):
    'From: http://wiki.scipy.org/Cookbook/Least_Squares_Circle'

    def calc_R(xc, yc):
        """ calculate the distance of each data points from the center (xc, yc) """
        return np.sqrt((x - xc)**2 + (y - yc)**2)

    def f_2b(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    def Df_2b(c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc = c
        df2b_dc = np.empty((len(c), x.size))

        Ri = calc_R(xc, yc)
        df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
        df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
        df2b_dc = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

        return df2b_dc

    center_estimate = x_m, y_m
    center_2b, ier = optimize.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)

    xc_2b, yc_2b = center_2b
    Ri_2b = calc_R(*center_2b)
    R_2b = Ri_2b.mean()
    #residu_2b    = sum((Ri_2b - R_2b)**2)
    return xc_2b, yc_2b, R_2b

def asymmetric_hanger_guess(fpts, vpts):
    scale_guess = (vpts[0] + vpts[-1]) / 2.
    vpts /= scale_guess
    i0 = np.argmin(vpts)
    f0_guess = fpts[i0]
    delta_i_guess = min(i0, len(fpts) - i0)
    delta_f_guess = delta_i_guess*(fpts[1] - fpts[0])/10
    q_total_guess = f0_guess/delta_f_guess
    qi_guess = q_total_guess/vpts[i0]
    qc_guess = 1/(1/q_total_guess - 1/qi_guess)
    vpts *= scale_guess

    assert all(q > 0 for q in (q_total_guess, qi_guess, qc_guess))
    assert not (asymmetric_complex_v_hanger(fpts, f0_guess, qi_guess, qc_guess, 0, scale_guess) == np.nan).any()

    return {
        'f0': param(f0_guess, fpts[0], fpts[-1]),
        'qi': param(qi_guess, qi_guess/10., qi_guess*10),
        'qcr': param(qc_guess, qc_guess/10., qc_guess*10),
        'qci': param(0, -qc_guess, qc_guess),
        'scale': param(scale_guess, scale_guess/2., 2*scale_guess)
    }
    
def lmfit_asymmetric_hanger_guess(fpts, vpts):
    scale_guess = (vpts[0] + vpts[-1]) / 2.
    vpts /= scale_guess
    i0 = np.argmin(vpts)
    f0_guess = fpts[i0]
    delta_i_guess = min(i0, len(fpts) - i0)
    delta_f_guess = delta_i_guess*(fpts[1] - fpts[0])/10
    q_total_guess = f0_guess/delta_f_guess
    qi_guess = q_total_guess/vpts[i0]
    qc_guess = 1/(1/q_total_guess - 1/qi_guess)
    vpts *= scale_guess

    assert all(q > 0 for q in (q_total_guess, qi_guess, qc_guess))
    assert not (asymmetric_complex_v_hanger(fpts, f0_guess, qi_guess, qc_guess, 0, scale_guess) == np.nan).any()

    p0 = lmfit.Parameters()
    p0.add('f0', f0_guess, vary=True, min=fpts[0], max=fpts[-1])
    p0.add('qi', qi_guess, vary=True, min=qi_guess/50., max=qi_guess*50.)
    p0.add('qcr', qc_guess, vary=True, min=qc_guess/10., max=qc_guess*10.)
    p0.add('qci', qc_guess/10., vary=True, min=-qc_guess/2., max=qc_guess*2.)
    p0.add('scale', scale_guess, vary=True, min=scale_guess/4., max=4.*scale_guess)

    return p0
    
def dumb_angle_guess(fpts, vpts):
    scale_guess = (vpts[0] + vpts[-1]) / 2.
    # vpts /= scale_guess
    i0 = np.argmax(vpts)
    f0_guess = fpts[i0]
    delta_i_guess = min(i0, len(fpts) - i0)
    delta_f_guess = delta_i_guess*(fpts[1] - fpts[0])/10
    q_total_guess = f0_guess/delta_f_guess
    qi_guess = q_total_guess/2
    qc_guess = q_total_guess/2
    # vpts *= scale_guess

    assert all(q > 0 for q in (q_total_guess, qi_guess, qc_guess))
    assert not (asymmetric_complex_v_hanger(fpts, f0_guess, qi_guess, qc_guess, 0, scale_guess) == np.nan).any()

    return {
        'f0': param(f0_guess, fpts[0], fpts[-1]),
        'qi': param(qi_guess, qi_guess/10., qi_guess*10),
        'qcr': param(qc_guess, qc_guess/10., qc_guess*10),
        'qci': param(0, -qc_guess, qc_guess),
        'scale': param(scale_guess, scale_guess/2., 2*scale_guess)
    }

def lmfit_asymmetric_v_hanger(fpts, vpts, show=False):
    fpts, vpts = np.array(fpts), np.array(vpts)
    params = lmfit_asymmetric_hanger_guess(fpts, vpts)
    return do_hangerfit(params, model_asymmetric_v_hanger, fpts, y=vpts)
    # return do_lmfit(fpts, vpts, model_asymmetric_v_hanger, params, show=show)

def fit_asymmetric_v_hanger(fpts, vpts, show=False):
    fpts, vpts = np.array(fpts), np.array(vpts)
    params = asymmetric_hanger_guess(fpts, vpts)
    return do_fit(fpts, vpts, asymmetric_v_hanger, params, show=show)

def fit_angle_asymmetric_v_hanger(fpts, vpts, show=False):
    fpts, vpts = np.array(fpts), np.array(vpts)
    params = dumb_angle_guess(fpts, vpts)
    return do_fit(fpts, vpts, angle_asymmetric_v_hanger, params, show=show)

def lmfit_asymmetric_db_hanger(fpts, db_pts, show=False):
    params, fitted_data = lmfit_asymmetric_v_hanger(fpts, from_dB(db_pts), show)
    return params, to_dB(fitted_data)

def fit_asymmetric_db_hanger(fpts, db_pts, show=False):
    params, fitted_data = fit_asymmetric_v_hanger(fpts, from_dB(db_pts), show)
    return params, to_dB(fitted_data)

def fit_angle_asymmetric_db_hanger(fpts, deg_pts, show=False):
    # params, fitted_data = fit_angle_asymmetric_v_hanger(fpts, from_dB(db_pts), show)
    params, fitted_data = fit_angle_asymmetric_v_hanger(fpts, deg_pts, show)
    # return params, to_dB(fitted_data)
    return params, fitted_data

def complex_v_lorentzian(f, f0, qi, qc):
    return 1/(1 + (qc/qi) - 2j*qc*(f - f0)/f0)

def dumb_v_lorentzian(f, f0, fwhm, a):
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    return a*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(f-f0)**2/(2*sigma**2))

def v_lorentzian(f, f0, qi, qc):
    return np.abs(complex_v_lorentzian(f, f0, qi, qc))


def db_lorentzian(f, f0, qi, qc):
    return to_dB(v_lorentzian(f, f0, qi, qc))


def lorentzian_guess(freqs, s21):
    i0 = np.argmax(s21)
    f0_guess = freqs[i0]
    fspan = freqs[-1] - freqs[0]
    q_guess = 10*f0_guess/fspan
    qc_guess = q_guess/s21[i0]
    #assert qc_guess > q_guess
    qi_guess = 1/(1/q_guess - 1/qc_guess)
    #check_relerror(s21[i0], power_lorentzian(f0_guess, f0_guess, qi_guess, qc_guess), .05)

    return {
        'f0': param(f0_guess, freqs[0], freqs[-1]),
        'qc': param(qc_guess, qc_guess/100, qc_guess*100),
        'qi': param(qi_guess, qi_guess/100, qi_guess*100),
    }

def dumb_lorentzian_guess(freqs, v):
    i0 = np.argmax(v)
    f0_guess = freqs[i0]
    fspan = freqs[-1] - freqs[0]
    fwhm_guess = fspan/10

    return {
        'f0': param(f0_guess, freqs[0], freqs[-1]),
        'fwhm': param(fwhm_guess, fwhm_guess/10, fwhm_guess*100),
        'a': param(fwhm_guess, fwhm_guess/1000, fwhm_guess*1000),
    }

def fit_v_lorentzian(freqs, s21, show=False):
    params = lorentzian_guess(freqs, s21)
    return do_fit(freqs, s21, v_lorentzian, params, show=show)

def fit_db_lorentzian(freqs, s21, show=False):
    params, fitted_data = fit_v_lorentzian(freqs, from_dB(s21), show)
    return params, to_dB(fitted_data)
    
def dumb_fit_v_lorentzian(freqs, y, show=False):
    params = dumb_lorentzian_guess(freqs, y)
    return do_fit(freqs, y, dumb_v_lorentzian, params, show=show)

def check_relerror(x0, x1, tolerance):
    error = (x0 - x1)/x0
    assert (x0 - x1)/x0 < tolerance, "error %s between %s,%s greater than %s"%(error, x0, x1, tolerance)


def test_fit_asymm_db_hanger():
    fpts = np.linspace(9.1e9, 9.101e9, 500)
    f0 = 9.1005e9
    qi = 4e5
    qcr = 1e6
    qci = 5e5
    offset = -3
    scale = from_dB(offset)
    fake = asymmetric_db_hanger(fpts, f0, qi, qcr, qci, offset) + random.normal(0, .05, len(fpts))
    result = fit_asymmetric_db_hanger(fpts, fake, show=True)
    check_relerror(f0, result['f0'].value, .08)
    check_relerror(qi, result['qi'].value, .08)
    check_relerror(qcr, result['qcr'].value, .08)
    check_relerror(qci, result['qci'].value, .08)
    check_relerror(scale, result['scale'].value, .08)

def test_fit_lorentzian():
    fpts = np.linspace(9.1e9, 9.101e9, 500)
    f0 = 9.1005e9
    qi = 4e5
    qc = 1e6

    fake = db_lorentzian(fpts, f0, qi, qc) + random.normal(0, .05, len(fpts))
    result = fit_db_lorentzian(fpts, fake, show=True)
    check_relerror(f0, result['f0'].value, .08)
    check_relerror(qi, result['qi'].value, .08)
    check_relerror(qc, result['qc'].value, .08)

def photon_number(power_dbm, q_int, t0, f0):
    hbar = 1.0546e-34
    s21 = np.sqrt(t0)
    omega = 2*np.pi*f0
    power_watts = from_dB(power_dbm)*1e-3
    return power_watts*q_int*s21*(1 - s21)/(hbar*omega**2)