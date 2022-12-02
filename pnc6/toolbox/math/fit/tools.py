import numpy as np
import lmfit

import fitter

DEBUG = False


def lmfit_params_str(params, fmt=':1.2E', stderr=True, delimiter='\n'):
    txt = ''
    for ip, p in enumerate(params):
        txt += ("{} = {" + fmt + "}").format(p, params[p].value)
        if stderr:
            txt += (" +/- {" + fmt + "}").format(params[p].stderr)

        if ip < len(params)-1:
            txt += delimiter

    return txt


def multi_fit(xs, ys, func_name, **kw):
    print_report = kw.pop('print_report', False)
    mask = kw.pop('mask', None)

    ret = {}
    ret['fit_yvals'] = np.zeros(ys.shape) * np.nan
    ret['fit_success'] = np.zeros(ys[...,0].shape, dtype=bool)

    it = np.nditer(ys[...,0], flags=['multi_index'])
    while not it.finished:
        idxs = it.multi_index
        y = ys[idxs]

        # check if we're supposed to fit this trace
        if mask is not None and not mask[idxs]:
            it.iternext()
            continue

        # simple check if data is actually available
        if np.isnan(y[0]):
            it.iternext()
            continue

        if xs.shape == ys.shape:
            x = xs[idxs]
        else:
            x = xs

        fit = fitter.Fitter(func_name)

        try:
            res = fit.perform_lmfit(x, y, print_report=print_report, **kw)
        except:
            if DEBUG:
                raise
            it.iternext()
            continue
        ret['fit_yvals'][idxs] = fit.eval_func()

        for pn in res.params:
            if pn not in ret:
                ret[pn] = np.zeros(ys[...,0].shape)
                ret[pn+"_stderr"] = np.zeros(ys[...,0].shape)
            ret[pn][idxs] = res.params[pn].value
            ret[pn+"_stderr"][idxs] = res.params[pn].stderr

        ret['fit_success'][idxs] = True
        it.iternext()

    return ret