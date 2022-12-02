import numpy as np

def func(xs, **kw):
    chi = kw['chi']
    chi_prime = kw.get('chi_prime', 0)
    x0 = kw.get('x0', 0)
    sigma = kw.get('sigma', 1)
    ofs = kw.get('ofs', 0)
    As = {}
    for k in kw:
        if k[0] == 'A':
            As[int(k[1:])] = kw[k]

    def gaussian(xs, x0, sigma):
        return np.exp(-(xs-x0)**2/(2*sigma**2))

    model = ofs
    for iA in As:
        model += As[iA] * gaussian(xs, x0 + iA*chi + iA*(iA-1.)*chi_prime/2., sigma)
    return model
