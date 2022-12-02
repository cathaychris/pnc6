# # -*- coding: utf-8 -*-
# """
# Created on Fri Jan 13 11:00:16 2017

# @author: Luke
# """
# import numpy as np
# twopi = 2*np.pi

# #units - GHz and ns
# MHz = 1e-3
# kHz = 1e-6


# kappa = 690*kHz*twopi
# chi_aa = -25.9*kHz*twopi
# chi_ab = -13.5*kHz*twopi
# chi_bb = -1.*kHz*twopi
# delta_g = -0*kHz*twopi
# delta_r = -304*kHz*twopi
# b_amp = 0.8#DAC units

# na = 8.56
# nb = 41.2


# # Bundle parameters for ODE solver
# params = {
#     'kappa': kappa,
#     'chi_aa': chi_aa,
#     'chi_ab': chi_ab,
#     'chi_bb': chi_bb,
#     'na': na,#photons per DAC unit
#     'nb': nb,#photons per DAC unit
#     'delta_g': delta_g,
#     'delta_r': delta_r,
#     'b_amp': b_amp
# }


import numpy as np
twopi = 2*np.pi

system_params = {}

system_params['kappa1'] = twopi*5.0
system_params['kappa2'] = twopi*5.0

system_params['chi_ab1'] = -twopi*0.04
system_params['chi_aa1'] = -twopi*0.02
system_params['chi_bb1'] = -twopi*0.02

system_params['chi_ab2'] = -twopi*0.04
system_params['chi_aa2'] = -twopi*0.02
system_params['chi_bb2'] = -twopi*0.02


system_params['xi_b1'] = np.sqrt(20)
system_params['xi_b2'] = np.sqrt(20)

#detuning between output modes: omega_b1 - omega_b2
system_params['delta_s'] = twopi*0.0 