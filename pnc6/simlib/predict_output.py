import numpy as np
import matplotlib.pyplot as plt

import scipy as sp

from scipy.integrate import ode, odeint
from scipy.interpolate import interp1d

twopi = 2*np.pi


class Setup_pitch(object):

    def __init__(self, 
                tlist = None, 
                system_params = None,
                ):
        if system_params is None:
            from pnc6.simlib import system_params # ; reload(system_params)
            self.system_params = system_params.system_params
        else:
            self.system_params = system_params

        self.tlist = tlist
        self.tstep = self.tlist[1] - self.tlist[0]
        a0 = 1/np.sqrt(2)
        b0 = 0.0
        self.y0 = [a0, b0]

        self.call_solver()
        self.initialize_parameters()
        # print 'Specify drives for all values in tlist'


    def interpolate(self):
        drive_a = interp1d(self.tlist, self.drive_a)
        self.get_drive_a = lambda t: drive_a(t)
        drive_b = interp1d(self.tlist, self.drive_b)
        self.get_drive_b = lambda t: drive_b(t)


    def get_xi_a(self, t):
        na = self.system_params['na']#photon number calibration
        return np.sqrt(na)*self.get_drive_a(t)

    def get_xi_b(self, t):
        nb = self.system_params['nb']#photon number calibration
        return np.sqrt(nb)*self.get_drive_b(t)


    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_aa*np.abs(self.get_xi_a(t))**2 + chi_ab*np.abs(self.get_xi_b(t))**2
        delta_g = self.system_params['delta_g']
        return stark_shift - delta_g

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2 + chi_ab*np.abs(self.get_xi_a(t))**2
        delta_r = self.system_params['delta_r']
        return stark_shift - delta_r

    def get_g(self, t):
        kappa = self.system_params['kappa']
        chi_ab = self.system_params['chi_ab']
        # delta_g = self.system_params['delta_g']
        # delta_r = self.system_params['delta_r']
        return chi_ab*self.get_xi_a(t)*self.get_xi_b(t).conjugate()#*np.exp(-1j*(delta_g - delta_r)*t)

    def eom(self, t, y, params):
        #equations of motion for the system
        a, b = y
        kappa = params['kappa']
        
        adot = -1*self.get_g(t)*b - 1j*self.get_delta_a(t)*a
        bdot = self.get_g(t).conjugate()*a - 1j*self.get_delta_b(t)*b - kappa/2*b
        
        derivs = [adot, bdot]      # list of derivatives
        return derivs

    def get_eom(self):
        #returns a function (not a bound method) that can be passed to ode
        return lambda t, y, args: self.eom(t, y, args)

    def call_solver(self):
        self.solver = ode(self.get_eom()).set_integrator('zvode')
        self.solver.set_initial_value(self.y0, self.tlist[0])

    def initialize_parameters(self):
        self.solver.set_f_params(self.system_params)


    def iterate_solver(self):
        sol_list = []
        sol_list.append(self.y0) #start with initial value

        self.interpolate()

        for t in self.tlist[1:]:
            sol_list.append( self.solver.integrate(t)) #append solution at next timestep
        self.sol = np.array(sol_list)
        self.a = self.sol[:,0]
        # delta_r = self.system_params['delta_r']
        kappa = self.system_params['kappa']
        self.b_out = np.sqrt(kappa)*self.sol[:,1]#*np.exp(-1j*delta_r*self.tlist[0:])




class Compute_pitch(Setup_pitch):

    def __init__(self, 
                nout = 0.5, 
                **kwargs
                ):
        self.nout = nout
        super(Compute_pitch, self).__init__(**kwargs)
        self.drive_a = None

    
    def interpolate(self):
        drive_b = interp1d(self.tlist, self.drive_b)
        self.get_drive_b = lambda t: drive_b(t)
        b_out = interp1d(self.tlist, self.b_out)
        self.get_b_out = lambda t: b_out(t)
        if self.drive_a is not None:
            raise Exception('no drive_a')

    def get_xi_a(self, t):
        index = int(t/self.tstep)
        ret = self.sol_list[index]

    def objective(self, ders, t, xi_a):
        #return expression to be minimized
        I_prime, Q_prime = ders
        xi_a_prime = I_prime + 1j*Q_prime
        I, Q = xi_a.real, xi_a.imag
        g = self.get_g(xi_a, xi_b, params)
        g_star = g.conjugate()
        g_star_prime = self.get_g(xi_a_prime, xi_b, params).conjugate()
        xi_b = self.get_xi_b(t)

    
        b, b_prime, b_prime_prime = get_bs(t, params)
    
        delta_a = self.get_delta_a(xi_a, xi_b, params)
        delta_b = self.get_delta_b(xi_a, xi_b, params)
        delta_b_prime = get_delta_b_prime(xi_a, xi_a_prime, params)
    
        expr =  -g_star*b_prime_prime \
            + (g_star_prime - 1j*delta_a*g_star)*( b_prime + (kappa/2 + 1j*delta_b)*b ) \
            -g_star*g*g_star*b \
            -(kappa/2 + 1j*delta_b)*b_prime*g_star \
            -1j*delta_b_prime*g_star*b
            
        return [expr.real, expr.imag]


    def get_objective(self):
        #handler for objective function
        return lambda t, ders, args: self.objective(t, ders, args)


    def ders(self, t, y, args):
        #function that returns the derivative to ode solver
        last_xi = self.sol_list[-1]
    
        guess = (y - last_xi)/self.tstep
        xi_a = y[0] + 1j*y[1]
        #numerically solve for derivatives
        result = root(self.get_objective, guess.real, args=(t, xi_a, params))

        der = result['x']
        return der


    def get_ders(self):
        #handler for differential equation derivatives
        return lambda t, y, args: self.ders(t, y, args)


    def call_solver(self):
        self.solver = ode(self.get_ders()).set_integrator('zvode')
        self.solver.set_initial_value(self.y0, self.tlist[0])

    def iterate_solver(self):
        self.sol_list = [self.y0] #start with initial value

        for t in self.tlist[1:]:
            self.sol_list.append( self.solver.integrate(t)) #append solution at next timestep
        self.sol = np.array(self.sol_list)
        self.a = self.sol[:,0]
        # delta_r = self.system_params['delta_r']
        kappa = self.system_params['kappa']
        self.b_out = np.sqrt(kappa)*self.sol[:,1]#*np.exp(-1j*delta_r*self.tlist[0:])
