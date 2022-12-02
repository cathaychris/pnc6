import numpy as np
import matplotlib.pyplot as plt

from pnc6.toolbox import h5tools

from pnc6.simlib.setup_helpers import *

import scipy as sp
import time
import copy

from scipy.integrate import ode, odeint, quad
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import root
import lmfit
import warnings



twopi = 2*np.pi

TIME_FMT = '%H:%M:%S'

class Setup_pitch(object):

    def __init__(self,
                tlist = None,
                system_params = None
                ):
        if system_params is None:
            from pnc6.simlib.system_params import system_params #; reload(system_params)
            self.system_params = system_params.system_params
        else:
            self.system_params = system_params

        self.tlist = tlist
        self.tstep = self.tlist[1] - self.tlist[0]

        if 'fit_g_params' in self.system_params:
            # print("Found fit params for g, no need to fit again.")
            # p = self.fit_g_params = self.system_params['fit_g_params']
            self.fit_g_params = self.system_params['fit_g_params']
            self.get_fit_g = lambda p, x, y: model_2dg(p, x, y)
        else:
            self.fit_g_info = self.system_params.get('fit_g_info', None)
            #fit_g_info takes the form (fn, date, time) from calibration experiment sweeper
            self.fit_g()
        self.call_solver()

    def fit_g(self):
        kappa = self.system_params['kappa']
        if self.fit_g_info is not None:
            fn, date, time = self.fit_g_info
            self.get_fit_g, self.fit_g_params = fit_g_from_file(fn, 'g_calibration', date, time, 'stamps', 'roamps', 'keffs', kappa)

    def get_xi_a(self, t):
        n = self.system_params['na']#photon number calibration
        xi = np.sqrt(n)*self.get_drive_a(t)
        return xi

    def get_xi_b(self, t):
        n = self.system_params['nb']#photon number calibration
        xi = np.sqrt(n)*self.get_drive_b(t)
        return xi

    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_aa*np.abs(self.get_xi_a(t))**2 + chi_ab*np.abs(self.get_xi_b(t))**2
        delta_a = self.system_params['delta_a']
        return stark_shift - delta_a

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2 + chi_ab*np.abs(self.get_xi_a(t))**2
        delta_b = self.system_params['delta_b']
        return stark_shift - delta_b

    def get_g(self, t):
        chi_ab = self.system_params['chi_ab']
        if hasattr(self, 'get_fit_g'):
            a = self.get_drive_a(t)
            b = self.get_drive_b(t)
            g = -self.get_fit_g(self.system_params['fit_g_params'], a, b)
            # amp_a = np.abs(a)
            # amp_b = np.abs(b)
            # phase_a = np.angle(a)
            # phase_b = np.angle(b)

            # mag_g = self.get_fit_g(amp_a, amp_b)
            # g = -np.exp(1j*(phase_a - phase_b))*mag_g
        else:
            g = chi_ab*self.get_xi_a(t)*self.get_xi_b(t).conjugate()
        return g



    def eom(self, t, y):
        pass

    def get_eom(self):
        #returns a function (not a bound method) that can be passed to ode
        return lambda t, y: self.eom(t, y)

    def call_solver(self):
        self.solver = ode(self.get_eom()).set_integrator('zvode')

    def initialize_parameters(self):
        self.solver.set_initial_value(self.sol_list[0], self.tlist[0])
        # self.solver.set_f_params(self.system_params)

    def iterate_solver(self):
        self.interpolate()
        self.get_initial_values()


        self.initialize_parameters()


        for t in self.tlist[1:]:
            self.sol_list.append(None)
            self.sol_list[-1] = self.solver.integrate(t)  #append solution at next timestep
        self.sol = np.array(self.sol_list)

        self.process()

    def process(self):
        pass






class Predict_pitch(Setup_pitch):
    def __init__(self, **kwargs):
        super(Predict_pitch, self).__init__(**kwargs)

    def get_initial_values(self):
        a0 = 1.
        b0 = 0.0
        y0 = [a0, b0]
        self.a0 = a0
        self.sol_list = [y0]

    def interpolate(self):
        pad = 200
        drive_a = interp_pad(self.tlist, self.drive_a, pad)
        self.get_drive_a = lambda t: drive_a(t)
        drive_b = interp_pad(self.tlist, self.drive_b, pad)
        self.get_drive_b = lambda t: drive_b(t)


    def eom(self, t, y):
        #equations of motion for the system
        a, b = y
        kappa = self.system_params['kappa']

        adot = -1*self.get_g(t)*b - 1j*self.get_delta_a(t)*a
        bdot = self.get_g(t).conjugate()*a - 1j*self.get_delta_b(t)*b - kappa/2*b

        derivs = [adot, bdot]      # list of derivatives
        return derivs


    def process(self):
        #give the solution sensible names
        self.a = self.sol[:,0]
        kappa = self.system_params['kappa']
        self.b_out = np.sqrt(kappa)*self.sol[:,1]



class Predict_catch(Setup_pitch):
    def __init__(self, **kwargs):
        super(Predict_catch, self).__init__(**kwargs)

    def get_initial_values(self):
        a0 = 0.0
        b0 = 0.0
        y0 = [a0, b0]
        self.sol_list = [y0]

    def interpolate(self):
        pad = 200
        drive_a = interp_pad(self.tlist, self.drive_a, pad)
        self.get_drive_a = lambda t: drive_a(t)
        drive_b = interp_pad(self.tlist, self.drive_b, pad)
        self.get_drive_b = lambda t: drive_b(t)
        b_in = interp_pad(self.tlist, self.b_in, pad)
        self.get_b_in = lambda t: b_in(t)


    def eom(self, t, y):
        #equations of motion for the system
        a, b = y
        kappa = self.system_params['kappa']

        adot = -1*self.get_g(t)*b - 1j*self.get_delta_a(t)*a
        bdot = self.get_g(t).conjugate()*a - 1j*self.get_delta_b(t)*b - kappa/2*b - np.sqrt(kappa)*self.get_b_in(t)

        derivs = [adot, bdot]      # list of derivatives
        return derivs


    def process(self):
        #give the solution sensible names
        self.a = self.sol[:,0]
        kappa = self.system_params['kappa']
        self.b = self.sol[:,1]
        self.b_out = np.sqrt(kappa)*self.b + self.b_in



class Compute_pitch(Setup_pitch):
    """
    needs to be fed an array of data points b_out_shape
    which is the (arbitrarily scaled) output field
    This class correctly scales it to be normalized to n_out photons (n_out <1)
    This means that the total energy out is n_out*n_total (for |0> + |1>, n_total is 0.5)

    Also fed an array of data points drive_b (usually constant)
    """
    def __init__(self,
                nout = 0.5,
                simulate=False,
                **kwargs
                ):
        self.nout = nout
        super(Compute_pitch, self).__init__(**kwargs)
        self.drive_a = None#np.zeros_like(self.tlist, dtype=np.complex128)

        self.simulate = simulate
        self.sim_type = Predict_pitch
        self.sim_kwargs = {'drive_a' : self.drive_a}


    def interpolate(self):
        pad = 200
        drive_b = interp_pad(self.tlist, self.drive_b, pad)
        self.get_drive_b = lambda t: drive_b(t)

        b_out = interp_pad(self.tlist, self.b_out_shape, pad)
        self.get_b_out_shape = lambda t: b_out(t)
        self.normalize_b_out()

        if self.drive_a is not None:#np.any(self.drive_a)
            drive_a = interp_pad(self.tlist, self.drive_a, pad)
            self.get_drive_a = lambda t: drive_a(t)



    def get_b_out_squared(self):
        return lambda t: np.abs(self.get_b_out_shape(t))**2

    def normalize_b_out(self):
        # self.get_b_out = normalize_b(self.get_b_out_shape(self.tlist), self.tlist, self.nout)
        scale = np.trapz(np.abs(self.get_b_out_shape(self.tlist))**2, self.tlist)
        self.get_b_out = lambda t: np.sqrt(self.nout)*self.get_b_out_shape(t)/np.sqrt(scale) #relative to a(0)


    def b_p(self, t):
        kappa = self.system_params['kappa']
        return self.get_b_out(t)/np.sqrt(kappa)

    def b_p_prime(self, t):
        tl = self.tlist[0]
        tr = self.tlist[-1]
        return f_der(self.b_p, t, self.tstep, tl, tr)

    def b_p_prime_prime(self, t):
        tl = self.tlist[0]
        tr = self.tlist[-1]
        return f_der(self.b_p_prime, t, self.tstep, tl, tr)

    def get_bs(self, t):
        return self.b_p(t), self.b_p_prime(t), self.b_p_prime_prime(t)

    def get_drive_a(self, t):
        # if this method is called, it has not get been overwritten, which means we are still integrating
        # the current value of xi_a lives in self.sol_list[-1]
        na = self.system_params['na']
        index = -1 # int(t/self.tstep)
        xi_v = self.sol_list[index]
        drive = (xi_v[0] + 1j*xi_v[1])/na**0.5
        return drive

    def get_g_prime(self, xi_a, xi_b, xi_a_prime):
        chi_ab = self.system_params['chi_ab']
        na = self.system_params['na']
        nb = self.system_params['nb']
        if hasattr(self, 'get_fit_g'):
            drive_a = xi_a/na**0.5
            drive_b = xi_b/nb**0.5
            drive_a_prime = xi_a_prime/na**0.5
            # amp_a = abs(xi_a_prime)/na**.5
            # amp_b = abs(xi_b)/nb**.5
            # phase_a = np.angle(xi_a_prime)
            # phase_b = np.angle(xi_b)
            g_prime = -model_g_prime(self.fit_g_params, drive_a, drive_b, drive_a_prime)
        else:
            g_prime = chi_ab*xi_a_prime*xi_b.conjugate()
        return g_prime


    # def get_xi_a(self, t):
    #     na = self.system_params['na']
    #     if self.drive_a is None:
    #         #integrating. The current value of xi_a lives in self.sol_list[-1]
    #         index = -1#int(t/self.tstep)
    #         xi_v = self.sol_list[index]
    #         xi = xi_v[0] + 1j*xi_v[1]
    #     else:
    #         xi = super(Compute_pitch, self).get_xi_a(t)
    #     return xi

    def get_delta_b_prime(self, t, xi_a_prime):
        chi_ab = self.system_params['chi_ab']
        xi_a = self.get_xi_a(t)
        I, Q = xi_a.real, xi_a.imag
        I_prime, Q_prime = xi_a_prime.real, xi_a_prime.imag
        stark_shift_derivative = 2*chi_ab*( I*I_prime + Q*Q_prime )
        return stark_shift_derivative

    def objective(self, ders, args):
        t, xi_a_v = args
        # xi_a = xi_a_v[0] + 1j*xi_a_v[0]
        # return expression to be minimized
        kappa = self.system_params['kappa']
        chi_ab = self.system_params['chi_ab']
        I_prime, Q_prime = ders
        xi_a_prime = I_prime + 1j*Q_prime
        xi_a = self.get_xi_a(t)
        # I, Q = xi_a.real, xi_a.imag
        xi_b = self.get_xi_b(t)
        g = self.get_g(t)#chi_ab*xi_a*xi_b.conjugate()#self.get_g(t)
        g_star = g.conjugate()
        # g_prime = chi_ab*xi_a_prime*xi_b.conjugate()
        g_prime = self.get_g_prime(xi_a, xi_b, xi_a_prime)
        g_star_prime = g_prime.conjugate()


        b, b_prime, b_prime_prime = self.get_bs(t)

        delta_a = self.get_delta_a(t)
        delta_b = self.get_delta_b(t)
        delta_b_prime = self.get_delta_b_prime(t, xi_a_prime)

        expr =  -g_star*b_prime_prime \
            + (g_star_prime - 1j*delta_a*g_star)*( b_prime + (kappa/2 + 1j*delta_b)*b ) \
            -g_star*g*g_star*b \
            -(kappa/2 + 1j*delta_b)*b_prime*g_star \
            -1j*delta_b_prime*g_star*b

        return [expr.real, expr.imag]


    def get_objective(self):
        # handler for objective function
        return lambda ders, args: self.objective(ders, args)


    def eom(self, t, y):
        # function that returns the derivative to ode solver
        last_xi = self.sol_list[-2]

        guess = (y - last_xi)/self.tstep
        xi_a = y[0] + 1j*y[1]
        self.sol_list[-1] = y
        # numerically solve for derivatives
        result = root(self.get_objective(), guess.real, args=[t, y])

        der = result['x']
        return der

    def get_initial_values(self):
        guess = [0.001,0.001]
        result = root(self.get_ini_fun(), guess)
        xi_a0 = result['x']
        self.sol_list = [xi_a0]


    def ini_fun(self, xi_a):
        kappa = self.system_params['kappa']
        # chi_ab = self.system_params['chi_ab']
        # xi_b0 = self.get_xi_b(0)
        xi_a_z = xi_a[0] + 1j*xi_a[1]
        self.sol_list = [xi_a]
        d_b0 = self.get_delta_b(0)
        b0, b_prime0, _ = self.get_bs(0)
        g0_star = self.get_g(0).conjugate() # chi_ab*xi_a_z.conjugate()*xi_b0
        expr = -g0_star + (kappa/2 + 1j*d_b0)*b0 + b_prime0 # takes a(0) = 1
        return [expr.real, expr.imag]

    def get_ini_fun(self):
        return lambda xi_a: self.ini_fun(xi_a)

    def process(self):
        na = self.system_params['na']
        self.drive_a = (self.sol[:,0] + 1j*self.sol[:,1])/np.sqrt(na)
        self.interpolate()
        if self.simulate:
            t0 = time.time()
            print('Running simulation at {}'.format(time.strftime(TIME_FMT, time.localtime(t0))))
            self.sim_kwargs['drive_a'] = self.drive_a
            sim_runner  = PncRunner(
                            self.sim_type,
                            tlist=self.tlist,
                            system_params=self.system_params,
                            **self.sim_kwargs)
            sim_runner.run()
            self.sim = sim_runner.setup


class Compute_catch(Compute_pitch):
    """
    needs to be fed an array of data points b_out_shape
    which is the (arbitrarily scaled) output field
    This class correctly scales it to be normalized to n_out photons (n_out <1)
    This means that the total energy out is n_out*n_total (for |0> + |1>, n_total is 0.5)

    Also fed an array of data points drive_b (usually constant)
    """
    def __init__(self,
                **kwargs
                ):
        nout = kwargs.pop('nout', 0.99)
        super(Compute_catch, self).__init__(nout = nout, **kwargs)
        self.sim_type = Predict_catch

    def objective(self, ders, args):
        t, xi_a_v = args
        # xi_a = xi_a_v[0] + 1j*xi_a_v[0]
        #return expression to be minimized
        kappa = -1*self.system_params['kappa']
        chi_ab = self.system_params['chi_ab']
        I_prime, Q_prime = ders
        xi_a_prime = I_prime + 1j*Q_prime
        xi_a = self.get_xi_a(t)
        # I, Q = xi_a.real, xi_a.imag
        xi_b = self.get_xi_b(t)
        g = self.get_g(t)#chi_ab*xi_a*xi_b.conjugate()#self.get_g(t)
        g_star = g.conjugate()
        # g_prime = chi_ab*xi_a_prime*xi_b.conjugate()
        g_prime = self.get_g_prime(xi_a, xi_b, xi_a_prime)
        g_star_prime = g_prime.conjugate()


        b, b_prime, b_prime_prime = self.get_bs(t)

        delta_a = self.get_delta_a(t)
        delta_b = self.get_delta_b(t)
        delta_b_prime = self.get_delta_b_prime(t, xi_a_prime)

        # this seems equivalent to Eq. S10
        expr =  -g_star*b_prime_prime \
            + (g_star_prime - 1j*delta_a*g_star)*( b_prime + (kappa/2 + 1j*delta_b)*b ) \
            -g_star*g*g_star*b \
            -(kappa/2 + 1j*delta_b)*b_prime*g_star \
            -1j*delta_b_prime*g_star*b

        return [expr.real, expr.imag]

    def eom(self, t, y):
        #function that returns the derivative to ode solver
        last_xi = self.sol_list[-2]

        guess = (y - last_xi)/self.tstep
        xi_a = y[0] + 1j*y[1]
        self.sol_list[-1] = y
        #numerically solve for derivatives
        result = root(self.get_objective(), guess.real, args=[t, y])

        der = result['x']
        return -1*der

    def ini_fun(self, xi_a):
        kappa = -1*self.system_params['kappa']
        # chi_ab = self.system_params['chi_ab']
        # xi_b0 = self.get_xi_b(0)
        T = self.tlist[-1]
        xi_a_z = xi_a[0] + 1j*xi_a[1]
        self.sol_list = [xi_a]
        d_b0 = self.get_delta_b(T) #the argument here doesn't actually change naything, just gets the current value
        b0, b_prime0, _ = self.get_bs(0)
        g0_star = self.get_g(T).conjugate()
        expr = -g0_star + (kappa/2 + 1j*d_b0)*b0 + b_prime0
        return [expr.real, expr.imag]


    def get_bs(self, t):
        T = self.tlist[-1]
        return super(Compute_catch, self).get_bs(T-t)

    def process(self):
        self.sol = np.flipud(self.sol)
        if self.simulate:
            self.sim_kwargs['b_in'] = self.b_out_shape#normalize_b(self.b_out_shape, self.tlist, self.nout)
        super(Compute_catch, self).process()


class Predict_pitch_simple(Predict_pitch):
    """
    No dynamic stark shifts
    """
    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift =  chi_ab*np.abs(self.get_xi_b(t))**2
        delta_a = self.system_params['delta_a']
        return stark_shift - delta_a

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2
        delta_b = self.system_params['delta_b']
        return stark_shift - delta_b



class Predict_catch_simple(Predict_catch):
    """
    No  dynamic stark shifts
    """
    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift =  chi_ab*np.abs(self.get_xi_b(t))**2
        delta_a = self.system_params['delta_a']
        return stark_shift - delta_a

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2
        delta_b = self.system_params['delta_b']
        return stark_shift - delta_b
    # def get_delta_a(self, t):
    #     delta_a = self.system_params['delta_a']
    #     return -delta_a

    # def get_delta_b(self, t):
    #     delta_b = self.system_params['delta_b']
    #     return -delta_b



class Compute_pitch_simple(Compute_pitch):
    """
    No stark shifts
    """

    def __init__(self, **kwargs):
        super(Compute_pitch_simple, self).__init__(**kwargs)
        self.sim_type = Predict_pitch_simple

    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift =  chi_ab*np.abs(self.get_xi_b(t))**2
        delta_a = self.system_params['delta_a']
        return stark_shift - delta_a

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2
        delta_b = self.system_params['delta_b']
        return stark_shift - delta_b

    def get_delta_b_prime(self, t, xi_a_prime):
        return 0.*t



class Compute_catch_simple(Compute_catch):
    """
    No stark shifts
    """

    def __init__(self, **kwargs):
        super(Compute_catch_simple, self).__init__(**kwargs)
        self.sim_type = Predict_catch_simple


    def get_delta_a(self, t):
        chi_aa = self.system_params['chi_aa']
        chi_ab = self.system_params['chi_ab']
        stark_shift =  chi_ab*np.abs(self.get_xi_b(t))**2
        delta_a = self.system_params['delta_a']
        return stark_shift - delta_a

    def get_delta_b(self, t):
        chi_bb = self.system_params['chi_bb']
        chi_ab = self.system_params['chi_ab']
        stark_shift = 2*chi_bb*np.abs(self.get_xi_b(t))**2
        delta_b = self.system_params['delta_b']
        return stark_shift - delta_b

    def get_delta_b_prime(self, t, xi_a_prime):
        return 0.*t


class PncRunner(object):
    '''
        setup_type is any of the extant children of the Setup_pitch class

        params is a dict of the form
        {
        param_name (string) : param__val,
        ...
        }

        '''
    def __init__(self, setup_type, initialize=True, **kwargs):
        self.setup_type = setup_type
        if initialize:
            self.initialize_setup(**kwargs)


    def initialize_setup(self, silent=True, **kwargs):
        #special attributes
        drive_a = kwargs.pop('drive_a', None)
        b_out_shape = kwargs.pop('b_out_shape', None)
        b_in = kwargs.pop('b_in', None)

        if not silent:
            t0 = time.time()
            print( 'Initializing setup at {}'.format(time.strftime(TIME_FMT, time.localtime(t0))) )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            setup = self.setup_type(
                                **kwargs
                                )

        if not silent:
            t1 = time.time()
            dt = t1-t0
            print( 'Initialization finished at {} in {} seconds'.format( time.strftime(TIME_FMT, time.localtime(t1)), int(dt) ) )

        if drive_a is not None:
            setup.drive_a = drive_a
        if b_out_shape is not None:
            setup.b_out_shape = b_out_shape
        if b_in is not None:
            setup.b_in = b_in

        setup.drive_b = setup.system_params['b_amp']*np.ones_like(setup.tlist)

        self.setup = setup
        self.silent = silent


    def run(self):
        if not self.silent:
            t0 = time.time()
            print( 'Starting solver at {}'.format(time.strftime(TIME_FMT, time.localtime(t0))) )

        self.setup.iterate_solver()

        if not self.silent:
            t1 = time.time()
            dt = t1-t0
            print( 'Solver finished at {} in {} seconds'.format( time.strftime(TIME_FMT, time.localtime(t1)), int(dt) ) )




class PncParamSweeper(object):


    def __init__(self, setup_type, params, silent_sweep=True, **kwargs):
        '''
        setup_type is any of the extant children of the Setup_pitch class

        params is a dict of the form
        {
        param_name (string) : param_range (list or one-dim numpy array),
        ...
        }
        If param_range is a single number or list of length 1, it will not be swept over

        '''
        self.setup_type = setup_type
        fit_g_info = params.pop('fit_g_info', None)
        self.params = params
        self.param_names = params.keys() #so that the parameters are always pulled in the same (arbitrary) order
        # self.b_out = b_out
        self.run_kws = kwargs
        self.run_kws['system_params'] = {'fit_g_info' : fit_g_info}
        self.silent_sweep = silent_sweep
        # self.run_kws['fit_g_info'] = fit_g_info



    def initialize_sweep(self):
        self.runner = PncRunner(self.setup_type, initialize=False)

        self.swept_params = {}
        for k,v in self.params.iteritems():
            v = np.array(v)
            self.swept_params[k] = v

        swp_vals = []
        for i_param, param_name in enumerate(self.param_names):
            # if param_name != 'fit_g_info':
            swp_vals.append(self.swept_params[param_name])
        grids = np.meshgrid(*swp_vals,indexing='ij')
        #make a list of 1-D arrays that contain the values for each parameter
        self.grids = [grid.reshape(-1) for grid in grids]

        self.n_steps = self.grids[0].size

        # self.swept_params = []
        # self.static_params = []
        # for k,v in self.params.iteritems():
        #     v = np.array(v)
        #     self.params[k] = v
        #     if len(v) > 1:
        #         self.swept_params.append(k)
        #     else:
        #         self.static_params.append

        # self.current_idxs = { param_name : 0  for param_name in self.swept_params }
        # self.final_idxs = { param_name : len(param_range)-1  for param_name, param_range in self.swept_params.iteritems() }




    # def params_from_idxs(self, idxs=None):
    #     if idxs = None:
    #         idxs = self.current_idxs

    #     params = {}
    #     for param_name, param_vals in self.params.iteritems():
    #         params[param] = param_vals[idxs[param_name]]

    #     return params
    def params_from_idx(self, idx=None):
        if idx is None:
            idx = self.current_step_idx

        _params = copy.deepcopy(self.params)
        for i_param, param_name in enumerate(self.param_names):
            # if param_name != 'fit_g_info':
            current_val = self.grids[i_param][idx]
            _params[param_name] = current_val

        return _params


    def initialize_step(self):
        current_params = self.params_from_idx()
        self.run_kws['system_params'].update(current_params)
        self.runner.initialize_setup(**self.run_kws)


    def do_step(self):
        t0 = time.time()
        if not self.silent_sweep:
            print( 'Running step {} of {}'.format(self.current_step_idx, self.n_steps) )
        self.runner.run()
        #save result somewhere
        res = self.runner.setup
        self.results.append(res)


    def run_next_step(self):
        self.initialize_step()
        self.do_step()


    def run_sweep(self):
        self.current_step_idx = 0
        self.results = []

        self.initialize_sweep()

        while self.current_step_idx < self.n_steps:
            self.run_next_step()
            self.current_step_idx += 1
        # while self.current_idxs != self.final_idxs:
        #     self.run_next()
