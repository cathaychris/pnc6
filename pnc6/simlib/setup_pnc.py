import numpy as np
from qutip import *



def rotating(t, det):
        return np.exp(-1j*t*det)

class Setup_pnc(object):
    #base class for generating hamiltonian coefficients for master equation
    '''
    TODO
     - implement and test corrections for detuning (and stark shifts)
    '''

    def __init__(self, 
                dim = 2,
                tlist = None, 
                det1 = 0.0,
                det2 = 0.0,
                amp_scale1 = None,
                amp_scale2 = None,
                waves = None,
                system_params = None,
                stark = True,
                chi = True,
                kerr = True
                ):

        self.dim_a = dim
        self.dim_b = 2
        self.tlist = tlist
        self.det1 = det1
        self.det2 = det2
        if system_params is None:
            from pnc6.simlib.system_params import system_params #; reload(system_params)
        self.system_params = system_params.system_params
        if stark:
            self.stark = 1.
        else:
            self.stark = 0

        if chi:
            self.chi = 1.
        else:
            self.chi = 0

        if kerr:
            self.kerr = 1.
        else:
            self.kerr = 0

        self.create_ops()
        # self.create_eops()
        
        
        self.amp_scale1 = amp_scale1
        self.amp_scale2 = amp_scale2
        self.waves = waves
        self.load_gs()


    def create_ops(self):

        a=destroy(self.dim_a)
        ad = a.dag()
        b=destroy(self.dim_b)
        bd = b.dag()
        Ia = qeye(self.dim_a)
        Ib = qeye(self.dim_b)

        self.a1 = tensor(
                    tensor(a, Ib),
                    tensor(Ia, Ib)
                    )
        self.b1 = tensor(
                    tensor(Ia, b),
                    tensor(Ia, Ib)
                    )
        self.a2 = tensor(
                    tensor(Ia, Ib),
                    tensor(a, Ib)
                    )
        self.b2 = tensor(
                    tensor(Ia, Ib),
                    tensor(Ia, b)
                    )

        #collape operator
        self.c1 = np.sqrt(self.system_params['kappa1'])*self.b1
        self.c2 = np.sqrt(self.system_params['kappa2'])*self.b2
 
    
        #hamiltonian terms
        #coupling
        self.Hc1 = 1j*np.sqrt(self.system_params['kappa1']*self.system_params['kappa2'])/2*(
                        tensor(
                            tensor(Ia, bd),
                            tensor(Ia, b)
                        )
                    )
        self.Hc2 = self.Hc1.dag()


        #beamsplitters - need also the hermitian conjugate of these
        self.H1a = -1j*tensor(
                tensor(a, bd),
                tensor(Ia, Ib)
            )
        self.H1b = self.H1a.dag()
        self.H2a = -1j*tensor(
                tensor(Ia, Ib),
                tensor(a, bd) 
            )
        self.H2b = self.H2a.dag()


        #physical detuning hamilotninans
        self.Hda1 = tensor(
                tensor(ad*a, Ib),
                tensor(Ia, Ib)
            )
        self.Hdb1 = tensor(
                tensor(Ia, bd*b),
                tensor(Ia, Ib)
            )
        self.Hda2 = tensor(
                tensor(Ia, Ib),
                tensor(ad*a, Ib)
            )
        self.Hdb2 = tensor(
                tensor(Ia, Ib),
                tensor(Ia, bd*b)
            )


        #expectation operators
        self.na1 = tensor(
                tensor(ad*a , Ib),
                tensor(Ia, Ib)
                )
        self.nb1 = tensor(
                tensor(Ia, bd*b ),
                tensor(Ia, Ib)
                )

        self.na2 = tensor(
                tensor(Ia, Ib),
                tensor(ad*a , Ib))
        self.nb2 = tensor(
                tensor(Ia, Ib),
                tensor(Ia, bd*b ))
        
    def load_gs(self):
        if self.waves is None:
            self.waves = np.load('gwaves.npz')
        self.g1w = self.waves['g1']
        self.g2w = self.waves['g2']
        if self.amp_scale1 is not None: self.g1w *= self.amp_scale1
        if self.amp_scale2 is not None: self.g2w *= self.amp_scale2
    
    def g1(self, t):
        return np.interp(t,self.tlist, self.g1w)*rotating(t, self.det1)
    def g2(self, t):
        return np.interp(t,self.tlist, self.g2w)*rotating(t, self.det2)
    
    
    def xi_a1(self, t):
        return self.g1(t)/(self.system_params['chi_ab1']*self.system_params['xi_b1'])
    def xi_a2(self, t):
        return self.g2(t)/(self.system_params['chi_ab2']*self.system_params['xi_b2'])



    def stark_a1(self, t):
        return 2*self.system_params['chi_aa1']*np.abs(self.xi_a1(t))**2 + self.system_params['chi_ab1']*np.abs(self.system_params['xi_b1'])**2
    def stark_a2(self, t):
        return 2*self.system_params['chi_aa2']*np.abs(self.xi_a2(t))**2 + self.system_params['chi_ab2']*np.abs(self.system_params['xi_b2'])**2
    def stark_b1(self, t):
        return self.system_params['chi_ab1']*np.abs(self.xi_a1(t))**2 + 2*self.system_params['chi_ab1']*np.abs(self.system_params['xi_b1'])**2
    def stark_b2(self, t):
        return self.system_params['chi_ab2']*np.abs(self.xi_a2(t))**2 + 2*self.system_params['chi_ab2']*np.abs(self.system_params['xi_b2'])**2



    def collapse_ops(self, t, args):
        return [self.c1*rotating(t, self.system_params['delta_s']) + self.c2*rotating(t, 0)]
    def get_collapse_ops(self):
        return lambda t, args: self.collapse_ops(t, args)

    def hamiltonian(self, t, args):
        return self.stark*(self.stark_a1(t)*self.a1.dag()*self.a1 + self.stark_b1(t)*self.b1.dag()*self.b1 + 
                    self.stark_a2(t)*self.a2.dag()*self.a2 + self.stark_b2(t)*self.b2.dag()*self.b2 ) + \
                self.chi*( self.system_params['chi_ab1']*self.a1.dag()*self.a1*self.b1.dag()*self.b1  + 
                    self.system_params['chi_ab2']*self.a2.dag()*self.a2*self.b2.dag()*self.b2 ) + \
                self.kerr*( self.system_params['chi_aa1']*self.a1.dag()**2*self.a1**2 + self.system_params['chi_bb1']*self.b1.dag()**2*self.b1**2 + 
                    self.system_params['chi_aa2']*self.a2.dag()**2*self.a2**2 + self.system_params['chi_bb2']*self.b2.dag()**2*self.b2**2 ) + \
                self.g1(t)*self.H1a + np.conj(self.g1(t))*self.H1b + self.g2(t)*self.H2a + np.conj(self.g2(t))*self.H2b  + \
                self.Hc1*rotating(-t, self.system_params['delta_s']) + self.Hc2*rotating(t, self.system_params['delta_s'])

    def get_hamiltonian(self):
        return lambda t, args: self.hamiltonian(t, args)


    def get_L(self):
        return lambda t, args: liouvillian(self.hamiltonian(t, args), self.collapse_ops(t, args)).data



class Setup_self_res_catch(Setup_pnc):
    def __init__(self, detuning, **kwargs):
        self.detuning_phys1 = 0
        self.detuning_phys2 = detuning
        self.det1 = 0
        self.det2 = 0
        super(Setup_self_res_catch, self).__init__(**kwargs)

        
class Setup_pitch_detuned_catcher(Setup_pnc):
    def __init__(self, detuning, **kwargs):
        self.detuning_phys1 = 0
        self.detuning_phys2 = detuning
        self.det1 = 0
        self.det2 = -detuning
        super(Setup_pitch_detuned_catcher, self).__init__(**kwargs)
    
class Setup_pitch_half_detuned(Setup_pnc):
    def __init__(self, detuning, **kwargs):
        self.detuning_phys1 = -detuning/2.0
        self.detuning_phys2 = detuning/2.0
        self.det1 = detuning/2
        self.det2 = -detuning/2
        super(Setup_pitch_half_detuned, self).__init__(**kwargs)


class Setup_pitch_detuned_pitcher(Setup_pnc):
    def __init__(self, detuning, **kwargs):
        self.detuning_phys1 = -detuning
        self.detuning_phys2 = 0
        self.det1 = detuning
        self.det2 = 0
        super(Setup_pitch_detuned_pitcher, self).__init__(**kwargs)
    
    
class Setup_catch_detuning_sweeper(Setup_pnc):
    def __init__(self, detuning_phys, detuning, **kwargs):
        self.detuning_phys1 = 0
        self.detuning_phys2 = detuning_phys
        self.det1 = 0
        self.det2 = detuning
        super(Setup_catch_detuning_sweeper, self).__init__(**kwargs)
    
    
class Setup_catch_stark_shifts(Setup_pnc):
    def __init__(self, detuning_phys, detuning, **kwargs):
        self.detuning_phys1 = 0
        self.detuning_phys2 = detuning_phys
        self.det1 = 0
        self.det2 = detuning
        super(Setup_catch_detuning_sweeper, self).__init__(**kwargs)


