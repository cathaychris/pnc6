"""
SSB tune up for both FPGA channels USING the SignalHound.
ATM, assumes that we're using the fpga_lib Calibratables.

Requires the following lab software:
* fpga_lib
* labpython (toolbox)

Todo/Comments:
* so far i've assumed only one AWG at most -- might need some adaption
* not checking if other generators are off (not sure that's a problem)
* so far we assume that we have at most 2 channels per FPGA card, and that for
  the second channel we use '2' in the sequences, whereas the instrument props
  for that channel have idx '1'.
"""

import time
import os
import numpy as np
from matplotlib import pyplot as plt

from fpga_lib.dsl import *
from fpga_lib.parameters import Calibratable, FloatParameter, IntParameter
from fpga_lib.experiments import FPGAExperiment

from toolbox import h5tools
from toolbox.math import minimizer
from toolbox.plotting import mpltools

from instrument_plugins import SA124B

#class SignalHound(object):
#    '''
#        Avoids using instrumentserver for the Signal Hound, since that can be troublesome
#    '''
#    def __init__(self):
#        self.open_sa()
#        self.first_sweep = True
#        self.f0 = 5e9
#        self.error_list = []
#        self.vbw, self.rbw = 250e3, 250e3
#        
#    def open_sa(self):
#        assert SA124B.saGetAPIVersion() == '3.0.11' # .15 doesn't play so nice
#        s_nums = SA124B.saGetSerialNumberList()
#        print s_nums
#        s_num = s_nums[0]
#        assert s_num != 0
#        err, did = SA124B.saOpenDeviceBySerialNumber(s_num)
#        print err
#        print "Device ID: ", did
#        assert did != 0
#        SA124B.saSetTimebase(did)
#        SA124B.saConfigAcquisition(did)
#        self.sa = did
#        
#    def close_sa(self):
#        SA124B.saCloseDevice(self.sa)
#
#    def _initiate(self):
#        SA124B.saInitiate(self.sa, SA124B.SA_SWEEPING)
#
#    def do_set_vbw(self, vbw):
#        err_code = SA124B.saConfigSweepCoupling(self.sa, self.rbw, vbw)
#        self._initiate()
#        if err_code == SA124B.NO_ERROR:
#            self.vbw = vbw
#        self.error_list.append(err_code)
#        return vbw
#    
#    def do_set_rbw(self, rbw):
#        err_code = SA124B.saConfigSweepCoupling(self.sa, rbw, self.vbw)
#        self._initiate()
#        if err_code == SA124B.NO_ERROR:
#            self.rbw = rbw
#        self.error_list.append(err_code)
#        return rbw
#        
#    def do_set_span(self, span):
#        center_frequency = self.do_get_center_frequency()
#        err_code = SA124B.saConfigCenterSpan(self.sa, center_frequency, span)
#        self._initiate()
#        self.error_list.append(err_code)
#        return span
#        
#    def do_get_center_frequency(self):
#        err_code, sweep_length, start_freq, bin_size = SA124B.saQuerySweepInfo(self.sa)
#        self.error_list.append(err_code)
#        return start_freq + round(sweep_length/2)*bin_size
#        
#    def do_get_span(self):
#        err_code, sweep_length, start_freq, bin_size = SA124B.saQuerySweepInfo(self.sa)
#        self.error_list.append(err_code)     
#        return sweep_length*bin_size
#        
#    def do_set_center_frequency(self, freq):
#        span = self.do_get_span()
#        err_code = SA124B.saConfigCenterSpan(self.sa, freq, span)
#        self._initiate()
#        self.error_list.append(err_code)
#        return freq
#        
#    def sweep(self, get_xs=True, sweeplen=None):
#        if get_xs:
#            err_code, xs, ys = SA124B.saGetSweep(self.sa, get_xs=True)
#            self.error_list.append(err_code)
#            return xs, ys
#        else:
#            assert not (sweeplen is None)
#            err_code, ys = SA124B.saGetSweep(self.sa, get_xs=False, sweeplen=sweeplen)
#            self.error_list.append(err_code)
#            return ys   
#    
#    def get_power(self):
#        z = self.first_sweep
#        val = SA124B.saReadSingleFreq(self.sa, self.f0, set_zeroif_params=z, verify_freq=z)
#        self.first_sweep = False
#        return val

#sh = SignalHound()
# sh = instruments['SignalHound']

class AWGSSBSettings(Calibratable):
    dphi = FloatParameter(0)
    Iamp = FloatParameter(1)
    Qamp = FloatParameter(1)
    Iof = FloatParameter(0)
    Qof = FloatParameter(0)
    if_period = IntParameter(-20)

    def __init__(self, name):
        super(AWGSSBSettings, self).__init__(name)

# class AWGSSBSettings


class YngwieCWToneExperiment(FPGAExperiment):
    chan = IntParameter(0)
    amp = FloatParameter(0.01)
    card = IntParameter(0)

    def sequence(self):
        constant_pulse((self.card, self.chan), 1e3, amp=self.amp,
                       phase=np.pi*0.0)

# class YngwieCWToneExperiment

class SSBPairCalibration(object):
    plot = True
    plot_vmin = -90 # 600
    plot_vmax = -50 # 2000
    plot_suffix = ''

    yng = None
    awg = None
    awg_loader = None
    spec = None
    LO = None
    vspecLO = None
    chan = None
    delay = 0.1
    setup_delay = 1.0
    if_period = -20
    target_frequency = 5e9 # just a placeholder -- will be changed dynamically


    def __init__(self, name):
        self.name = name
        
        
    def get_power(self, **kwargs):
        typ = self.spec.get_type()
        if typ == 'SA124B':
            return self.sh_get_power(**kwargs)
        elif typ == 'HP_E4407B':
            return self.hp_get_power(**kwargs)
        else:
            raise Exception('Spectrum analyzer {} not a known type.'.format(typ,))

    def set_target_freq(self, freq, **kwargs):
        typ = self.spec.get_type()
        if typ == 'SA124B':
            return self.sh_set_target_freq(freq, **kwargs)
        elif typ == 'HP_E4407B':
            return self.hp_set_target_freq(freq, **kwargs)
        else:
            raise Exception('Spectrum analyzer {} not a known type.'.format(typ,))


    def hp_get_power(self, first_run=False, new_freq=False):
        # import time
        # start = time.time()
        set_zeroif_params = True if first_run else False
        verify_freq = True if new_freq else False
        if new_freq:
            if self.LO.get_type() == 'LabBrick_RFSource': # bricks can be off by a bit
                self.spec.set_center_freq(self.target_frequency)                
                print 'LabBrick detected'
                # set_zeroif_params = True
            else:
                self.spec.set_center_freq(self.target_frequency)
        if set_zeroif_params:
            self.spec.set_span(400e3)
            self.spec.set_points(101)
#            self.spec.set_span(0)
#            self.spec.set_points(2)
            self.spec.set_if_bandwidth(5e3)
            self.spec.set_reference_level(-40)
            self.spec.write('INIT:CONT 0')
        self.spec.single_meas()
        # pwr = np.array(self.spec.do_get_data())[1].mean()
        pwr = np.array(self.spec.do_get_yaxes()).max()
        if set_zeroif_params:
            print "Setting HP mode to zero-IF (two points min)."
        if verify_freq:
            print "Verifying HP zero-IF frequency."
            actual_freq = self.spec.get_center_freq()
            diff = np.abs(self.target_frequency-actual_freq)
            if diff > 100e3:
                raise Exception('Frequency inaccurate')
        # end = time.time()
        # print "%.3fs to execute single frequency acquisition."%(end - start,)
        return pwr
    
    def hp_set_target_freq(self, frequency):
        self.target_frequency = frequency
        self.hp_get_power(new_freq=True)


    def sh_get_power(self, first_run=False, new_freq=False):
        # import time
        # start = time.time()
        set_zeroif_params = True if first_run else False
        verify_freq = True if new_freq else False
        pwr = self.spec.read_single_freq(self.target_frequency, set_zeroif_params=set_zeroif_params, verify_freq=verify_freq)
        if set_zeroif_params:
            print "Setting SignalHound mode to zero-IF (kludged)."
        if verify_freq:
            print "Verifying SignalHound zero-IF frequency."
        # end = time.time()
        # print "%.3fs to execute single frequency acquisition."%(end - start,)
        return pwr

    def sh_set_target_freq(self, frequency):
        self.target_frequency = frequency
        self.sh_get_power(new_freq=True)

        
    def load_awg_sequence(self):
        awg = self.awg
        awg.clear_sequence()
        ps = [pulselib.Square(w=1000, a=0, chan=i) for i in range(1,5)]
        s = sequencer.Sequence(sequencer.Combined(ps))
        s = sequencer.Sequencer(s)
        for i in range(1,5):
            s.add_required_channel(i)
        seqs = s.render()
        self.awg_loader.load(seqs)
        awg.run()

    def setup(self, amp=None, verbose=True):
        # self.spec.set_rfsource(self.vspecLO.get_name())
        # self.spec.set_rf_on(True)
        self.LO.set_rf_on(True)
        self.play_SSB_tone(amp=amp)
        self.get_power(first_run=True)
        self.verbose = verbose
        time.sleep(self.setup_delay)


    def finish(self):
        typ = self.spec.get_type()
        if typ == 'HP_E4407B':
            self.spec.write('INIT:CONT 1')
            self.spec.set_span(100e6)
            self.spec.set_points(401)
        else:
            pass


    def minimize(self, func, (val1, val2), (vrange1, vrange2), n_it=4, n_eval=13):
        m = minimizer.Minimizer(func, n_it=n_it, n_eval=n_eval, verbose=True)
        m.add_parameter(minimizer.Parameter('v1', value=val1, vrange=vrange1))
        m.add_parameter(minimizer.Parameter('v2', value=val2, vrange=vrange2))
        m.minimize()
        return m.params


    def _sideband_func(self, params):
        self.set_phase_params(params['v1'].value, params['v2'].value)
        time.sleep(self.delay)
        pwr = self.get_power()
        if self.verbose:
            print 'Power = {:.2f}'.format(pwr)

        if self.plot:
            c = mpltools.get_color(pwr, vmin=self.plot_vmin, vmax=self.plot_vmax)
            self.ax_sideband.scatter([params['v1'].value], [params['v2'].value],
                                     s=20, c=c, alpha=0.4, lw=0.5)
            mpltools.process_events(self.fig_sideband, activate_window=False)

        return pwr


    def minimize_sideband(self, (vrange1, vrange2), **kw):
        if self.plot:
            self.fig_sideband = plt.figure()
            self.ax_sideband = self.fig_sideband.add_subplot(111)
            self.ax_sideband.set_title("{}: sideband".format(self.name+self.plot_suffix))

        v1, v2 = self.get_phase_params()
        self.set_target_freq(self.LO.get_frequency() - self.SSBfrq)
        params = self.minimize(self._sideband_func, (v1, v2), (vrange1, vrange2), **kw)
        v1 = params['v1'].value
        v2 = params['v2'].value

        if self.plot:
            minpwr = self.get_power()
            self.ax_sideband.text(0.02, 0.98, "Min pwr = {}\nat ({:.3f},{:.3f})".format(minpwr, v1, v2),
                                  size='small', ha='left', va='top',
                                  transform=self.ax_sideband.transAxes)

            self.ax_sideband.axvline(v1, lw=0.5, dashes=[1,1])
            self.ax_sideband.axhline(v2, lw=0.5, dashes=[1,1])
            mpltools.process_events(self.fig_sideband)

        return params['v1'].value, params['v2'].value


    def _leakage_func(self, params):
        self.set_DC_offsets(params['v1'].value, params['v2'].value)
        time.sleep(self.delay)
        pwr = self.get_power()
        if self.verbose:
            print 'Power = {:.2f}'.format(pwr)

        if self.plot:
            c = mpltools.get_color(pwr, vmin=self.plot_vmin, vmax=self.plot_vmax)
            self.ax_leakage.scatter([params['v1'].value], [params['v2'].value],
                                     s=20, c=c, alpha=0.4, lw=0.5)
            mpltools.process_events(self.fig_leakage, activate_window=False)

        return pwr


    def minimize_leakage(self, (vrange1, vrange2), **kw):
        if self.plot:
            self.fig_leakage = plt.figure()
            self.ax_leakage = self.fig_leakage.add_subplot(111)
            self.ax_leakage.set_title("{}: leakage".format(self.name+self.plot_suffix))

        v1, v2 = self.get_DC_offsets()
        self.set_target_freq(self.LO.get_frequency())
        params = self.minimize(self._leakage_func, (v1, v2), (vrange1, vrange2), **kw)
        v1 = params['v1'].value
        v2 = params['v2'].value

        if self.plot:
            minpwr = self.get_power()
            self.ax_leakage.text(0.02, 0.98, "Min pwr = {}\nat ({:.3f},{:.3f})".format(minpwr, v1, v2),
                                  size='small', ha='left', va='top',
                                  transform=self.ax_leakage.transAxes)

            self.ax_leakage.axvline(v1, lw=0.5, dashes=[1,1])
            self.ax_leakage.axhline(v2, lw=0.5, dashes=[1,1])
            mpltools.process_events(self.fig_leakage)

        return params['v1'].value, params['v2'].value

# class SSBPairCalibration


class YngwieSSBPairCalibration(SSBPairCalibration):
    card = 0
    chan = 0

    def __init__(self, name):
        super(YngwieSSBPairCalibration, self).__init__(name)

    def setup(self, amp=None, **kwargs):
        if self.yng.get_nmodes()==4:
            self.chan2 = int(self.chan) # 0, 1, 2, or 3
        else:
            self.chan2 = int(self.chan>0) # 0 or 1
        self.SSBfrq = self.yng.get('ssbfreq{}'.format(self.chan2)) # ssbfreq 0..3
        if self.awg is not None:        
            self.load_awg_sequence()
        super(YngwieSSBPairCalibration, self).setup(amp=amp, **kwargs)

    def play_SSB_tone(self, amp=None):
        self.yng.stop()
        self.yng.set_unlimited(True)
        self.yng.update_modes()
        if self.yng.get_nmodes()==4:
            self.chan2 = int(self.chan) # 0, 1, 2, or 3
        else:
            self.chan2 = int(self.chan>0) # 0 or 1

        exp = YngwieCWToneExperiment(self.name)
        exp.chan = self.chan
        exp.card = self.card
        exp.amp = 0.2 if amp is None else amp
        exp.run()

    def get_DC_offsets(self):
        f = getattr(self.yng, 'get_offset{}'.format(int(self.chan>0))) # offset in 0,1 only
        return f()

    def set_DC_offsets(self, Iof, Qof):
        if self.verbose:
            print 'setting DC offsets: {}, {};'.format(int(Iof), int(Qof)),
        f = getattr(self.yng, 'set_offset{}'.format(int(self.chan>0))) # offset in 0,1 only
        f([int(Iof),int(Qof)])

    def get_phase_params(self):
        theta0 = self.yng.get('ssbtheta{}'.format(self.chan2)) # ssbtheta, ssbratio 0..3
        ratio0 = self.yng.get('ssbratio{}'.format(self.chan2))
        return theta0, ratio0

    def set_phase_params(self, theta, ratio):
        if self.verbose:
            print 'setting SSB angle and ratio: {:.3f}, {:.3f};'.format(theta, ratio),
        self.yng.set('ssbtheta{}'.format(self.chan2), theta)
        self.yng.set('ssbratio{}'.format(self.chan2), ratio)
        self.yng.update_modes()

# class YngwieSSBPairCalibration


class AWGSSBPairCalibration(SSBPairCalibration):
    Qskew = 0
    chan = (1, 2)
    delay = 0.3
    pulseamp = 1

    def __init__(self, name, cal=None):
        super(AWGSSBPairCalibration, self).__init__(name)
        self.cal = cal

        self.if_period = -20
        self.dphi = 0
        self.amps = (1, 1)

        if self.cal is not None:
            self.if_period = cal.if_period
            self.dphi = cal.dphi
            self.amps = cal.Iamp, cal.Qamp
            print 'amps:', self.amps
            
    def run_yngwie(self):
        self.yng.stop()
        self.yng.set_unlimited(True)
        self.yng.update_modes()

        exp = YngwieCWToneExperiment(self.name)
        exp.card = 0
        exp.chan = 2
        exp.amp = 0
        exp.run()


    def load_awg_sequence(self):
        awg = self.awg
        awg.clear_sequence()

        awg.set('ch{:d}_skew'.format(self.chan[1]), 0.0)
        awg.set('ch{:d}_amplitude'.format(self.chan[1]), self.AWGQamp0)

        # remember the opposite sign convention in the SSB setup of the sequencer
        print 'amps:', self.amps
        ssb = sequencer.SSB(if_period=self.if_period,
                            chans=('hay', 'boo'), dphi=self.dphi,
                            amps=self.amps,
                            outchans=self.chan)

        pulse = pulselib.AmplitudeRotation(base=pulselib.Square,
                                           w=100*abs(self.if_period),
                                           pi_amp=self.pulseamp,
                                           chans=('hay', 'boo'),
                                           pad=0)

        pulse = sequencer.Combined(pulse(alpha=np.pi, phase=0))
        s = sequencer.Sequence()
        s.append(pulse)

        # load
        s = sequencer.Sequencer(s)
        s.add_ssb(ssb)
        for i in range(1,5):
            s.add_required_channel(i)

        seqs = s.render()
        self.awg_loader.load(seqs)
        awg.run()

    def setup(self):
        self.run_yngwie()
        self.AWGIamp = self.awg.get('ch{:d}_amplitude'.format(self.chan[0]))
        self.AWGQamp = self.awg.get('ch{:d}_amplitude'.format(self.chan[1]))
        self.AWGQamp0 = self.awg.get('ch{:d}_amplitude'.format(self.chan[1]))
        self.Iof = self.awg.get('ch{:d}_offset'.format(self.chan[0]))
        self.Qof = self.awg.get('ch{:d}_offset'.format(self.chan[0]))

        self.SSBfrq = 1e9/self.if_period
        super(AWGSSBPairCalibration, self).setup()

    def play_SSB_tone(self, amp=None):
        self.load_awg_sequence()

    def finish(self):
        super(AWGSSBPairCalibration, self).finish()

        dphi = - self.Qskew * 1e-3 / self.if_period * 2 * np.pi
        awgampratio = self.AWGQamp/self.AWGIamp
        cfgampratio = self.amps[1]/self.amps[0]
        ampratio = awgampratio * cfgampratio

        if ampratio <= 1:
            Qamp = ampratio
            Iamp = 1.
        else:
            Qamp = 1.
            Iamp = 1./ampratio

        self.dphi += dphi
        self.Iamp = Iamp
        self.Qamp = Qamp
        self.amps = (self.Iamp, self.Qamp)
        self.Qskew = 0

        if self.cal is not None:
            self.cal.dphi = self.dphi
            self.cal.Iamp = self.Iamp
            self.cal.Qamp = self.Qamp
            self.cal.Iof = self.Iof
            self.cal.Qof = self.Qof
            self.cal.save_params()

        self.play_SSB_tone()

    def get_DC_offsets(self):
        Iof = self.awg.get('ch{:d}_offset'.format(self.chan[0]))
        Qof = self.awg.get('ch{:d}_offset'.format(self.chan[1]))
        return Iof, Qof

    def set_DC_offsets(self, Iof, Qof):
        print 'setting DC offsets: {:.3f}, {:.3f};'.format(Iof, Qof),
        self.Iof = Iof
        self.Qof = Qof
        self.awg.set('ch{:d}_offset'.format(self.chan[0]), Iof)
        self.awg.set('ch{:d}_offset'.format(self.chan[1]), Qof)

    def get_phase_params(self):
        Qskew = self.Qskew
        # Qskew = - self.cal.dphi * 1e3 * self.if_period / (2*np.pi)
        Qamp  = self.awg.get('ch{:d}_amplitude'.format(self.chan[1]))
        return Qskew, Qamp

    def set_phase_params(self, Qskew, Qamp):
        print 'setting Q skew, amp: {:.3f}, {:.3f};'.format(Qskew, Qamp),
        self.Qskew = Qskew
        self.AWGQamp = Qamp
        self.awg.set('ch{:d}_skew'.format(self.chan[1]), Qskew)
        self.awg.set('ch{:d}_amplitude'.format(self.chan[1]), Qamp)

# class AWGSSBPairCalibration


class AWGSSBPairWidebandCalibration(AWGSSBPairCalibration):
    ssb_frqs = None

    def __init__(self, name, cal_fn):
        super(AWGSSBPairWidebandCalibration, self).__init__(name, cal=None)
        self.cal_grp = name
        self.cal_fn = cal_fn

    def run_calibration(self, ranges, n_it=3, n_eval=11):
        self.cal_table = np.zeros(len(self.ssb_frqs),
                                  dtype=[('IF period (ns)', 'f4'),
                                           ('SSB frq (MHz)', 'f4'),
                                           ('I amp', 'f4'),
                                           ('Q amp', 'f4'),
                                           ('phase', 'f4')])
        self.cal_table.fill(np.nan)

        self.if_periods = (1e9/self.ssb_frqs).astype(int)
        self.ssb_frqs = 1e9/self.if_periods
        for i, ifp in enumerate(self.if_periods):
            self.plot_suffix = " {:.1f}".format(self.ssb_frqs[i])
            self.if_period = ifp
            self.setup()
            self.minimize_sideband(ranges, n_it=n_it, n_eval=n_eval)
            self.finish()
            self.cal_table[i] = (ifp, self.ssb_frqs[i]*1e-6, self.amps[0], self.amps[1], self.dphi)

        h5tools.save_msmt_data(self.cal_fn, grpn=self.cal_grp, datetimesubgrp=True,
                               datasets={'calibration_table' : self.cal_table})

# class AWGSSBPairWidebandCalibration


# def get_AWG_wideband_calibration_table(name, cal_fn):
#     with h5py.File(cal_fn, 'r') as f:
#         daygrpn = f[name].keys()[-1]
#         timegrpn = f[name][daygrpn].keys()[-1]
#         grpn = '{}/{}/{}'.format(name, daygrpn, timegrpn)
#         return f[grpn]['calibration_table'].value

# def get_AWG_wideband_calibration_table


# def get_AWG_wideband_calibration(name, cal_fn, ssb_frq_MHz=50):
#     table = get_AWG_wideband_calibration_table(name, cal_fn)
#     idx = np.argmin(np.abs(table['SSB frq (MHz)']-ssb_frq_MHz))
#     return table[idx]

# def get_AWG_wideband_calibration