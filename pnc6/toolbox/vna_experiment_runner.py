import time
import os
import copy
import logging
from pprint import pprint
import posixpath

import localconfig

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from toolbox.plotting import mpltools
import h5py
import h5tools
import sweeper

try:
    from qtpy import API
    if API == 'pyqt5':
        from qtpy.QtWidgets import QApplication
        from qtpy.QtGui import QPixmap
    else:
        from PyQt4.QtGui import QApplication, QPixmap
except:
    pass


MEASUREMENTS = ['S11', 'S21', 'S12', 'S22']

DATADIR = os.path.join(localconfig.datadir, 'vna')

OVERHEAD = 1.25 #overhead on waiting for data
ms = 1e3

#parameters to always save from the vna
SAVEPARAMS = ['if_bandwidth', 'power', 
                'average_factor', 'averaging_state',
                'center_freq', 'span', 'points', 
                'electrical_delay', 'correction', 'phase_offset',
                ]
FMT = 'PLOG'


def getset_vna_param(vna, getset, name, *args):
    fun_name = '{}_{}'.format(getset, name)
    if hasattr(vna, fun_name):
        fun = getattr(vna, fun_name)
        ret = fun(*args)
    else:
        raise NameError('VNA has no attribute {}'.format(fun_name))
    return ret

def get_vna_param(vna, name):
    ret = getset_vna_param(vna, 'get', name)
    if isinstance(ret, unicode):
        ret = np.string_(ret)
    return ret

def set_vna_param(vna, name, val):
    ret = getset_vna_param(vna, 'set', name, val)
    return ret


class VNAExperiment(object):

    def __init__(self, vna, fridge=None,
                plot=True, printlog=True,
                atten=0, **kwargs):
    
        # atten = e.g. -10 dB :: negative if attenuating; enter powers you actually want AFTER attenuation (not VNA output)
        self.vna = vna
        self.ys = None
        self.plot = plot
        self.printlog = printlog
        self.atten = atten

        self.fmt = kwargs.pop('fmt', FMT)
        
        self.datadir = kwargs.pop('datadir', DATADIR)
        if not os.path.isdir(self.datadir):
            print('Making data directory {}'.format(self.datadir))
            os.makedirs(self.datadir)
        self.setup_file()

        self.saveparams = SAVEPARAMS
        self.vna_params = self.vna.get_parameter_names()
        self.measurements_done = []


    def setup_file(self):
        # find and create h5 file (if necessary)
        date = time.strftime('%Y%m%d')
        self.timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        self.fn = os.path.join(self.datadir, '{}.h5'.format(date))

        '''
        get next group name
        this only gets run once so multiple realizations 
        of the experiment (with, say, different measurmenets)
        will be saved in the same group
        '''
        with h5py.File(self.fn, 'a') as f:
            grps = f.keys()
            grps_n = [int(g) for g in grps]
            if len(grps_n) > 0:
                last = max(grps_n)
                self.gn = str(last+1)
            else:
                self.gn = '0'
        self.group_created = False
            

    def setup_vna(self, **kwargs):
        for k, v in kwargs.items():
            # set things like power, IF BW, span, etc
            # if not supplied, they will be unset, keeping whatever values they were before
            if k in self.vna_params:
                set_vna_param(self.vna, k, v)

        self.vna.do_enable_averaging()
        self.vna.set_average_factor(self.average_factor)
        #get some info
        self.xs = self.vna.do_get_xaxis()
        # a few useful vna quantities
        params = ['points', 'sweep_time', 'power']
        for param_name in params:
            val = get_vna_param(self.vna, param_name)
            setattr(self, param_name, val)


    def run(self, **kwargs):
        self.average_factor = kwargs.pop('average_factor', 1)
        
        self.measurement = kwargs.get('measurement', get_vna_param(self.vna, 'measurement'))
        if self.measurement not in MEASUREMENTS:
            raise NameError('{} is not a recognized VNA measurement'.format(self.measurement))
        elif self.measurement in self.measurements_done:
            raise NameError('Measurment {} already saved for this experiment'.format(self.measurement))

        if self.printlog and kwargs != {}:
            print('Updating vna settings')
        #set up the VNA
        self.setup_vna(**kwargs)

        if self.printlog:
            print("Getting timing information and preparing trigger...")
        time.sleep(0.1)

        tot_time = self.average_factor * self.sweep_time
        
        if self.printlog:
            print('Total duration of this experiment will be {:.2f} minutes.'.format(tot_time/60))
            
        if self.plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        try:
            # an approximation based on 10ms/pt for transfer, 25% overhead on waiting/capture
            timeout = (self.sweep_time * self.average_factor * OVERHEAD * ms) + self.points*10 # ms
            '''
            this is the objectsharwer timeout,
            and the VNA timeout if short enough
            otherwise we want to know that the VNA is not responding sooner, 
            and so we trigger every average. The VNA timeout now corresponds to the time of one scan
            '''

            timeout = np.max(np.array([timeout, 5000.])) # 5s timeout minimum
            if timeout < 30e3:   
                self.vna.set_timeout(timeout)
                self.ys = self.vna.do_get_data(fmt=self.fmt, opc=True, timeout=timeout)
            else:
                if self.printlog:
                    print("[NOTICE] Triggering each average.")
                to = np.max(np.array([self.sweep_time * OVERHEAD * ms, 5000.]))
                self.vna.set_timeout(to)

                self.ys = self.vna.do_get_data(fmt=self.fmt, opc=True, trig_each_avg=True, timeout=timeout)
         
            try:
                self.temp=self.fridge.do_get_temperature()*1000
            except:
                self.temp=-1.0
                
            try:
                if self.plot:
                    ax.plot(self.xs,self.ys[0],label='{:.2f} dBm, {:.2f} mK'.format(self.power, self.temp))
                    fig.canvas.draw()
                    # QtGui.QApplication.processEvents()
            except:
                raise
            
            self.save_data()
                
            if self.plot:
                ax.set_title(self.measurement)
                ax.set_xlabel("Frequency (Hz)")
                ax.set_ylabel("Magnitude (dB)")
                ax.legend(loc='best')

        except Exception as e:
            print "EXCEPTION", e
            raise


    def save_data(self):
        if self.printlog:
            print('Saving data')
        with h5py.File(self.fn) as f:
            if not self.group_created:
                g = f.create_group(self.gn)
                self.group_created = True
            else:
                g = f[self.gn]
            g.attrs.create('attenuation', self.atten)

            sg = g.create_group(self.measurement)
            for param_name in self.saveparams:
                param_val = get_vna_param(self.vna, param_name)
                sg.attrs.create(param_name, param_val)
            sg.attrs.create('timestamp', self.timestamp)
            ds_freqs = sg.create_dataset('frequencies', data=self.xs)
            ds_mag = sg.create_dataset('magnitude', data=self.ys[0])
            ds_phase = sg.create_dataset('phase', data=self.ys[1])

            if self.temp>0:
                sg.attrs.create('fridge_temperature', self.temp)
            
        if self.printlog:
            print('Data saved at {}_{}_{}'.format(self.fn, self.gn, self.measurement))
        self.measurements_done.append(self.measurement)


class VNAExperimentRunner(object):
    debug = False

    def __init__(self, vna, measurements=None, additional_settings={}, plot=False, verbose=False, **expkwargs):
        self.vna = vna
        #measurements is which data the VNA will acquire
        if measurements is None:
            measurements = [get_vna_param(self.vna, 'measurement'), ]
        if isinstance(measurements, str):
            #just one measurement
            measurements = [measurements,]
        for measurement in measurements:
            if measurement not in MEASUREMENTS:
                raise ValueError('{} is not a recognized VNA measurement'.format(measurement))

        self.measurements = measurements
        self.n_measurements = len(self.measurements)
        
        # additional settings is a dict of param names and lists
        # the lists should be the same length as measurements
        # these are the params that are set differently for each measurement (span, average_factor, etc)
        for additional_k, additional_vs in additional_settings.items():
            if len(additional_vs) != self.n_measurements:
                raise ValueError('Additional setting list {} must have same length as measurements'.format(additional_k))
        self.additional_settings = additional_settings

        self.plot = plot
        self.name = 'VNAExperiment'
        
        self.expkwargs = expkwargs
        self.verbose = verbose
        self.exp = VNAExperiment(self.vna, plot=self.plot, printlog=self.verbose)


    def run(self, **kw):
        #add any kwargs
        for k in kw:
            self.expkwargs[k] = kw[k]

        for im, measurement in enumerate(self.measurements):
            self.expkwargs['measurement'] = measurement
            for additional_k, additional_vs in self.additional_settings.items():
                self.expkwargs[additional_k] = additional_vs[im]

            if self.debug:
                print('Running experiment...')
            self.exp.run(**self.expkwargs)

            t0 = time.time()
            tplot = time.time()

            if self.debug:
                pprint(self.expkwargs)


class VNAExperimentSweep(sweeper.Sweeper):
    def __init__(self, vna, name, verbose=False):
        super(VNAExperimentSweep, self).__init__(name)

        self.vna = vna
        self.expt_filenames = []
        self.expt_grps = []
        self.expkwargs = {}
        self.verbose = verbose


    def run_sweep(self, **kw):
        # allow static experiment kwargs to be supplied here
        # for k in kw:
        #     self.expkwargs[k] = kw[k]
        self.expkwargs.update(kw)

        super(VNAExperimentSweep, self).run_sweep()


    def update_expkwargs(self, cur_sweep_vals):
        #update the experiment kwargs
        this_kwargs = {}
        for i_p, param_name in enumerate(self.sweep_axes_names):
            param_val = cur_sweep_vals[i_p]
            this_kwargs[param_name] = param_val
        self.expkwargs.update(this_kwargs)


    def msmt_loop(self):
        e = VNAExperimentRunner(self.vna, **self.expkwargs)
        e.run()

        fn = copy.deepcopy(e.exp.fn)
        gn = copy.deepcopy(e.exp.gn)

        self.expt_filenames.append(fn)
        self.expt_grps.append(gn)

        return ({},
                {'expt_filenames' : self.expt_filenames,
                 'expt_grps' : self.expt_grps, },
                {})


class VNAExperimentSweepData(object):
    replace_basedir = None

    def __init__(self, fn, grp):
        self.sweep = h5tools.get_content_from_file(fn, grp)
        self.expt_name = os.path.splitext(os.path.split(fn)[1])[0]
        self.expt_fns = self.sweep['attrs']['expt_filenames']
        self.expt_grps = self.sweep['attrs']['expt_grps']
        self.expt_datestamps = [os.path.split(fn)[1][:-3] \
            for fn in self.expt_fns]
        self.tag = "{}/{}".format(self.expt_name, grp)

        self.sweep_axes_names = self.sweep['attrs']['sweep_axes_names']
        if 'sweep_axes_lens' in self.sweep['attrs']:
            self.sweep_shape = tuple(self.sweep['attrs']['sweep_axes_lens'])
        else:
            self.sweep_shape = []
            for a in self.sweep['attrs']['sweep_axes_names']:
                self.sweep_shape.append(self.sweep[a]['value'].size)
            self.sweep_shape = tuple(self.sweep_shape)

        self.nsweep_axes = len(self.sweep_shape)
        self.sweep_idxs = [tuple(idx) for idx in self.sweep['completed_sweep_idxs']['value']]
        self.expt_idxs = np.arange(np.prod(self.sweep_shape)).reshape(
            self.sweep_shape)

        self.sweep_vals = [self.sweep[a]['value'] for a in self.sweep_axes_names]
        

    def get_result(self, idxs):
        idx = self.expt_idxs[idxs]
        if idx > len(self.sweep_idxs)-1:
            return None
        fn, grp = self.expt_fns[idx], self.expt_grps[idx]
        if self.replace_basedir is not None:
            fn = posixpath.join(*fn.split('\\'))
            relfnpath = os.path.relpath(fn, self.replace_basedir[0])
            fn = os.path.join(self.replace_basedir[1], relfnpath)
        return self.datasets_from_file(fn, grp)
    
    
    def get_results(self):
        # return a dictionary of n-dim arrays of all vna data
        first_idx = self.sweep_idxs[0]
        first_result = self.get_result(first_idx)
        self.msmts = first_result.keys()
        msmt_lens = [len(first_result[msmt]['fs']) for msmt in self.msmts]
        data_shapes = tuple(np.concatenate((self.sweep_shape, [msmt_len,])) for msmt_len in msmt_lens)
        results = {msmt : np.zeros(data_shapes[im], dtype=np.complex128) for im, msmt in enumerate(self.msmts)}
        fs = {}

        for msmt in self.msmts:
            fs[msmt] = first_result[msmt]['fs']
            for idxs in self.sweep_idxs:
                result = self.get_result(idxs)
                results[msmt][idxs] = 10**(result[msmt]['ms']/20) * np.exp(1j*(np.pi/180)*result[msmt]['ps'])
        return results, fs  


    def get_sweep_vals(self, idxs):
        ret = []
        for i,a in enumerate(self.sweep_axes_names):
            ret.append(self.sweep[a]['value'][idxs[i]])
        return ret


    def get_sweep_vals_str(self, idxs, numfmt='', sep='; '):
        ret = ''
        for i,a in enumerate(self.sweep_axes_names):
            ret += ('{} = {' + numfmt + '} ({}/{})').format(
                a, self.sweep[a]['value'][idxs[i]], idxs[i]+1, self.sweep_shape[i])

            if i < (len(self.sweep_axes_names)-1):
                ret += sep

        return ret

    
    @staticmethod
    def datasets_from_file(fn, gn):
        datasets = {}
        with h5py.File(fn) as f:
            g = f[gn]
            dataset_names = g.keys()
            for dataset_name in dataset_names:
                dataset = {}
                sg = g[dataset_name]
                the_attrs = sg.attrs
                attrs = {}
                for a in the_attrs:
                    attrs[a] = the_attrs.get(a)
                dataset['attrs'] = attrs
                dataset['fs'] = sg['frequencies'][:]
                dataset['ms'] = sg['magnitude'][:]
                dataset['ps'] = sg['phase'][:]
                datasets[dataset_name] = dataset
        return datasets