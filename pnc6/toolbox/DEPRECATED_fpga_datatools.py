import os
import time
import numpy as np
import numpy.ma as ma

import h5tools
import localconfig

datadir = localconfig.datadir

### tools for data generated with fpga_lib
def get_experiment(exp_name, datestamp, group):
    fn = os.path.join(datadir, "exp", exp_name, "archive", "{}.h5".format(datestamp))
    content = h5tools.get_content_from_file(fn, group)
    content['time'] = time.strftime("%H%M%S", time.strptime(content['attrs']['timestamp'], 
                                                            "%Y-%m-%d %H:%M:%S"))
    return content


def get_threshold(exp, readout_name='readout', threshold='thresh0'):
    return exp['calibrations'][readout_name]['attrs'][threshold]


def condition_on_readouts(condition_data, criteria, threshold=0):
    """
    Create an array that contains elements that are True if all conditions are met,
    or False if any condition is not met. 

    Conditions are specified as follows:
        condition_data : list
            list of arrays of same shape with (complex) readout SE values
        criteria : list
            list of criteria (one per condition data); 
            either 'left' or 'right' for each element.
            'left' corresponds to success if the value is smaller than the
            threshold.
        thresholds = 0 : list or float
            if not a list, use the same threshold for all criteria.

    Note that all comparison is only done with the real parts right now.
    """
    success = np.ones(condition_data[0].shape, dtype=np.bool)
    for ic, c in enumerate(condition_data):
        if type(threshold) == list:
            th = threshold[ic]
        else:
            th = threshold
            
        if criteria[ic] == 'left':
            success = success & (c.real < th)
        else:
            success = success & (c.real >= th)
    
    return success


def get_experiment_data(experiment, data_names, threshold=None, avg_axis=0,
    condition_data_names=[], condition_criteria=[], condition_threshold=0, real_data=True):

    ret = {}
    success = None
    nsuccess = experiment[data_names[0]]['data']['value'].size
    nall = experiment[data_names[0]]['data']['value'].size
    ret['shots'] = nall
    ret['condition'] = {}
    ret['condition']['conditioned'] = False
    
    if len(condition_data_names) > 0:
        try:
            condition_data = [experiment[cdn]['data']['value'] for cdn in condition_data_names]
            success = condition_on_readouts(condition_data, condition_criteria, condition_threshold)
            nsuccess = len(np.where(success)[0])
            nall = success.size
            ret['condition']['conditioned'] = True

        except KeyError:
            success = None
            nsuccess = experiment[data_names[0]]['data']['value'].size

    ret['condition']['successful shots'] = nsuccess
    ret['condition']['success rate'] = float(nsuccess)/nall

    for idn, dn in enumerate(data_names):
        data = experiment[dn]['data']['value'].copy()
        naxes = len(data.shape)
        if threshold is not None:
            if type(threshold) == list:
                th = threshold[idn]
            else:
                th = threshold

            iszero = data.real < th
            data[iszero] = 0.
            data[np.logical_not(iszero)] = 1.

            if success is not None:
                data = ma.masked_array(data, mask=np.logical_not(success))

        if avg_axis is not None:
            if real_data:
                data = data.mean(axis=avg_axis).real
            else:
                data = data.mean(axis=avg_axis)

        ret[dn] = {}
        ret[dn]['value'] = data
        ret[dn]['axes'] = []
        for i in range(naxes):
            ax = {}
            if avg_axis != i:
                ax['value'] = experiment[dn]['ax_data_{}'.format(i)]['value'].copy()
                ax['name'] = experiment[dn]['data']['attrs']['DIMENSION_LABELS'][i]
                ret[dn]['axes'].append(ax)

    return ret


### tools for use with the sweeper
def get_sweep(exp_name, datestamp, timestamp):
    fn = os.path.join(datadir, "{}.h5".format(exp_name))
    grp = '{}/{}'.format(datestamp, timestamp)
    return h5tools.get_content_from_file(fn, grp)


class SweepData(object):
    def __init__(self, exp_name, datestamp, timestamp):
        self.sweep = get_sweep(exp_name, datestamp, timestamp)
        self.expt_name = exp_name
        self.expt_fns = self.sweep['attrs']['expt_filenames']
        self.expt_grps = self.sweep['attrs']['expt_grps']
        self.expt_datestamps = [os.path.split(fn)[1][:-3] for fn in self.expt_fns]
        self.tag = "{}/{}/{}".format(exp_name, datestamp, timestamp)
        
        if 'sweep_axes_lens' in self.sweep['attrs']:
            self.sweep_shape = tuple(self.sweep['attrs']['sweep_axes_lens'])
        else:
            self.sweep_shape = []
            for a in self.sweep['attrs']['sweep_axes_names']:
                self.sweep_shape.append(self.sweep[a]['value'].size)
            self.sweep_shape = tuple(self.sweep_shape)

        self.nsweep_axes = len(self.sweep_shape)
        self.sweep_idxs = self.sweep['completed_sweep_idxs']['value']
        self.expt_idxs = np.arange(np.prod(self.sweep_shape)).reshape(self.sweep_shape)

        # self.completed_slice = np.array([])
        # for i in range(self.nsweep_axes)[::-1]:
        #     if maxes[i]
    
    def get_attrs(self, attrs):
        ret = {}

        for ds, gn, idxs in zip(self.expt_datestamps, self.expt_grps, self.sweep_idxs):
            exp = get_experiment(self.expt_name, ds, gn)
            
            for a in attrs:
                if a not in ret:
                    ret[a] = np.empty(self.sweep_shape, dtype=exp['attrs'][a].dtype)
                    ret[a].fill(np.nan)
                
                ret[a][tuple(idxs)] = exp['attrs'][a]

        return ret


    def get_data(self, data_names, **kw):
        ret = {}
        
        for ds, gn, idxs in zip(self.expt_datestamps, self.expt_grps, self.sweep_idxs):
            exp = get_experiment(self.expt_name, ds, gn)
            data = get_experiment_data(exp, data_names, **kw)

            for idatan, datan in enumerate(data_names):
                if datan not in ret:
                    datashape = data[datan]['value'].shape
                    sweepshape = tuple(list(self.sweep_shape) + list(datashape))
                    ret[datan] = {'value' : np.empty(sweepshape, dtype=data[datan]['value'].dtype)}
                    ret[datan]['value'].fill(np.nan)
                    ret[datan]['axes'] = data[datan]['axes']
                ret[datan]['value'][tuple(idxs)] = data[datan]['value']

        return ret
