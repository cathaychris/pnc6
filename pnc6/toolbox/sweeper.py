import os
import time
import numpy as np
from matplotlib import pyplot as plt

from toolbox import h5tools
from toolbox.plotting import mpltools
import localconfig

import copy

class Sweeper(object):
    """
    TODO:
    * should maybe save data/images/settings like phil does
    * logger
    * plotting: more than 1 dimension? traces? error bars?
    """
    sweep_axes_names = ['']
    sweep_axes_vals = np.array([[0]])
    pre_msmt = None
    post_msmt = None
    
    def __init__(self, name):
        self.name = name
        self.filename = os.path.join(localconfig.datadir, self.name+'.h5')
        self.grpn = None
        self.data_location = None
        
        self.do_plot = False
        self.fig = None
        self.fig_num = None
        self.datasets = {}
        self.timings = []
        
    def create_empty_dataset(self, name, shape, ax_data=None, ax_labels=None, ylabel=''):
        data = np.ones(shape) * np.nan
        self.set_dataset(name, data, ax_data, ax_labels, ylabel)


    def set_dataset(self, name, data, ax_data=None, ax_labels=None, ylabel=''):
        self.datasets[name] = {}
        self.datasets[name]['data'] = data
        if ax_data is None:
            self.datasets[name]['ax_data'] = [np.arange(a) for a in data.shape]
        else:
            self.datasets[name]['ax_data'] = ax_data
        if ax_labels is None:
            self.datasets[name]['ax_labels'] = ['axis {}'.format(i) for i in range(len(data.shape))]
        else:
            self.datasets[name]['ax_labels'] = ax_labels
        self.datasets[name]['ylabel'] = ylabel
            
    
    def set_datapoint(self, name, idxs, val):           
        self.datasets[name]['data'][idxs] = val
     

    def update_plot(self):
        if self.fig is None:
            self.fig = plt.figure(num=self.fig_num)
            self.fig.canvas.manager.window.activateWindow()
            self.fig.canvas.manager.window.raise_()
        self.plot()
        self.fig.suptitle(self.data_location, size='x-small')
        self.fig.tight_layout()
        mpltools.process_events(fig=self.fig)
        
    
    def plot(self):
        self.fig.clear()
        naxes = len(self.datasets.keys())
        for i,d in enumerate(self.datasets):
            ax = self.fig.add_subplot(naxes, 1, i+1)
            ax.plot(self.datasets[d]['ax_data'][0], self.datasets[d]['data'], 'o-')
            ax.set_xlabel(self.datasets[d]['ax_labels'][0])
            ax.set_ylabel(self.datasets[d]['ylabel'])
    

    def save_data(self, datasets={}, attrs={}, dataset_attrs={}):
        self.grpn = h5tools.save_msmt_data(self.filename, self.grpn,
                                           datasets=datasets, attrs=attrs,
                                           dataset_attrs=dataset_attrs)
        self.data_location = "{}/{}".format(os.path.split(self.filename)[1], self.grpn)


    def run_next(self):
        print
        print 'Dataname: {}'.format(self.data_location)
        print
        print 'Run msmt {} of {}'.format(self.cur_sweep_pt, self.nsweep_pts)
        print 'Sweep values:'
        for idim, dim in enumerate(self.sweep_axes_names):
            print '  {} = {} ({}/{})'.format(dim,
                self.sweep_axes_vals[idim][self.cur_sweep_idxs[idim]], self.cur_sweep_idxs[idim]+1,
                self.sweep_axes_lens[idim])
        
        if self.pre_msmt is not None:
            print
            print "Run pre-msmt tasks..."
            self.pre_msmt(self, self.cur_sweep_vals)

        t0 = time.time()
        ret = self.msmt_loop()
        if self.do_plot:
            self.update_plot()
        t1 = time.time()
        elapsed_run_time = t1-t0
        self.timings.append(elapsed_run_time)

        self.completed_sweep_idxs.append(self.cur_sweep_idxs)
        self.completed_pts += 1
        self.save_data(datasets={'completed_sweep_idxs' : self.completed_sweep_idxs})

        if ret:
            datasets, attrs, dattrs = ret
            attrs['elapsed_run_times'] = self.timings
            self.save_data(datasets=datasets, attrs=attrs, dataset_attrs=dattrs)

        print "... finished data point in {} seconds.".format(elapsed_run_time)
        print
        
        if self.post_msmt is not None:
            print "Run post-msmt tasks..."
            print
            self.post_msmt(self, self.cur_sweep_vals)

    
    def pre_sweep_tasks(self):
        pass


    def run_sweep(self):
        self.sweep_axes_lens = tuple([len(a) for a in self.sweep_axes_vals])
        self.nsweep_pts = np.product(np.array(self.sweep_axes_lens))
        self.cur_sweep_idxs = tuple([0 for d in self.sweep_axes_names])
        self.cur_sweep_vals = tuple([a[0] for a in self.sweep_axes_vals])
        self.cur_sweep_pt = 1
        self.completed_sweep_idxs = []
        self.completed_pts = 0
        self.start_time = time.time()

        axesdsets = {}
        for i,n in enumerate(self.sweep_axes_names):
            axesdsets[n] = self.sweep_axes_vals[i]
        self.save_data(datasets=axesdsets, 
                       attrs={'sweep_axes_names' : self.sweep_axes_names,
                              'sweep_axes_lens' : self.sweep_axes_lens, })

        self.pre_sweep_tasks()        
        
        try:
            self.run_next()
            while self.cur_sweep_idxs != tuple([l-1 for l in self.sweep_axes_lens]):
                for idim in range(len(self.sweep_axes_names)):
                    dimidx = len(self.sweep_axes_names) - idim - 1
                    if self.cur_sweep_idxs[dimidx] < self.sweep_axes_lens[dimidx]-1:
                        _cur_idxs = list(self.cur_sweep_idxs)
                        _cur_idxs[dimidx] += 1
                        self.cur_sweep_idxs = tuple(_cur_idxs)
                        # print 'Incrementing {} to idx {}'.format(self.sweep_dims[dimidx], _cur_idxs[dimidx])
                        break
                    else:
                        _cur_idxs = list(self.cur_sweep_idxs)
                        _cur_idxs[dimidx] = 0
                        self.cur_sweep_idxs = tuple(_cur_idxs)
                        # print 'Resetting {} to idx 0'.format(self.sweep_dims[dimidx])

                self.cur_sweep_vals = tuple([copy.copy(a[idx]) for a,idx in zip(self.sweep_axes_vals, self.cur_sweep_idxs)])
                self.cur_sweep_pt += 1
                self.run_next()

        except KeyboardInterrupt:
            print "Sweep interrupted by user."

        finally:
            pass

    def msmt_loop(self):
        return None


if __name__ == '__main__':
    s = Sweeper('sweeptest')
    s.sweep_axes_names = ('axis a', 'axis b')
    s.sweep_axes_vals = (np.arange(3), np.arange(4))
    s.run_sweep()
