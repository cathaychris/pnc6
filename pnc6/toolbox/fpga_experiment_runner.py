import time
import os
import copy
import logging
from pprint import pprint
import posixpath

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from fpga_lib.dsl import result

from toolbox.plotting import mpltools
import h5tools
import sweeper


class ExperimentRunner(object):
    debug = True

    plot_defaults = {
        'indices' : None,
        'thresh' : True,
        'axes' : None,
    }

    def __init__(self, expclass, name, **expkwargs):
        self.do_plot = True
        self.plot_interval = 0.1
        self.fig = None
        self.fig_num = None
        self.plots = {
            'default': {
                'datasets' : ['default', ],
                }
            }
        self.note = None

        self.name = name
        self.exp = expclass(name)
        self.expkwargs = {
            'fit_threshold' : True,
            'n_blocks' : 10,
            'averages_per_block' : 100, }
        for k in expkwargs:
            self.expkwargs[k] = expkwargs[k]


    def run(self, **kw):
        for k in kw:
            self.expkwargs[k] = kw[k]

        for k in self.expkwargs:
            setattr(self.exp, k, self.expkwargs[k])
        self.exp.run()

        t0 = time.time()
        tplot = time.time()

        cur_blocks = self.exp.n_blocks_available()
        try:
            while cur_blocks < self.exp.n_blocks:
                if self.exp.n_blocks_available() > cur_blocks:
                    cur_blocks = self.exp.n_blocks_available()
                    self.exp.process_available_blocks()

                    if self.do_plot and (time.time()-tplot) > self.plot_interval:
                        try:
                            self.update_plot()
                            tplot = time.time()
                        except:
                            logging.warning('Could not plot...')

                mpltools.wait_and_update(0.1, fig=self.fig)

        except KeyboardInterrupt:
            logging.info("User interrupt")
            self.exp.process_available_blocks()
            self.exp.save_data()
            self.exp.stop_running()
            if self.do_plot:
                try:
                    self.update_plot()
                except:
                    logging.warning('Could not plot...')


        if self.fig is not None:
            mpltools.update_and_copy(self.fig)

            datafn = self.exp.last_file_name
            datagrp = self.exp.last_file_group
            p, fn = os.path.split(datafn)
            p, _ = os.path.split(p)

            figfn = '{}__{}'.format(os.path.splitext(fn)[0], datagrp)
            if self.note is not None:
                figfn += '__{}'.format(self.note)

            self.fig.savefig(os.path.join(p, 'images', "{}.png".format(figfn)))

        if self.debug:
            pprint(self.expkwargs)


    def update_plot(self):
        if self.fig is None:
            self.fig = plt.figure(num=self.fig_num)
            self.fig.canvas.manager.window.activateWindow()
            self.fig.canvas.manager.window.raise_()

        self.fig.clear()
        self.plot()
        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.9)


    def _fit_label(self, fitparams):
        lbl = ""
        for p in fitparams:
            if p[0] != '_':
                lbl += '{} = {:1.3e} +/- {:1.3e}'.format(p, fitparams[p], fitparams['_{}_err'.format(p)])
                lbl += '\n'
        return lbl


    def plot(self):
        nplots = len(self.plots)

        for ik, k in enumerate(self.plots):
            ax = self.fig.add_subplot(ik+1, 1, nplots)

            p = self.plots[k]
            for k2 in self.plot_defaults:
                if k2 not in p:
                    p[k2] = self.plot_defaults[k2]

            for idset, dset in enumerate(p['datasets']):
                r = self.exp.results[dset]
                if p['thresh']:
                    r = r.thresh_mean()
                else:
                    r = r.axis_mean()

                data = r.data.copy().real
                if p['indices'] is not None:
                    data = data[p['indices']]

                if len(data.shape) == 1:
                    if p['axes'] is None:
                        xvals = np.arange(data.size)
                        xlabel = 'x point'
                    else:
                        xvals = r.ax_data[p['axes'][0]]
                        xlabel = r.labels[p['axes'][0]]
                    ax.plot(xvals, data, 'o', label=dset)

                    ax.set_xlabel(xlabel)
                    ax.legend(loc='best', fontsize='x-small')

                elif len(data.shape) == 2:
                    im = ax.imshow(data, cmap=cm.bwr)
                    cb = self.fig.colorbar(im)

            ax.set_title(k, size='x-small')

        figtitle = "{} [{}]".format(self.exp.last_file_name,
                                    self.exp.last_file_group)
        if self.note is not None:
            figtitle += " ({})".format(self.note)

        self.fig.suptitle(figtitle, size='x-small')

class FPGAExperimentSweep(sweeper.Sweeper):
    expclass = None
    do_exp_plot = False
    exp_plots = None
    exp_fig_num = None
    # exp_plot_results = ['default']
    # exp_plot_thresholded = True
    # exp_plot_title = None

    fit_params_datasets = []

    def __init__(self, expclass, name):
        super(FPGAExperimentSweep, self).__init__(name)

        self.expclass = expclass
        self.expt_filenames = []
        self.expt_grps = []
        self.expkwargs = {}

        if self.do_plot:
            self.fig = plt.figure(self.fig_num)


    def pre_sweep_tasks(self):
        for p in self.fit_params_datasets:
            self.create_empty_dataset(p, shape=self.sweep_axes_lens,
                                      ax_data=self.sweep_axes_vals,
                                      ax_labels=self.sweep_axes_names,
                                      ylabel=p)


    def run_sweep(self, **kw):
        for k in kw:
            self.expkwargs[k] = kw[k]

        super(FPGAExperimentSweep, self).run_sweep()


    def msmt_loop(self):
        e = ExperimentRunner(self.expclass, self.name)
        e.do_plot = self.do_exp_plot
        if self.exp_plots is not None:
            e.plots = self.exp_plots
        e.note = str(self.cur_sweep_idxs).replace(' ', '')
        e.fig_num = self.exp_fig_num
        e.run(**self.expkwargs)

        fn = copy.deepcopy(e.exp.last_file_name)
        gn = copy.deepcopy(e.exp.last_file_group)

        self.expt_filenames.append(fn)
        self.expt_grps.append(gn)

        for p in self.fit_params_datasets:
            self.set_datapoint(p, self.cur_sweep_idxs, e.exp.fit_params[p])


        return ({},
                {'expt_filenames' : self.expt_filenames,
                 'expt_grps' : self.expt_grps, },
                {})


### class to easily get the data from a sweep
class FPGAExperimentSweepData(object):
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
        self.sweep_idxs = self.sweep['completed_sweep_idxs']['value']
        self.expt_idxs = np.arange(np.prod(self.sweep_shape)).reshape(
            self.sweep_shape)


    def get_result(self, idxs):
        idx = self.expt_idxs[idxs]
        if idx > len(self.sweep_idxs)-1:
            return None
        fn, grp = self.expt_fns[idx], self.expt_grps[idx]
        if self.replace_basedir is not None:
            fn = posixpath.join(*fn.split('\\'))
            relfnpath = os.path.relpath(fn, self.replace_basedir[0])
            fn = os.path.join(self.replace_basedir[1], relfnpath)
        return result.Results.create_from_file(fn, grp)


    def get_sweep_vals(self, idxs):
        ret = []
        for i,a in enumerate(self.sweep_axes_names):
            ret.append(self.sweep[a]['value'][idxs[i]])


    def get_sweep_vals_str(self, idxs, numfmt='', sep='; '):
        ret = ''
        for i,a in enumerate(self.sweep_axes_names):
            ret += ('{} = {' + numfmt + '} ({}/{})').format(
                a, self.sweep[a]['value'][idxs[i]], idxs[i]+1, self.sweep_shape[i])

            if i < (len(self.sweep_axes_names)-1):
                ret += sep

        return ret


