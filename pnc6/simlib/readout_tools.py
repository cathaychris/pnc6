import numpy as np
from scipy.integrate import quad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pnc6.toolbox.plotting import mpltools


twopi = 2*np.pi


def square_pulse(t0, l):
    return lambda t: (np.sign(t-t0) +1)/2 -  (np.sign((t-t0-l)) +1)/2
    

def refl_filter(freqs, k, d):
    return (k/2+1j*(freqs-d))/(k/2-1j*(freqs-d))


def trans_filter(freqs, k, d):
    return k/2/(k/2-1j*(freqs-d))


def apply_filters_ge(incoming, kappa, chi, delta_rel, dt=0.020, filt=trans_filter):
    ft = np.fft.fft(incoming)
    freqs = np.fft.fftfreq(len(incoming), dt)

    delta_g = delta_rel*chi
    delta_e = -(1-delta_rel)*chi
    
    g_filt = filt(-freqs, kappa, delta_g)
    gft_filt = ft*g_filt
    e_filt = filt(-freqs, kappa, delta_e)
    eft_filt = ft*e_filt
    
    traj_g = np.fft.ifft(gft_filt)
    traj_e = np.fft.ifft(eft_filt)
    
    return np.array([traj_g, traj_e])

def get_env(traj_g, traj_e):
    env = np.conjugate(traj_g - traj_e)
    return env


def plot_trajectories(tlist, trajs, datas, shots_axis=1, labels = ['g','e']):
    avg_datas = datas.mean(axis=shots_axis)
    ntrajs = len(trajs)
    fig = plt.figure(figsize=(ntrajs*2, 3))
    gs = gridspec.GridSpec(2, ntrajs, height_ratios=[2, 1])
    axes = []
    for ix in range(ntrajs):
        ax1kw, ax2kw = {}, {}
        if ix > 0:
            ax1kw['sharex'] = axes[0][0]
            ax1kw['sharey'] = axes[0][0]
            ax2kw['sharey'] = axes[0][1]

        ax1 = plt.subplot(gs[0, ix], **ax1kw)
        ax2 = plt.subplot(gs[1, ix], sharex=ax1, **ax2kw)
        axes.append( ( ax1, ax2 ) )

    for ix, (ax, rax) in enumerate(axes):
        traj = trajs[ix]
        avg_data = avg_datas[ix]
        ax.plot(tlist, traj.real,'k')
        ax.plot(tlist, traj.imag,'k')
        ax.plot(tlist, avg_data.real)
        ax.plot(tlist, avg_data.imag)
        resid = avg_data - traj
        rax.plot(tlist, resid.real)
        rax.plot(tlist, resid.imag)
        plt.setp(ax.get_xticklabels(), visible=False)
        if ix > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.setp(rax.get_yticklabels(), visible=False)
        ax.set_title(labels[ix])

    return fig


# def plot_shots(data_g, data_e, shots_axis=0):
#     avg_g, avg_e = data_g.mean(axis=shots_axis), data_e.mean(axis=shots_axis)

#     datas = [data_g, data_e, np.vstack((data_g, data_e))]
    
#     fig, axes = plt.subplots(1,3,figsize=(6,3), sharex=True,sharey=True)
#     for ix, ax in enumerate(axes):
#         integrated = datas[ix].sum(axis=1)
#         histo, x, y = np.histogram2d(integrated.real, integrated.imag, bins=40)
#         ax.pcolormesh(x, y, histo.transpose())

#     #apply envelope
#     histos = []
#     integrated_shots = []
#     env = get_env(avg_g, avg_e)
#     fig, axes = plt.subplots(1,3,figsize=(6,3), sharex=True,sharey=True)
#     for ix, ax in enumerate(axes):
#         enveloped = datas[ix]*env
#         integrated = enveloped.sum(axis=1)
#         integrated_shots.append(integrated)
#         histo, x, y = np.histogram2d(integrated.real, integrated.imag, bins=40)
#         histos.append((histo, x, y))
#         ax.pcolormesh(x, y, histo.transpose())

#     fig, ax = plt.subplots()
#     for integrated in integrated_shots:
#         ax.hist(integrated.real, normed=True, alpha=0.3, bins = 40)


def plot_shots(datas, shots_axis=1, env_map=[0,1], labels=['g','e']):
    labels.append('combined')
    ntrajs = datas.shape[0]
    avgs = datas.mean(axis=shots_axis)

    datas = list(datas)
    datas.append(np.vstack(tuple(datas)))
    
    fig, axes = plt.subplots(1, ntrajs+1, figsize=(2*(ntrajs+1), 3), sharex=True, sharey=True)
    for ix, ax in enumerate(axes):
        integrated = datas[ix].sum(axis=1)
        histo, x, y = np.histogram2d(integrated.real, integrated.imag, bins=40)
        ax.pcolormesh(x, y, histo.transpose())
        ax.set_title(labels[ix])

    #apply envelope
    histos = []
    integrated_shots = []
    env = get_env(avgs[env_map[0]], avgs[env_map[1]])
    fig, axes = plt.subplots(1, ntrajs+1, figsize=(2*(ntrajs+1), 3), sharex=True, sharey=True)
    for ix, ax in enumerate(axes):
        enveloped = datas[ix]*env
        integrated = enveloped.sum(axis=1)
        integrated_shots.append(integrated)
        histo, x, y = np.histogram2d(integrated.real, integrated.imag, bins=40)
        histos.append((histo, x, y))
        ax.pcolormesh(x, y, histo.transpose())
        ax.set_title(labels[ix])

    fig, [iax, qax] = plt.subplots(1,2, sharey=True)
    for i, integrated in enumerate(integrated_shots):
        iax.hist(integrated.real, normed=True, alpha=0.3, bins = 40)
        qax.hist(integrated.imag, normed=True, alpha=0.3, bins = 40, label=labels[i])
    iax.set_title('I')
    qax.set_title('Q')
    qax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)



