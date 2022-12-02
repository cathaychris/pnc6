import time
import numpy as np
import math
import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from matplotlib.colors import rgb2hex
from matplotlib.ticker import MaxNLocator

try:
    from qtpy import API
    if API == 'pyqt5':
        from qtpy.QtWidgets import QApplication
        from qtpy.QtGui import QPixmap
    else:
        from PyQt4.QtGui import QApplication, QPixmap
        # from qtpy.QtGui import QApplication, QPixmap
except:
	pass

# import colormaps
# viridis = colormaps.viridis
default_cmap = cm.viridis
cmap_div = cm.PRGn
cmap_seq = cm.viridis


# tools for figure handling
def copy(fig=None):
    if fig is None:
        canvas = plt.gcf().canvas
    else:
        canvas = fig.canvas

    canvas.draw()
    try:
        pixmap = QPixmap.grabWidget(canvas)
    except: # PyQt5+
        pixmap = canvas.grab()
    QApplication.clipboard().setPixmap(pixmap)
    print("Image copied to clipboard.")


def update_and_copy(fig=None, tight=True):
    if fig is None:
        fig = plt.gcf()
    if tight:
        fig.tight_layout()
    copy(fig)


def process_events(fig=None, activate_window=True):
    if fig is not None:
        fig.canvas.draw()
        if activate_window:
            fig.canvas.manager.window.activateWindow()
            fig.canvas.manager.window.raise_()
    QApplication.processEvents()


def wait_and_update(t, sleeptime=0.05, fig=None):
    t0 = time.time()
    while time.time() < t0 + t:
        time.sleep(sleeptime)
        process_events(fig)


# Creating and formatting figures
def get_fig(widths, heights, margins=0.5, dw=0.2, dh=0.2, make_axes=True):
    """
    Create a figure and grid where all dimensions are specified in inches.
    Arguments:
        widths: list of column widths
        heights: list of row heights
        margins: either a scalar or a list of four numbers (l, r, t, b)
        dw: white space between subplots, horizontal
        dh: white space between subplots, vertical
        make_axes: bool; if True, create axes on the grid and return,
                   else return the gridspec.
    """
    wsum = sum(widths)
    hsum = sum(heights)
    nrows = len(heights)
    ncols = len(widths)
    if type(margins) == list:
        l, r, t, b = margins
    else:
        l = r = t = b = margins

    figw = wsum + (ncols-1) * dw + l + r
    figh = hsum + (nrows-1) * dh + t + b

    # margins in fraction of the figure
    top = 1.-t/figh
    bottom = b/figh
    left = l/figw
    right = 1.-r/figw

    # subplot spacing in fraction of the subplot size
    wspace = dw / np.average(widths)
    hspace = dh / np.average(heights)

    fig = plt.figure(figsize=(figw, figh))
    gs = gridspec.GridSpec(nrows, ncols,
                           height_ratios=heights, width_ratios=widths)
    gs.update(top=top, bottom=bottom, left=left, right=right,
              wspace=wspace, hspace=hspace)

    if make_axes:
        axes = []
        for i in range(nrows):
            for j in range(ncols):
                axes.append(fig.add_subplot(gs[i,j]))

        return fig, axes

    else:
        return fig, gs


# Colors
def get_color_cycle(n, colormap=default_cmap, start=0., stop=1., format='hex'):
    pts = np.linspace(start, stop, n)
    if format == 'hex':
        colors = [rgb2hex(colormap(pt)) for pt in pts]
    return colors


def get_color(val, colormap=default_cmap, vmin=0, vmax=1, format='hex'):
    delta = float(vmax)-float(vmin)
    cval = val/delta - vmin/delta
    if cval > 1:
        cval = 1
    if cval < 0:
        cval = 0

    return rgb2hex(colormap(cval))


# tools for prettier plotting
def pplot(ax, x, y, yerr=None, linex=None, liney=None, color='k', fmt='o',
          alpha=1, mew=0.5, **kw):

    zorder = kw.pop('zorder', 0)
    line_dashes = kw.pop('line_dashes', [])
    line_lw = kw.pop('line_lw', 2)
    line_alpha = kw.pop('line_alpha', 0.5)
    line_color = kw.pop('line_color', color)
    line_zorder = kw.pop('line_zorder', -1)
    line_from_ypts = kw.pop('line_from_ypts', False)
    elinewidth = kw.pop('elinewidth', 0.5)
    label = kw.pop('label', None)
    label_x = kw.pop('label_x', x[-1])
    label_y_ofs = kw.pop('label_y_ofs', 0)
    label_kw = kw.pop('label_kw', {})
    fill_color = kw.pop('fill_color', 'w')

    syms = []

    if linex is None:
        linex = x

    if type(liney) == str:
        if liney == 'data':
            liney = y

    if yerr is not None:
        err = ax.errorbar(x, y, yerr=yerr, fmt='none', ecolor=color, capsize=0,
                          elinewidth=elinewidth, zorder=zorder)
        syms.append(err)

    if liney is None and line_from_ypts:
        liney = y.copy()

    if liney is not None:
        line, = ax.plot(linex, liney, dashes=line_dashes, lw=line_lw,
                        color=line_color, zorder=line_zorder, alpha=line_alpha)
        syms.append(line)
    if fill_color == 'same':
        fill_color = color
    fill, = ax.plot(x, y, fmt, mec='none', mfc=fill_color, alpha=alpha,
                    zorder=zorder, **kw)
    edge, = ax.plot(x, y, fmt, mec=color, mfc='None', mew=mew,
                    zorder=zorder, **kw)
    syms.append(fill)
    syms.append(edge)

    if label is not None:
        label_idx = np.argmin(np.abs(x-label_x))
        ax.annotate(label, (label_x, y[label_idx] + label_y_ofs),
                    color=color, **label_kw)

    return tuple(syms)


def ppcolormesh(ax, x, y, z, cmap=default_cmap, make_grid=True, set_lims=True, **kw):
    if make_grid:
        _x, _y = pcolorgrid(x, y)
    else:
        _x, _y = x, y

    im = ax.pcolormesh(_x, _y, z, cmap=cmap, **kw)
    if set_lims:
        ax.set_xlim(_x.min(), _x.max())
        ax.set_ylim(_y.min(), _y.max())

    return im


def waterfall(ax, xs, ys, offset=None, style='pplot', **kw):
    cmap = kw.pop('cmap', default_cmap)
    linex = kw.pop('linex', xs)
    liney = kw.pop('liney', None)

    ntraces = ys.shape[0]
    if offset is None:
        offset = ys.max()-ys.min()

    colorseq = None
    if 'color' not in kw:
        colorseq = get_color_cycle(ntraces, colormap=cmap)

    for iy, yvals in enumerate(ys):
        x = xs if len(xs.shape) == 1 else xs[iy]
        y = yvals + iy * offset
        lx = linex if len(linex.shape) == 1 else linex[iy]
        ly = None if liney is None else liney[iy] + iy * offset
        if 'color' in kw:
            color = kw['color']
        else:
            color = colorseq[iy]

        if style == 'pplot':
            pplot(ax, x, y, linex=lx, liney=ly, color=color, **kw)
        elif style == 'lines':
            ax.plot(x, y, '-', color=color, **kw)


def phist1d(ax, centers, hist, horizontal=False, **kw):
    fc = kw.pop('fc', 'silver')
    lw = kw.pop('lw', 0)
    edgecolor = kw.pop('edgecolor', 'w')
    width_scale = kw.pop('width_scale', 1)

    w = (centers[1]-centers[0]) * width_scale
    if horizontal:
        ax.barh(centers, hist, align='center', height=w, fc=fc,
                lw=lw, edgecolor=edgecolor)
    else:
        ax.bar(centers, hist, align='center', width=w, fc=fc, lw=lw,
               edgecolor=edgecolor)


# tools for better color plotting
def centers2edges(arr):
    e = (arr[1:] + arr[:-1])/2.
    e = np.concatenate(([arr[0]-(e[0]-arr[0])], e))
    e = np.concatenate((e, [arr[-1]+(arr[-1]-e[-1])]))
    return e

def pcolorgrid(xaxis, yaxis):
    xedges = centers2edges(xaxis)
    yedges = centers2edges(yaxis)
    xx, yy = np.meshgrid(xedges, yedges)
    return xx, yy


def cmap_powerlaw_adjust(cmap, a):
    '''
    returns a new colormap based on the one given
    but adjusted via power-law:

    newcmap = oldcmap**a
    '''
    if a < 0.:
        return cmap
    from copy import copy
    cdict = copy(cmap._segmentdata)
    fn = lambda x : (x[0]**a, x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), \
            "Resulting indices extend out of the [0, 1] segment."
    return mpl.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_center_adjust(cmap, center_ratio):
    '''
    returns a new colormap based on the one given
    but adjusted so that the old center point higher
    (>0.5) or lower (<0.5)
    '''
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    return cmap_powerlaw_adjust(cmap, a)


def cmap_center_point_adjust(cmap, range, center):
    '''
    converts center to a ratio between 0 and 1 of the
    range given and calls cmap_center_adjust(). returns
    a new adjusted colormap accordingly
    '''
    if not ((range[0] < center) and (center < range[1])):
        return cmap
    return cmap_center_adjust(cmap,
        abs(center - range[0]) / abs(range[1] - range[0]))


### plotting of IQ distributions
def plot_IQdist(ax, xedges, hist, norm=None, center=None, div=False,
                cmap=None, epsilon=0.01, ticklabels=None, make_edges=False,
                cax=None,
                contourlines=False, contourlines_color='k',
                contours=False, ncontours=8,
                **kw):

    rasterized = kw.pop('rasterized', True)
    cax_fontsize = kw.pop('cax_fontsize', 'x-small')
    grid = kw.pop('grid', True)
    gridcolor = kw.pop('gridcolor', '0.5')

    if cmap is None:
        if div:
            cmap = cmap_div
        else:
            cmap = cmap_seq

    if xedges is None:
        xedges = np.linspace(-1, 1, hist.shape[0]+1)
    elif make_edges:
        xedges = centers2edges(xedges)

    if norm is None:
        vmin = hist.min()
        vmax = hist.max()
        if vmin == vmax:
            vmax += epsilon
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if center is not None:
        if norm.vmin >= center:
            norm.vmin = center - epsilon
        if norm.vmax <= center:
            norm.vmax = center + epsilon
        cmap = cmap_center_point_adjust(cmap, (norm.vmin, norm.vmax), center)

    if not contours:
        im = ax.pcolormesh(xedges, xedges, hist, cmap=cmap, norm=norm,
                           rasterized=rasterized, **kw)
    else:
        dx = xedges[1]-xedges[0]
        xpts = xedges[:-1] + dx/2.
        im = ax.contourf(xpts, xpts, hist, ncontours, norm=norm, cmap=cmap)

    if contourlines:
        dx = xedges[1]-xedges[0]
        xpts = xedges[:-1] + dx/2.
        cntr = ax.contour(xpts, xpts, hist, ncontours, linewidths=0.5,
                          alpha=0.25, norm=norm, colors=contourlines_color)

    if grid:
        ax.grid(True, color=gridcolor, lw=0.5, dashes=[0.5, 2], alpha=1)

    if hasattr(ax, 'cax'):
        cax = ax.cax
    if cax is not None:
        cb = cax.colorbar(im)
        cax.tick_params(axis='both', labelsize=cax_fontsize)
        cb.solids.set_rasterized(rasterized)

    if contours:
        xedges = xpts

    ax.set_xlim(xedges.min(), xedges.max())
    ax.set_ylim(xedges.min(), xedges.max())

    if ticklabels is None:
        # ax.set_xticks([xedges.min(), 0, xedges.max()])
        # ax.set_yticks([xedges.min(), 0, xedges.max()])
        # _max = int(np.floor(xedges.max()))
        # ax.set_xticks(range(-_max//2, _max+1, 2))
        # ax.set_yticks(range(-_max//2, (_max+1)//2, 2))
        # ax.set_xticklabels([])
        # ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(3))
    elif ticklabels == 'off':
        ax.set_xticklabels([]); ax.set_yticklabels([])
    else:
        ax.set_xticks(ticklabels); ax.set_yticks(ticklabels)

    return im


def plot_IQdist_projections(xpts, hist, figsize=(3, 3), im_aspect='equal',
                            width_ratios=[2, 1], height_ratios=[1, 2],
                            cb_width=0.03, vcut=None, hcut=None,
                            prj_kw={}, **kw):

    xedges = centers2edges(xpts)

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 2, width_ratios=width_ratios,
                           height_ratios=height_ratios)
    gs.update(left=0.2, bottom=0.2)

    ax_img = fig.add_subplot(gs[1, 0], aspect=im_aspect)
    ax_hcut = fig.add_subplot(gs[0, 0], sharex=ax_img)
    ax_vcut = fig.add_subplot(gs[1, 1], sharey=ax_img)

    # construct the colorbar
    cbx = ax_vcut.get_position().bounds[0]
    cby = ax_hcut.get_position().bounds[1]
    cbw = cb_width
    cbh = ax_hcut.get_position().bounds[3]
    ax_cb_rect = [cbx, cby, cbw, cbh]
    ax_cb = fig.add_axes(ax_cb_rect)

    ax_cb.tick_params(axis='x', bottom='off', top='off', labelbottom='off')
    ax_cb.tick_params(axis='y', left='off', labelleft='off', labelright='on')
    ax_hcut.tick_params(axis='x', labelbottom='off')
    ax_vcut.tick_params(axis='y', labelleft='off')

    im = plot_IQdist(ax_img, xedges, hist, **kw)
    cb = fig.colorbar(im, cax=ax_cb)

    _prj_kw = {
        'lw': 0.,
        'fc': 'silver',
    }
    for k in prj_kw:
        _prj_kw[k] = prj_kw[k]

    if vcut is None:
        vcut = np.trapz(hist, x=xpts, axis=1)
    if hcut is None:
        hcut = np.trapz(hist, x=xpts, axis=0)

    phist1d(ax_vcut, xpts, vcut, horizontal=True, **_prj_kw)
    phist1d(ax_hcut, xpts, hcut, **_prj_kw)
    ax_vcut.set_xlim(0, None)
    ax_hcut.set_ylim(0, None)

    return fig, (ax_img, ax_hcut, ax_vcut, ax_cb), (im, cb), (hcut, vcut)
