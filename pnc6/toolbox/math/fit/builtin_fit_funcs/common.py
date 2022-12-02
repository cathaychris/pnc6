import numpy as np

def find_index_of(ar, v):
    '''
    Find the index in <ar> of the value closest to <v>.
    '''
    return np.argmin(np.abs(ar - v))

def determine_offset(ar):
    '''
    Determine offset of <ar>, by histogramming and returning the value
    that occurs most often.
    Mostly useful for peak fitting (i.e. Gaussian or Lorentzian)
    '''

    nbins = max(10, round(len(ar)/100.0))
    hist, edges = np.histogram(ar, bins=nbins)
    idx = np.argmax(hist)
    return (edges[idx] + edges[idx+1]) / 2.0

def determine_peak_area(xs, ys):
    '''
    Determine area under curve
    '''
    spacing = np.abs(np.average(xs[1:] - xs[:-1]))
    wid = np.sum(ys) * spacing
    if np.sum(ys)<0.001*np.max(ys):
        wid = 0.5 * np.max(ys) * spacing
    return wid

def determine_peak_width(xs, ys, idx, factor=0.5):
    '''
    Determine width of a peak at <idx>, offset should already be subtracted.
    '''

    peakval = ys[idx]
    if idx == 0:
        idx +=1
    if idx == len(xs):
        idx -=1
    left = find_index_of(ys[:idx], peakval*factor)#hi
    right = find_index_of(ys[idx:], peakval*factor) + idx
    return np.abs(xs[right] - xs[left])

