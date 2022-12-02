import numpy as np

### tools
def demod(arr, dt=1e-9, nbins=20, frq=50e6, phase=0):
    tvals = np.arange(1,arr.shape[-1]+1)*dt
    isig = arr * np.cos(np.pi*2*tvals*frq + phase)
    qsig = arr * np.sin(np.pi*2*tvals*frq + phase)
    newshape = list(arr.shape)
    assert arr.shape[-1] % nbins == 0, \
        "number of samples needs to be an integer multiple of number of bins"
    newshape[-1] /= nbins
    newshape.append(nbins)
    newshape = tuple(newshape)
    return tvals.reshape((-1,nbins))[:,-1], \
        isig.reshape(newshape).sum(axis=-1)/float(nbins), \
        qsig.reshape(newshape).sum(axis=-1)/float(nbins)

def weighted_moving_average(x,y,step_size=None,width=1,width_points=None):
    ''' 
    step_size sets the bin width, and therefore how many points we'll get out
    width is gaussian width for the per-bin filter
    '''
    if step_size is None: # set step size to return same-sized array
        step_size = 0.5/((x[1]-x[0])*(y.shape[0]-1))
    if width_points is not None: # use width in points instead of in units
        width = (x[1]-x[0])*width_points
    bin_centers  = np.arange(np.min(x),np.max(x)-0.5*step_size,step_size)+0.5*step_size
    bin_avg = np.zeros(len(bin_centers))

    #We're going to weight with a Gaussian function
    def gaussian(x,amp=1,mean=0,sigma=1):
        return amp*np.exp(-(x-mean)**2/(2*sigma**2))

    for index in range(0,len(bin_centers)):
        bin_center = bin_centers[index]
        weights = gaussian(x,mean=bin_center,sigma=width)
        bin_avg[index] = np.average(y,weights=weights)

    return (bin_centers,bin_avg)