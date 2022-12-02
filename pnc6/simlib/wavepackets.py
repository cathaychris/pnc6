import numpy as np





def gauss(tlist, chop=6.):
    T = tlist[-1]-tlist[0]
    sigma = T/chop
    t0 = T/2 + tlist[0]
    env = np.exp( - (tlist - T/2)**2/(2 * sigma**2) )
    return env



ENV_MAP = {
            'gaussian' : gauss
            }



def get_wavepacket_shape(name, tlist, detuning, **kwargs):
    try:
        env_func = ENV_MAP[name]
    except KeyError:
        raise NameError('{} not implemented!'.format(name))

    env = env_func(tlist, **kwargs)
    rotation = np.exp(-1j*detuning*tlist)

    return env*rotation 