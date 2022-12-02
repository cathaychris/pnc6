import numpy as np


def centers2edges(arr):
    e = (arr[1:] + arr[:-1])/2.
    e = np.concatenate(([arr[0]-(e[0]-arr[0])], e))
    e = np.concatenate((e, [arr[-1]+(arr[-1]-e[-1])]))
    return e
