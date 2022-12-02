import os
import numpy as np

from fpga_lib import dsl
from fpga_lib.constants import STREAM_INFO

BLOCK_FMT = 'Data_0x%04d[0x00]B%d_#%04d.bin'


def load_se_block(directory, block, card):
    se0 = load_block(directory, 'se0', block, card)
    se1 = load_block(directory, 'se1', block, card)
    return se0 + 1j*se1


def load_block(directory, stream, block=0, card=0):
    if stream is 'se':
        return load_se_block(directory, block, card)

    stream_num, dtype = STREAM_INFO[stream]
    block_name = BLOCK_FMT % (stream_num, card, block)
    path = os.path.join(directory, block_name)

    n_attempts = 5
    n = 0
    while True:
        try:
            arr = np.fromfile(path, dtype=dtype)
            break
        except IOError:
            n += 1
            if n >= n_attempts:
                logging.error("Cannot load block, abort.")
                raise
            else:
                logging.warn("Cannot load block, will retry ({}/{})".format(n, n_attempts))
    
    if stream == 'rel':
        I = ((arr[::2] & 0x1ffff) - (arr[::2] & 0x20000)).astype(np.float)
        Q = ((arr[1::2] & 0x1ffff) - (arr[1::2] & 0x20000)).astype(np.float)
        return (I + 1j*Q) / float(2**16)
    if stream in ('se0', 'se1'):
        return arr / float(2**15)
    if stream == 'rec':
        return arr.reshape((len(arr)/8, 8))
    return arr