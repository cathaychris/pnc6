import numpy as np
import h5py
import time
import os


### constants
TIMESTAMP_FORMAT = '%Y%m%d%H%M%S'
TIMESTR_FORMAT = '%Y-%m-%d %H:%M:%S'
INS_SETTINGS_DIR = r'C:\Data\settings\instruments'


### convert lmfit parameters to table
def lmfit_params_to_data(params):
    dtype = []
    for p in params:
        dtype.append((p, 'f8'))
    dset = np.zeros(2, dtype=dtype)
    dset[0] = tuple([params[p].value for p in params])
    dset[1] = tuple([params[p].stderr for p in params])
    return dset


### tools for nested dictionaries
def get_from_dict(data_dict, map_list):
    return reduce(lambda d, k: d[k], map_list, data_dict)


### tools for getting data from hdf5 containers
def get_attrs(link):
    attrs = {}
    for a in link.attrs:
        attrs[a] = link.attrs[a]
    return attrs


def get_content(link):
    contents = {}
    for k in link.keys():
        contents[k] = {}
        if type(link[k]) == h5py._hl.dataset.Dataset:
            contents[k]['value'] = link[k].value
        if type(link[k]) == h5py._hl.group.Group:
            contents[k] = get_content(link[k])
        
        contents[k]['attrs'] = get_attrs(link[k])
    
    contents['attrs'] = get_attrs(link)
    return contents


def get_content_from_file(fn, grp=None):
    with h5py.File(fn, 'r') as f:
        link = f[grp] if grp is not None else f
        return get_content(link)


def get_latest_from_file(fn, grp=None, level=1, exclude=['attrs']):
    content = get_content_from_file(fn, grp=grp)
    path = []

    for l in range(level):
        keys = np.array([k for k in content.keys() if k not in exclude])
        srt = np.argsort(keys)
        keys = keys[srt]
        path.append(keys[-1])
        content = content[keys[-1]]

    return path, content


### tools for getting hdf5 objects. use carefully -- files will be open
# def get_grp(fn, grpn, mode='r'):
#     f = h5py.File(fn, mode)
#     return f, f[grpn]


### tools for writing in hdf5 containers
def init_file(fn):
    if not os.path.exists(fn):
        with h5py.File(fn, 'w') as f:
            pass


def init_group(fn, grpn, **attrs):
    init_file(fn)
    with h5py.File(fn, 'a') as f:
        if grpn not in f:
            grp = f.create_group(grpn)
            grp.attrs['time_created'] = time.time()
            grp.attrs['time_created_str'] = time.strftime(TIMESTR_FORMAT)
        else:
            grp = f[grpn]
        for a in attrs:
            grp.attrs[a] = attrs[a]


def save_dataset(fn, grpn, name, data, attrs={}, *args, **kws):
    """
    Overwrites any existing array, sets data to <arr>.
    args and kws are passed on to create_dataset.
    if group does not exist, create.
    """
    init_group(fn, grpn)
    with h5py.File(fn, 'a') as f:
        grp = f[grpn]
        if name in grp:
            del grp[name]
        dset = grp.create_dataset(name, data=data, *args, **kws)
        dset.attrs['time_created'] = time.time()
        dset.attrs['time_created_str'] = time.strftime(TIMESTR_FORMAT)
        for a in attrs:
            dset.attrs[a] = attrs[a]


def save_msmt_data(fn, grpn=None, datasets={}, attrs={}, dataset_attrs={}, 
    msmt_title='', datetimesubgrp=False):
    
    if grpn is None or datetimesubgrp == True:
        daten = time.strftime("%Y%m%d")
        timen = time.strftime("%H%M%S")
        _grpn = "{}/{}".format(daten, timen)
        
        if grpn is None:
            grpn = _grpn
        else:
            grpn += '/'+_grpn

    if msmt_title != '':
        grpn += '_{}'.format(msmt_title)
    
    for d in datasets:
        if d not in dataset_attrs:
            dattrs = {}
        else:
            dattrs = dataset_attrs[d]
        save_dataset(fn, grpn, d, datasets[d], attrs=dattrs)

    init_group(fn, grpn, **attrs)
    return grpn



