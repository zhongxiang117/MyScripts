#!/usr/bin/env python3

import numpy as np

import os
import io
import tarfile
import hashlib
import pickle
from PIL import Image


# https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
def md5(fname):
    size = 8192
    hash_md5 = hashlib.md5()
    with open(fname, 'rb') as f:
        #for chunk in iter(lambda: f.read(size), b''):
        #    hash_md5.update(chunk)
        chunk = f.read(size)
        while chunk:
            hash_md5.update(chunk)
            chunk = f.read(size)
    return hash_md5.hexdigest()


def extract_files_from_tgz(file,name_startswith=None):
    """Extract all files starting with given string to bytes from tar.gz archive"""
    if not isinstance(file,str) or not (file.endswith('.gz') or file.endswith('.tgz') \
        or file.endswith('.tar.gz')):
        print(f'Warning: not supported file format: {file}')
        return []
    myfiles = []
    with tarfile.open(file, 'r:gz') as tar:
        for member in tar:
            if member.isfile():
                filename = os.path.basename(member.name)
                if name_startswith is None:
                    myfiles.append(tar.extractfile(member).read())
                elif filename.startswith(name_startswith):
                    myfiles.append(tar.extractfile(member).read())
    return myfiles


def unpickle(dataset):
    d = None
    if hasattr(dataset,'read'):
        fname = '<fileobj>'
        d = pickle.load(dataset,encoding='bytes')
    else:
        fname = dataset
        with open(dataset, 'rb') as f:
            d = pickle.load(f, encoding='bytes')
    return d, fname


def load_cifar_10_dataset(dataset):
    """Load CIFAR-10 Dataset
    
    Returns:
        data (np.ndarray): [[pic,channel,row,column], ...]
        labels (np.int32): correspondent labels [0,9]
    """
    d, fname = unpickle(dataset)
    if (not d) or (b'data' not in d) or (b'labels' not in d):
        print(f'Warning: no data extracted for input: {fname}')
        return [], []
    dataset = np.zeros((10000, 3*32*32), dtype=np.int32)
    labels = np.zeros((10000), dtype=np.int32)
    for j in range(10000):
        dataset[j] = d[b'data'][j]      # broad assignment
        labels[j] = d[b'labels'][j]
    data = dataset.reshape(10000, 3, 32, 32)
    return data, labels


def load_cifar_100_dataset(dataset):
    """Load CIFAR-100 Dataset
    
    Returns:
        data (np.ndarray): [[pic,channel,row,column], ...]
        fine_labels (np.int32): correspondent labels [0,99]
        coarse_labels (np.int32): correspondent labels [0,19]
    """
    d, fname = unpickle(dataset)
    if (not d) or (b'data' not in d):
        print(f'Warning: no data extracted for input: {fname}')
        return [], [], []
    if (b'fine_labels' not in d) or (b'coarse_labels' not in d):
        print(f'Warning: mising labels for input: {fname}')
        return [], [], []
    dataset = np.zeros((10000, 3*32*32), dtype=np.int32)
    fine_labels = np.zeros((10000), dtype=np.int32)
    coarse_labels = np.zeros((10000), dtype=np.int32)
    for j in range(10000):
        dataset[j] = d[b'data'][j]      # broad assignment
        fine_labels[j] = d[b'fine_labels'][j]
        coarse_labels[j] = d[b'coarse_labels'][j]
    data = dataset.reshape(10000, 3, 32, 32)
    return data, fine_labels, coarse_labels


def load_cifar_labels(dataset,label):
    """Load META labels from dataset"""
    d, fname = unpickle(dataset)
    if isinstance(label,bytes):
        pass
    elif isinstance(label,str):
        label = str.encode(label)
    else:
        print(f'Warning: unsupported label type: {type(label)} --> (str,bytes)')
        return []
    if (not d) or (label not in d):
        print(f'Warning: no data extracted for input: {fname}')
        return []
    return [i.decode('utf-8') for i in d[label]]


def read_cifar_10_labels_from_tgzfile(tgzfile):
    datasets = extract_files_from_tgz(tgzfile,'batches.meta')
    if not datasets or len(datasets) != 1:
        print(f'Fatal: no data extracted from file: {tgzfile}')
        return []
    bio = io.BytesIO()
    bio.write(datasets[0])
    bio.seek(0)
    return load_cifar_labels(bio,'label_names')


def read_cifar_100_labels_from_tgzfile(tgzfile):
    """
    Returns:
        fine_labels   : List[str]
        coarse_labels : List[str]
    """
    datasets = extract_files_from_tgz(tgzfile,'meta')
    if not datasets or len(datasets) != 1:
        print(f'Fatal: no data extracted from file: {tgzfile}')
        return [],[]
    bio = io.BytesIO()
    bio.write(datasets[0])
    bio.seek(0)
    d, fname = unpickle(bio)
    if (b'fine_label_names' not in d) or (b'coarse_label_names' not in d):
        print(f'Fatal: missing data labels for file: {tgzfile}')
        return [],[]
    fine = [i.decode('utf-8') for i in d[b'fine_label_names']]
    coarse = [i.decode('utf-8') for i in d[b'coarse_label_names']]
    return fine, coarse


def read_cifar_tgzfile(tgzfile,name_startswith=None,cifar100=False):
    datasets = extract_files_from_tgz(tgzfile,name_startswith)
    if not datasets:
        print(f'Fatal: no data extracted from file: {tgzfile}')
        return []
    mydata = []
    for d in datasets:
        bio = io.BytesIO()
        bio.write(d)
        bio.seek(0)     # important
        if cifar100:
            mydata.append(load_cifar_100_dataset(bio))
        else:
            mydata.append(load_cifar_10_dataset(bio))
    return mydata


def read_cifar_10_datasets_from_tgzfile(tgzfile,test=False):
    """
    Returns:
      2-elements:
        data: List[ numpy.ndarray( (picutures, channel, row, column) ) ]
    """
    if test:
        return read_cifar_tgzfile(tgzfile, 'test_batch')
    return read_cifar_tgzfile(tgzfile,'data_batch')


def read_cifar_100_datasets_from_tgzfile(tgzfile,test=False):
    """
    Returns:
      3-elements:
        data: List[ numpy.ndarray( (picutures, channel, row, column) ) ]
    """
    if test:
        return read_cifar_tgzfile(tgzfile, 'test', True)
    return read_cifar_tgzfile(tgzfile,'train', True)


def dataset2image(data,filename=None):
    """
    Args:
        data : np.ndarray( [channel, h, w] )
    """
    if not filename: filename = 'image.png'
    md = np.zeros((*data.shape[1:],3),dtype=np.uint8)   # dtype is important
    md[:,:,0] = data[0,:,:]
    md[:,:,1] = data[1,:,:]
    md[:,:,2] = data[2,:,:]
    img = Image.fromarray(md)
    img.show()


# format: [(filename,md5), ...]
CIFAR_10 = ('cifar-10-python.tar.gz', 'c58f30108f718f92721af3b95e74349a')
CIFAR_100 = ('cifar-100-python.tar.gz', 'eb9058c3a382ffc7106e4002c42a8d85')



