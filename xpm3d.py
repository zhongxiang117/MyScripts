#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.interpolate import griddata

import os


def plot_pm3d(
        x,y,z,bins=None,xlabel=None,ylabel=None,title=None,fout=None,
        xmin=None,xmax=None,ymin=None,ymax=None,colormap=None,interpolate=None
    ):
    """
    works like the contour map `gnuplot pm3d`, http://www.gnuplot.info/docs/loc6259.html

    Args:
        x | y | z (np.array)) : 1D correspondent array, has the relation `z = f(x,y)`
    """
    assert len(x) == len(y)
    assert len(y) == len(z)

    if xmin is None: xmin = x.min()
    if xmax is None: xmax = x.max()
    if ymin is None: ymin = y.min()
    if ymax is None: ymax = y.max()
    if xlabel is None: xlabel = 'x'
    if ylabel is None: ylabel = 'y'
    if fout is None or not isinstance(fout,str): fout = 'plot_im3d.png'
    if colormap is None or colormap not in ['jet','rainbow','gist_rainbow', 'magma','brg']:
        colormap = 'rainbow'
    if interpolate is None or interpolate not in ['nearest','linear','cubic']:
        interpolate = 'linear'
    if bins is None: bins = int(len(x)/2*0.8)

    xx = np.linspace(xmin, xmax, bins)
    yy = np.linspace(ymin, ymax, bins)
    grid = np.array(np.meshgrid(xx, yy.T))
    grid = grid.reshape(2, grid.shape[1]*grid.shape[2]).T
    points = np.array([x, y]).T     # transpose for `griddata`
    zz = griddata(points, z, grid, method=interpolate)
    zz = zz.reshape(xx.shape[0], yy.shape[0])
    fig, ax = plt.subplots(1, 1, figsize=(30, 30))
    cax = ax.imshow(zz, extent=[xmin,xmax,ymin,ymax], origin='lower', cmap=colormaps[colormap])
    if title is not None: ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.colorbar(cax)
    plt.show()
    #fig.savefig(fout)


file = 'fes_2D.dat'
# format:   #  phi  psi  fes  *
if os.path.isfile(file):
    print(f'Fatal: not a file: {file}')
    exit()

data = np.loadtxt(file,comments='#',dtype=np.float_)
x = data[:,0]
y = data[:,1]
z = data[:,2]

plot_pm3d(x,y,z,xmin=-3.14,xmax=3.14,ymin=-3.14,ymax=3.14)




