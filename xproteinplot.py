#!/usr/bin/env python3

from rdkit import Chem
import numpy as np
from scipy.spatial import distance_matrix
import matplotlib
import matplotlib.pyplot as plt

import os
import sys
import argparse

FEATURES = [
    'version 0.1.0  : plot protein interactions',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION

DEFAULT_COLORMAPS = ['gnuplot2', 'brg', 'rainbow', 'jet', 'turbo']
DEFAULT_CHOOSES = ['CA', 'CB', 'N', 'C', 'O', '1', '2']


def chordimage(distmat, filename=None, colormap=None, labels=None):
    """
    Reference:
        http://deep.cs.umsl.edu/disteval
    Args:
        distmat List[List[float,float]]:    distance matrix
        colormap (str): https://matplotlib.org/stable/users/explain/colors/colormaps.html
        labels List[str|int]: the labels used for text
    """
    if colormap is None or colormap not in DEFAULT_COLORMAPS: colormap = 'jet'
    L = len(distmat)
    dp = 2. * np.pi / L
    if L < 10:
        stride = 1
    elif L <= 20:
        stride = 3
    else:
        stride = 5
    cdistmat = np.copy(distmat)
    cdistmat[np.isnan(cdistmat)] = 1000.0
    cdistmat = 4.0 / (cdistmat+0.000001)        # avoid being divided by zero
    cdistmat[cdistmat>4] = 4.0
    pair_intensity = {}
    for i in range(L):
        for j in range(i+1, L):
            if cdistmat[i,j] < 0.1: continue    # only show strong interaction
            if abs(i-j) <= stride: continue     # only show strides interaction
            pair_intensity[(i,j)] = cdistmat[i,j] * cdistmat[j,i]
    if len(pair_intensity) < 2:
        print('Fatal: too less interactions')
        return

    fig, ax = plt.subplots(1, 1, figsize=(16, 16))
    lxy = [(np.cos((i+0.5)*dp),np.sin((i+0.5)*dp)) for i in range(L)]
    n = 5 * L       # only keep those number of interactions
    cmap = matplotlib.colormaps[colormap]
    v = max(pair_intensity.values())    # used to normalize
    for p, intensity in sorted(pair_intensity.items(), key=lambda x: x[1], reverse=True)[:n]:
        gx, gy = lxy[p[0]], lxy[p[1]]
        c = (p[0]+p[1]) / 2 / L     # color in between
        ax.plot(
            [gx[0],gy[0]], [gx[1],gy[1]], color=cmap(c),
            linewidth=5*intensity/v, alpha=intensity/v, zorder=-1
        )

    stride = L // 36 + 1    # maximum 36 points/arcs
    dt = 360. / L
    if not (isinstance(labels,(list,tuple)) and len(labels) == L):
        labels = [i for i in range(1,L+1)]
    for i in range(L):
        mycolor = cmap(i/L)
        angle = i * dt      # make the center of arc in where text positioned
        arc = matplotlib.patches.Arc(
            xy=(0,0), width=2, height=2, linewidth=12, angle=angle, theta2=dt, color=mycolor
        )
        ax.add_patch(arc)
        if i % stride != 0: continue
        a = (i+0.5) * dt / 180. * np.pi
        x = (0.52*np.cos(a) + 0.04*np.cos(i * dp)) * 2 - 0.04
        y = (0.52*np.sin(a) + 0.02*np.sin(i * dp)) * 2 - 0.02
        ax.text(x, y, labels[i], fontsize=32, color=cmap(i/L))
    ax.axis('off')
    if filename is None: filename = 'chord.png'
    print(f'Note: saving figure to: {filename}')
    fig.savefig(filename)


class ProteinInteractionPlot:
    def __init__(self, file, choose=None, usemass=None, colormap=None, outfile=None, *args,**kws):
        """
        Args:
            file (str|Chem.Mol) : pdb type object
            choose (str): {'CA', 'CB', 'N', 'C', 'O', '1', '2'}
                -> 1: residue geometric center
                -> 2: residue backbone geometric center
                => if `usemass=True`, then mass-weighted center will be computed
                => if 'CB' does not exist, e.g. GLY, 'CA' will be used instead
        """
        default_chooses = DEFAULT_CHOOSES
        if choose is None:
            choose = ['CA']
        elif isinstance(choose,str):
            choose = choose.replace('+',' ').replace(',',' ').replace(';',' ')
            choose = [i.upper() for i in choose.split()]
        elif isinstance(choose,(list,tuple)):
            choose = list(choose)
        else:
            choose = ['CA']
        choose = [i for i in choose if i in default_chooses]
        usemass = True if usemass is True else False
        self.colormap = colormap
        self.outfile = outfile
        self.parse_file(file,choose,usemass)

    def parse_file(self,file,choose,usemass=None,colormap=None,outfile=None):
        if not file: return
        if not colormap: colormap = self.colormap
        if not outfile: outfile = self.outfile
        mol = None
        if isinstance(file,str):
            if os.path.isfile(file) and file.endswith('.pdb'):
                mol = Chem.MolFromPDBFile(file)
        elif isinstance(file,Chem.Mol):
            mol = file
        if not mol: return
        resdict = {}
        for a in mol.GetAtoms():
            r = a.GetPDBResidueInfo()
            if not r: continue
            l = resdict.setdefault(r.GetResidueNumber(),[])
            l.append(a)

        positions = []
        labels = []
        conf = mol.GetConformer()
        backbone = ['CA','N','C','O']
        for k in sorted(resdict.keys()):
            name, pos, mass = [], [], []
            for a in resdict[k]:
                name.append( a.GetPDBResidueInfo().GetName().strip().upper() )
                pos.append( conf.GetAtomPosition(a.GetIdx()) )
                mass.append( a.GetMass() )
            if '1' in choose:
                if usemass:
                    x = sum([p[0]*m for p,m in zip(pos,mass)]) / len(pos)
                    y = sum([p[1]*m for p,m in zip(pos,mass)]) / len(pos)
                    z = sum([p[2]*m for p,m in zip(pos,mass)]) / len(pos)
                else:
                    x = sum([p[0] for p in pos]) / len(pos)
                    y = sum([p[1] for p in pos]) / len(pos)
                    z = sum([p[2] for p in pos]) / len(pos)
            elif '2' in choose:
                have = set(name).intersection(backbone)
                idx = [name.index(i) for i in have]
                if usemass:
                    x = sum([pos[i][0]*mass[i] for i in idx]) / len(idx)
                    y = sum([pos[i][1]*mass[i] for i in idx]) / len(idx)
                    z = sum([pos[i][2]*mass[i] for i in idx]) / len(idx)
                else:
                    x = sum([pos[i][0] for i in idx]) / len(idx)
                    y = sum([pos[i][1] for i in idx]) / len(idx)
                    z = sum([pos[i][2] for i in idx]) / len(idx)
            else:
                have = set(name).intersection(choose)
                if have:
                    idx = [name.index(i) for i in have]
                else:
                    # consideration when choose atom does not exist
                    if 'CA' in name:
                        idx = [name.index('CA')]
                    elif 'C' in name:
                        idx = [name.index('C')]
                    elif 'N' in name:
                        idx = [name.index('N')]
                    elif 'O' in name:
                        idx = [name.index('O')]
                    else:
                        idx = [0]
                if usemass:
                    x = sum([pos[i][0]*mass[i] for i in idx]) / len(idx)
                    y = sum([pos[i][1]*mass[i] for i in idx]) / len(idx)
                    z = sum([pos[i][2]*mass[i] for i in idx]) / len(idx)
                else:
                    x = sum([pos[i][0] for i in idx]) / len(idx)
                    y = sum([pos[i][1] for i in idx]) / len(idx)
                    z = sum([pos[i][2] for i in idx]) / len(idx)
            positions.append([x,y,z])
            labels.append(k)        # `k` has been sorted

        print('Note: valid atoms choose:', ' '.join(choose))
        if usemass:
            print('Note: plot is on mass-weighted geometric average')
        else:
            print('Note: plot is on geometric average')
        distmat = distance_matrix(positions,positions)
        chordimage(distmat, filename=outfile, colormap=colormap, labels=labels)


def main():
    parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         usage=f"""
PROGRAM {VERSION}. plot protein interactions
""",
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--file',
        dest='file',
        help='input PDB file'
    )
    parser.add_argument(
        '-c','--choose',
        metavar='v',
        nargs='+',
        default=['CA'],
        help=(
            'atoms to be analyzed, multiple, currently supports are: ' +
            ' '.join(DEFAULT_CHOOSES) + ', special key `1` is for whole residues, ' +
            '`2` is for backbone, default is `CA`'
        )
    )
    parser.add_argument(
        '-m', '--usemass',
        action='store_true',
        help='calculate mass-weighted center instead'
    )
    parser.add_argument(
        '-cmap','--colormap',
        choices=DEFAULT_COLORMAPS,
        help='colormap to be plotted'
    )
    parser.add_argument(
        '-o','--outfile',
        help='file to be output, default: chord.png'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return
    w = parser.parse_args(sys.argv[1:])
    if w.features:
        for i in FEATURES: print(i)
        return
    w.choose = ' '.join(w.choose)
    ProteinInteractionPlot(**vars(w))


if __name__ == '__main__':
    main()



