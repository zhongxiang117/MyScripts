#!/usr/bin/env python3


from rdkit import Chem
from rdkit.Chem import Draw

import os
import sys
import argparse


FEATURES = [
    'version 0.1.0    : file to images',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


class ReadQueue:
    def __init__(self, files, insmiles=None):
        self.files = files if files else []
        self.insmiles = insmiles if insmiles else []

    def read_from_file(self,file=None):
        if not (file and os.path.isfile(file)): return
        if file.endswith('.pdb'):
            try:
                mol = Chem.MolFromPDBFile(file, removeHs=False)
            except:
                pass
            else:
                if not (mol.HasProp('_Name') and len(mol.GetProp('_Name').strip())):
                    mol.SetProp('_Name',os.path.basename(file))
                yield mol
        elif file.endswith('.sdf'):
            sup = Chem.SDMolSupplier(file, removeHs=False)
            n = 0
            name = os.path.splitext(os.path.basename(file))[0]
            while True:
                n += 1
                try:
                    mol = next(sup)
                except StopIteration:
                    break
                else:
                    if not (mol.HasProp('_Name') and len(mol.GetProp('_Name').strip())):
                        mol.SetProp('_Name',name+'_n'+str(n))
                    yield mol
        elif file.endswith('.mol2'):
            name = os.path.splitext(os.path.basename(file))[0]
            with open(file, 'rt') as f:
                doc = f.readlines()
            ndxlist = [i for (i,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
            ndxlist.append(len(doc))
            for i,b in enumerate(ndxlist[:-1]):
                e = ndxlist[i+1]
                block = ''.join(doc[b:e])
                mol = Chem.MolFromMol2Block(block)
                if not mol: continue
                if not (mol.HasProp('_Name') and len(mol.GetProp('_Name').strip())):
                    mol.SetProp('_Name', name+'_'+str(i+1))
                yield mol
        elif file.endswith('.smi') or file.endswith('.smiles'):
            name = os.path.splitext(os.path.basename(file))[0]
            with open(file,'rt') as f:
                for i,line in enumerate(f):
                    l = line.replace(',',' ').split()       # consider `csv`
                    if not l: continue
                    try:
                        mol = Chem.MolFromSmiles(l[0])
                        if not mol:
                            raise ValueError
                    except:
                        pass
                    else:
                        if len(l) >= 2:
                            mol.SetProp('_Name',l[1])
                        else:
                            mol.SetProp('_Name',name+'_'+str(i+1))
                        yield mol

    def read_from_insmiles(self,insmiles=None):
        if not isinstance(insmiles,(tuple,list)): return
        for i,s in enumerate(insmiles,start=1):
            try:
                mol = Chem.MolFromSmiles(s)
            except:
                pass
            else:
                mol.SetProp('_Name',str(i))
                yield mol

    def __iter__(self):
        for f in self.files:
            print(f'Note: read from file: {f}')
            it = self.read_from_file(f)
            if not it: continue
            for m in it: yield m
        it = self.read_from_insmiles(self.insmiles)
        if it:
            for m in it: yield m

    __next__ = __iter__


def mols2img_gridfile(mols,file=None,row=None,col=None,label=True,title=True,size=None):
    mols = [Chem.Mol(m) for m in mols if m]      # make a copy
    if not mols: return

    for m in mols: m.RemoveAllConformers()
    if size is None: size = [500,500]
    if title:
        legends = [m.GetProp('_Name') for m in mols]
    else:
        legends = None
    if label:
        for m in mols:
            for a in m.GetAtoms():
                a.SetProp('atomNote',str(a.GetIdx()+1))

    tot = len(mols)
    if row:
        col = tot // row
        if col*row != tot: col += 1
    elif col:
        pass
        # row = tot // col
        # if col*row != tot: row += 1
    else:
        row = int(tot / pow(tot,0.5))
        col = tot // row
        if col*row != tot: col += 1

    im = Draw.MolsToGridImage(mols,molsPerRow=col,subImgSize=size,legends=legends)
    if file:
        if not (file.endswith('.jpg') or file.endswith('.png')):
            file = os.path.splitext(file)[0] + '.jpg'
    else:
        file = 'imgs.jpg'
    print('> writing to file: ',file)
    im.save(file)


def mols2img_files(mols,filebasename=None,label=True,title=True,size=None):
    mols = [Chem.Mol(m) for m in mols if m]      # make a copy
    if not mols: return

    for m in mols: m.RemoveAllConformers()
    if size is None: size = [500,500]
    if title:
        legends = [m.GetProp('_Name').strip() for m in mols]
    else:
        legends = ['' for m in mols]
    if label:
        for m in mols:
            for a in m.GetAtoms():
                a.SetProp('atomNote',str(a.GetIdx()+1))
    if not filebasename: filebasename = 'img'
    for i,t,m in zip(range(1,len(mols)+1),legends,mols):
        if t:
            file = filebasename + '_' + str(i) + '_' + t + '.jpg'
        else:
            file = filebasename + '_' + str(i) + '.jpg'
        print('> writing to file: ',file)
        Draw.MolToImageFile(m,file,size=size)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=f"""files2img {VERSION}""",
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--files',
        metavar='f',
        nargs='+',
        help='input files, [sdf, mol2, pdb, smiles/smi(first column will be SMILES)]'
    )
    parser.add_argument(
        '-is','--insmiles',
        metavar='smi',
        nargs='+',
        help='input SMILES'
    )
    parser.add_argument(
        '-T', '--not-show-title',
        action='store_false',
        dest='title',
        help='do not show title, (default:False)'
    )
    parser.add_argument(
        '-L', '--not-show-label',
        action='store_false',
        dest='label',
        help='do not show label, (default:False)'
    )
    parser.add_argument(
        '-s', '--size',
        metavar='i',
        type=int,
        nargs=2,
        help='image size, (default:[500,500])'
    )
    g = parser.add_argument_group('save to one combined grid file')
    g.add_argument(
        '-nr', '--row',
        metavar='v',
        type=int,
        help='rows for image (higher priority)'
    )
    g.add_argument(
        '-nc', '--col',
        metavar='v',
        type=int,
        help='images per column'
    )
    parser.add_argument(
        '-M',
        action='store_true',
        help='save to individual image file'
    )
    parser.add_argument(
        '-o', '--filename',
        help='(png,jpg), if save to individual file `-M`, it will be the basename'
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

    rq = ReadQueue(files=w.files, insmiles=w.insmiles)
    nrq = rq.__next__()     # fire up
    if w.M:
        mols = [m for m in nrq]
        mols2img_files(mols,filebasename=w.filename,label=w.label,title=w.title,size=w.size)
    else:
        if w.filename:
            base,ext = os.path.splitext(w.filename)
        else:
            base,ext = 'imgs', '.jpg'
        if w.row and w.col:
            r,c = w.row,w.col
        elif w.row:
            r,c = w.row,6
        elif w.col:
            r,c = 6,w.col
        else:
            r,c = 6,6
        n = 0
        first = True
        while True:
            mols = []
            bo = False
            for i in range(r*c):
                try:
                    m = nrq.__next__()
                except StopIteration:
                    bo = True
                    break
                else:
                    mols.append(m)
            if bo and first:
                name = base + ext
            else:
                n += 1
                name = base + '_' + str(n) + ext
            first = False
            if mols:
                if len(mols) < c: c = len(mols)
                mols2img_gridfile(mols,file=name,col=c,label=w.label,title=w.title,size=w.size)
            if bo: break


if __name__ == '__main__':
    main()



