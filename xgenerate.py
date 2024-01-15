#!/usr/bin/env python3

import numpy as np
import rmsd
from rdkit import Chem
from rdkit.Chem import AllChem

import os
import sys
import random
import threading
import queue
import argparse


FEATURES = [
    'version 0.3.0    : add option `--insmiles`',
    'version 0.4.0    : add option `centroid_xyz`',
    'version 0.5.0    : make `--insmiles` to accept multiple inputs',
    'version 0.6.0    : `--autofiles`',
    'version 0.7.0    : filtration on pattern mol',
    'version 0.8.0    : refactor',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION

SETTINGS = {
    'cutoff_energy' : 0.5,
    'cutoff_rmsd'   : 0.5,
    'seed'          : -1,
    'forcefiled'    : 'uff',
    'forcefield_choices': ('uff', 'mmff94', 'mmff94s'),
    'outdir'        : '_xconformers',
    'num_confs'     : 5,
    'total_cpu'     : os.cpu_count(),
    'write_hydrogens': False,
}
SETTINGS['num_procs'] = SETTINGS['total_cpu'] - 2
if SETTINGS['num_procs'] <= 0: SETTINGS['num_procs'] = -1


class ConformerGenerator:
    """Generate conformers using RDKit

    Args:
        num_confs(int)      : number of conformers to generate
        cutoff_rmsd(float)  :
        cutoff_energy(float):
        forcefield(str)     :
        seed(int)           : random seed
        xyz(List[float])    : whether centroid onto this given coordinates
    """
    filebase = 'abcdefghijklmnopqrstuvwxyz'
    def __init__(
        self, num_confs=None, cutoff_rmsd=None, cutoff_energy=None,
        forcefield=None, seed=None, xyz=None, outdir=None, write_hydrogens=None,
        *args,**kws
    ):
        self.num_confs = num_confs
        self.cutoff_rmsd = cutoff_rmsd if cutoff_rmsd else 0.5
        self.cutoff_energy = cutoff_energy if cutoff_energy else 0.5
        self.forcefield = forcefield        # uff, mmff94
        self.seed = seed
        self.xyz = xyz
        self.outdir = outdir if outdir else '.'
        if not os.path.isdir(self.outdir): os.mkdir(self.outdir)
        self.write_hydrogens = True if write_hydrogens is True else False

    def generate_conformers(self, mol, num_confs=None, seed=None):
        if not num_confs: num_confs = self.num_confs
        if not seed: seed = self.seed
        file = self._get_filename(mol)
        if not file: return     # check before everything
        mol = self.embed_molecule(mol, num_confs, seed)
        energies = self.minimize_conformers(mol)
        self.filter_conformers(mol, num_confs, energies, self.cutoff_energy, self.cutoff_rmsd)
        if not self.write_hydrogens:
            mol = Chem.RemoveHs(mol)
        self.centroid_onto_xyz(mol,xyz=self.xyz)
        self.write(mol,file)

    def embed_molecule(self, mol, num_confs, seed=None):
        if not isinstance(mol,Chem.Mol): return
        print('  -> embedding')
        if seed is None: seed = -1
        brot = AllChem.CalcNumRotatableBonds(mol)
        tot = num_confs * 5
        if brot < 8:
            n = max(50,tot)
        elif brot <= 12:
            n = max(100,tot)
        elif brot <= 18:
            n = max(200,tot)
        else:
            n = max(300,tot)
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
        AllChem.EmbedMultipleConfs(
            mol, numConfs=tot, maxAttempts=n, clearConfs=True,
            pruneRmsThresh=-1.0, randomSeed=seed, ignoreSmoothingFailures=True,
        )
        return mol

    def _get_forcefield(self, mol, confId):
        if self.forcefield == 'uff':
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
        else:       # MMFF94, MMFF94s
            AllChem.MMFFSanitizeMolecule(mol)
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=self.forcefield)
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=confId)
        return ff

    def minimize_conformers(self, mol):
        if not isinstance(mol,Chem.Mol): return
        print('  -> minimizing')
        energies = []
        for conf in mol.GetConformers():
            ff = self._get_forcefield(mol,confId=conf.GetId())
            ff.Minimize()
            energies.append(ff.CalcEnergy())
        return energies

    def filter_conformers(self, mol, num_confs, energies=None, cutoff_energy=None, cutoff_rmsd=None):
        """Filter conformers based on Energy and RMSD threshold"""
        if not isinstance(mol,Chem.Mol): return
        if mol.GetNumConformers() < num_confs: return
        print('  -> filtering')
        ndxes = sorted(range(len(energies)),key=lambda x: energies[x])
        back = ndxes[:]
        get = [ndxes.pop(0), ]          # real index
        mark = {}
        while True:
            # energy difference
            for i in ndxes:
                ediff = [abs(energies[v]-energies[i]) for v in get]
                if all([v>=cutoff_energy for v in ediff]):
                    # RMSD difference
                    rdiff = []
                    for j in get:
                        if (i,j) in mark:
                            v = mark[(i,j)]
                        else:
                            posi = mol.GetConformer(i).GetPositions()
                            posj = mol.GetConformer(j).GetPositions()
                            v = rmsd.rmsd(posi, posj)
                            mark[(i,j)] = mark[(j,i)] = v
                        rdiff.append(v)
                    if all([v>=cutoff_rmsd for v in rdiff]):
                        get.append(i)
                        if len(get) >= num_confs:
                            break

            if len(get) < num_confs:
                cutoff_energy -= 0.05
                cutoff_rmsd -= 0.02
                ndxes = [i for i in back if i not in get]
            else:
                break

        for i in range(mol.GetNumConformers()-1,-1,-1):     # reversely update
            if i in get:
                mol.GetConformer(i).SetDoubleProp('energy',energies[i])
            else:
                mol.RemoveConformer(i)
        return get

    def centroid_onto_xyz(self, mol, xyz=None):
        if not xyz: return
        print('  -> onto xyz')
        n = mol.GetNumAtoms()
        for conf in mol.GetConformers():
            poses = conf.GetPositions()
            cent = np.sum(poses,axis=0) / n
            poses += np.array(xyz) - cent
            for i in range(n):
                conf.SetAtomPosition(i,poses[i])

    def _get_filename(self, mol):
        if mol.HasProp('_Name'):
            name = mol.GetProp('_Name')
            name = '_'.join(name.split())
            p = os.path.join(self.outdir,name)
            if not p.endswith('.sdf'): p += '.sdf'
            if os.path.isfile(p):
                s = os.stat(p)
                if s.st_mtime == s.st_ctime:
                    print(f'Warning: skip: file already exist: {p}')
                    return
        else:
            while True:
                n = random.randint(1,len(self.filebase))
                c = random.choices(self.filebase,k=n)
                name = ''.join(c)
                p = os.path.join(self.outdir,name)
                if not p.endswith('.sdf'): p += '.sdf'
                if not os.path.isfile(p): break
        return p

    def write(self, mol, file):
        print(f'  => writing to file: {file}')
        w = Chem.SDWriter(file)
        for conf in mol.GetConformers():     # confId may not sequentially starts from 1
            w.write(mol,confId=conf.GetId())
        w.close()


class ReadQueue(queue.Queue):
    def __init__(self, files, insmiles=None, maxsize=None):
        if not maxsize: maxsize = 10
        super().__init__(maxsize=maxsize)
        self.files = files
        self.insmiles = insmiles if insmiles else []

    def _read(self,file):
        if not (file and os.path.isfile(file)): return
        if file.endswith('.pdb'):
            try:
                yield Chem.MolFromPDBFile(file, removeHs=True)
            except:
                pass
        elif file.endswith('.sdf'):
            sup = Chem.SDMolSupplier(file, removeHs=True)
            while True:
                try:
                    yield next(sup)
                except StopIteration:
                    break
        elif file.endswith('.mol2'):
            try:
                yield Chem.MolFromMol2File(file, removeHs=True)
            except:
                pass
        elif file.endswith('.smi') or file.endswith('.smiles'):
            with open(file,'rt') as f:
                for line in f:
                    l = line.split()
                    if not l: continue
                    try:
                        mol = Chem.MolFromSmiles(l, removeHs=True)
                        if not mol:
                            raise StopIteration
                    except:
                        pass
                    else:
                        if len(l) >= 2:
                            mol.SetProp('_Name',l[1])
                        yield mol

    def read(self):
        for s in self.insmiles:
            try:
                mol = Chem.MolFromSmiles(s)
            except:
                pass
            else:
                self.put(mol)
        for f in self.files:
            print(f'Note: read file: {f}')
            it = self._read(f)
            for m in it:
                self.put(m)


def xmain():
    parser = argparse.ArgumentParser(
        description='Generate conformers from mol2 or SMILES',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=f"""
xconformer {__version__}:

# generate conformers from `mol2` file
>>>xconformer -m file.mol2 [-n 10] [--seed 20230905]

# generate conformers from `smiles` file
# the first two columns will be used for SMILES and the name of molecules
# if the name is not provided, a random name will be generated
>>>xconformer -s file.smiles [-n 10] [--seed 20230905]

"""
    )
    parser.add_argument(
        '-f', '--files',
        metavar='f',
        nargs='+',
        help='Input files, format determined by extension, (pdb, mol2, sdf, smi/smiles)',
    )
    parser.add_argument(
        '-is','--insmiles',
        nargs='+',
        metavar='smiles',
        help='Input SMILES string, final file name will be randomly decided'
    )
    parser.add_argument(
        '-xyz','--xyz',
        type=float,
        nargs=3,
        metavar='v',
        help='Whether centroid results onto xyz coordinates'
    )
    parser.add_argument(
        '-n', '--num-conf',
        default=SETTINGS['num_confs'],
        type=int,
        help='Number of conformations to be generated'
    )
    parser.add_argument(
        '--seed',
        default=SETTINGS['seed'],
        type=int,
        help='Random seed for conformer generation',
    )
    parser.add_argument(
        '-r', '--cutoff-rmsd',
        default=SETTINGS['cutoff_rmsd'],
        type=float,
        help='RMSD cutoff between conformers',
    )
    parser.add_argument(
        '-e', '--cutoff-energy',
        default=SETTINGS['cutoff_energy'],
        type=float,
        help='Energy cutoff between conformers',
    )
    parser.add_argument(
        '-ff', '--forcefield',
        default=SETTINGS['forcefiled'],
        choices=SETTINGS['forcefield_choices'],
        help='Type of forcefield to be used',
    )
    parser.add_argument(
        '-o', '--outdir',
        default=SETTINGS['outdir'],
        help='Folder to save conformers',
    )
    parser.add_argument(
        '-p','--num-procs',
        default=SETTINGS['num_procs'],
        type=int,
        help=(
            f'Number of processors to use, total cpus: {SETTINGS["total_cpu"]}'
        )
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='Show development features',
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        return

    args = parser.parse_args()
    if args.features:
        for i in FEATURES: print(i)
        return

    kws = SETTINGS.copy()       # shallow
    kws.update(dict(args._get_kwargs()))
    return kws


if __name__ == '__main__':
    kws = xmain()
    if not kws: sys.exit()

    def worker():
        while True:
            mol = q.get()
            if mol is None:
                break
            work.generate_conformers(mol)
            q.task_done()

    q = ReadQueue(files=kws['files'], insmiles=kws['insmiles'], maxsize=kws['num_procs']*2)
    q.read()

    threads = []
    work = ConformerGenerator(**kws)
    for _ in range(kws['num_procs']):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for i in range(len(threads)): q.put(None)
    for t in threads: t.join()


