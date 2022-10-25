#!/usr/bin/env python3

from rdkit import Chem

import re
import sys
import os

FEATURES = [
    'version 0.1.0  : Start',
    'version 0.2.0  : deal with error occured on `Chem.MolFromPDBBlock`',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


Usage = """
Usage: AutoDock dock-dlg file process:
    ./this.script.py  dock.dlg
"""
if len(sys.argv) != 2:
    print(Usage)
    exit()
File = sys.argv[1]
if os.path.splitext(File)[1].lower() != '.dlg':
    print(Usage)
    print('Fatal: only dlg file is supported')
    exit()


Pat_origin = re.compile(r'INPUT-LIGAND-PDBQT: (.+)(?:\n|\r\n?)')
Pat_models = re.compile(r'DOCKED: MODEL.*?DOCKED: ENDMDL',re.DOTALL)

Autodock_atom_types = {
    'H' : ['H','HD','HS'],
    'C' : ['C','A'],
    'N' : ['N','NA','NS'],
    'O' : ['O','OA','OS'],
    'F' : ['F'],
    'Mg': ['Mg','MG'],
    'P' : ['P'],
    'S' : ['S','SA'],
    'Cl': ['Cl','CL'],
    'Br': ['Br','BR'],
    'Ca': ['Ca','CA'],
    'Mn': ['Mn','MN'],
    'Fe': ['Fe','FE'],
    'Zn': ['Zn','ZN'],
    'I' : ['I'],
}

Contents = None
with open(File,'rt') as f:
    Contents = f.read()

if not Contents:
    print('Fatal: File Contents: no inputs')
    exit()


# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
#'HETATM    1  C1  AQ4 A 999      29.552  -1.989  53.496  1.00 82.82     0.018 C '
def func_pdbqt2pdb(lines):
    prolines = []
    ebinding = []
    for l in lines:
        if 'Free Energy of Binding' in l:
            g = l.split()
            if len(g) > 8:
                try:
                    v = float(g[7])
                except ValueError:
                    pass
                else:
                    ebinding.append(v)
        if len(l) >= 78 and l[:4] in ['HETA','ATOM']:
            at = l[77:79].strip()
            ra = None
            for k,v in Autodock_atom_types.items():
                if at in v:
                    ra = k
                    break
            if not ra:
                print(f'Error: Atomtype not defined: {l}')
                continue
            new = l[0:12] + ' {:2} '.format(ra) + l[16:54] + ' '*26
            prolines.append(new)
    if not ebinding:
        ebinding = [0.0, ]
    elif len(ebinding) > 1:
        print('Warning: multiple binding values')
    return prolines,min(ebinding)


full = []
titles = []
origin = Pat_origin.findall(Contents)
if not origin:
    print('Fatal: no found original structure')
    exit()
oripdb,e = func_pdbqt2pdb(origin)
if not oripdb:
    print('Error: no found original structure')
full.append(oripdb)
titles.append('Original')


models = Pat_models.findall(Contents)
if not models:
    print('Fatal: No models found')
    exit()

for m in models:
    m = m.replace('DOCKED: ','').split('\n')
    t,e = func_pdbqt2pdb(m)
    if not t:
        print('Error: error on found model')
        continue
    full.append(t)
    titles.append('BindingEnergy: {:}'.format(e))

base = 'out-' + os.path.splitext(File)[0] + '.sdf'
w = Chem.SDWriter(base)
for m,t in zip(full,titles):
    mol = Chem.MolFromPDBBlock('\n'.join(m))
    if not mol:
        print('\nError:')
        for i in m:
            print(i)
        print('\n')
        continue
    mol.SetProp('_Name',t)
    w.write(mol)
w.close()

print('DONE')



