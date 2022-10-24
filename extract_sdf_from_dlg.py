#!/usr/bin/env python3

from rdkit import Chem

import re


File = 'out-aq4.dlg'
Pat_origin = re.compile(r'INPUT-LIGAND-PDBQT: (.+)(?:\n|\r\n?)')
Pat_models = re.compile(r'DOCKED: MODEL.*?DOCKED: ENDMDL',re.DOTALL)

Autodock_atom_types = {
    'H':['H','HD','HS'],
    'C':['C','A'],
    'N':['N','NA','NS'],
    'O':['O','OA','OS'],
    'F':['F'],
    'Mg':['Mg','MG'],
    'P':['P'],
    'S':['S','SA'],
    'Cl':['Cl','CL'],
    'Br':['Br','BR'],
    'Ca':['Ca','CA'],
    'Mn':['Mn','MN'],
    'Fe':['Fe','FE'],
    'Zn':['Zn','ZN'],
    'I':['I'],
}

Contents = None
with open(File,'rt') as f:
    Contents = f.read()

if not Contents:
    print('Fatal: File Contents: no inputs')
    exit()



# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
#'HETATM    1  C1  AQ4 A 999      29.552  -1.989  53.496  1.00 82.82     0.018 C '
def read_pdbqt(lines):
    prolines = []
    for l in lines:
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
            try:
                x = float(l[30:38])
                y = float(l[38:46])
                z = float(l[46:54])
            except ValueError:
                print(f'Error: Line: {l}')
                continue
            else:
                prolines.append((at,x,y,z))
    return prolines

full = []
origin = Pat_origin.findall(Contents)
if not origin:
    print('Fatal: no found original structure')
    exit()
full.append(origin)



models = Pat_models.findall(Contents)
if not models:
    print('Fatal: No models found')
    exit()

for m in models:
    t = read_pdbqt(m)
    if not t:
        print('Error: error on found model')
        continue
    full.append(t)


mol = Chem.Mol()
conf = Chem.Conformer(len(full[0]))







