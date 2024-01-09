#!/usr/bin/env python3
"""My Chemoinformatics Utilities"""

import os
import math
import numpy as np


FEATURES = [
    'version 0.1.0  : Chemoinfo Utilities, Nov 2nd, 2022',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


PDB_STD_RESIDUES = {
    'ALA': ['N', 'C', 'O', 'CA', 'CB'],
    'CYS': ['N', 'C', 'O', 'CA', 'CB', 'SG'],
    'ASP': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'OD1', 'OD2'],
    'GLU': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
    'PHE': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'GLY': ['N', 'C', 'O', 'CA'],
    'HIS': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE': ['N', 'C', 'O', 'CA', 'CB', 'CG1', 'CG2', 'CD1'],
    'LYS': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ'],
    'LEU': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2'],
    'MET': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'SD', 'CE'],
    'ASN': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'OD1', 'ND2'],
    'PRO': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD'],
    'GLN': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
    'ARG': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'SER': ['N', 'C', 'O', 'CA', 'CB', 'OG'],
    'THR': ['N', 'C', 'O', 'CA', 'CB', 'OG1', 'CG2'],
    'VAL': ['N', 'C', 'O', 'CA', 'CB', 'CG1', 'CG2'],
    'TRP': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'C', 'O', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
}


ATOMIC_MASS = {
    'H' : 1.00794,   'He': 4.00260,   'Li': 6.94100,   'Be': 9.01218,   'B' : 10.81100,
    'C' : 12.01070,  'N' : 14.00670,  'O' : 15.99940,  'F' : 18.99840,  'Ne': 20.17970,
    'Na': 22.98977,  'Mg': 24.30500,  'Al': 26.98154,  'Si': 28.08550,  'P' : 30.97376,
    'S' : 32.06500,  'Cl': 35.45300,  'Ar': 39.94800,  'K' : 39.09830,  'Ca': 40.07800,
    'Sc': 44.95591,  'Ti': 47.86700,  'V' : 50.94150,  'Cr': 51.99610,  'Mn': 54.93805,
    'Fe': 55.84500,  'Co': 58.93320,  'Ni': 58.69340,  'Cu': 63.54600,  'Zn': 65.40900,
    'Ga': 69.72300,  'Ge': 72.64000,  'As': 74.92160,  'Se': 78.96000,  'Br': 79.90400,
    'Kr': 83.79800,  'Rb': 85.46780,  'Sr': 87.62000,  'Y' : 88.90585,  'Zr': 91.22400,
    'Nb': 92.90638,  'Mo': 95.94000,  'Tc': 98.00000,  'Ru': 101.07000, 'Rh': 102.90550,
    'Pd': 106.42000, 'Ag': 107.86820, 'Cd': 112.41100, 'In': 114.81800, 'Sn': 118.71000,
    'Sb': 121.76000, 'Te': 127.60000, 'I' : 126.90447, 'Xe': 131.29300, 'Cs': 132.90545,
    'Ba': 137.32700, 'La': 138.90550, 'Ce': 140.11600, 'Pr': 140.90765, 'Nd': 144.24000,
    'Pm': 145.00000, 'Sm': 150.36000, 'Eu': 151.96400, 'Gd': 157.25000, 'Tb': 158.92534,
    'Dy': 162.50000, 'Ho': 164.93032, 'Er': 167.25900, 'Tm': 168.93421, 'Yb': 173.04000,
    'Lu': 174.96700, 'Hf': 178.49000, 'Ta': 180.94790, 'W' : 183.84000, 'Re': 186.20700,
    'Os': 190.23000, 'Ir': 192.21700, 'Pt': 195.07800, 'Au': 196.96655, 'Hg': 200.59000,
    'Tl': 204.38330, 'Pb': 207.20000, 'Bi': 208.98038, 'Po': 209.00000, 'At': 210.00000,
    'Rn': 222.00000, 'Fr': 223.00000, 'Ra': 226.00000, 'Ac': 227.00000, 'Th': 232.03810,
    'Pa': 231.03588, 'U' : 238.02891, 'Np': 237.00000, 'Pu': 244.00000, 'Am': 243.00000,
    'Cm': 247.00000, 'Bk': 247.00000, 'Cf': 251.00000, 'Es': 252.00000, 'Fm': 257.00000,
    'Md': 258.00000, 'No': 259.00000, 'Lr': 262.00000, 'Rf': 261.00000, 'Db': 262.00000,
    'Sg': 266.00000, 'Bh': 264.00000, 'Hs': 277.00000, 'Mt': 268.00000, 'Ds': 281.00000,
    'Rg': 272.00000, 'Cn': 285.00000, 'Uuq': 289.00000, 'Uuh' : 292.00000
 }
"""Atomic Mass, dict, key is sensitive"""


def read_pdb(file):
    """read PDB file
    
    Return:
        models (dict): models with their atoms inside, for each model, atom in a format:
            [seq,atomname,altLoc,resname,chainid,resnum,x,y,z,occ,tempfactor,element,charge, lineno],
            sepcially, resnum,x,y,z and line-number are numbers, else are strings
    """
    if isinstance(file,list):
        filelines = file
    elif isinstance(file,str):
        if os.path.isfile(file) and file.endswith('.pdb'):
            with open(file,'rt') as f: filelines = f.readlines()
        else:
            print(f'Fatal: not a valid pdb file: {file}')
            return []
    else:
        print('Fatal: not a valid file/list')
        return []
    models = {-1:[], }
    model = -1
    for idx,line in enumerate(filelines):
        n = len(line)
        if n > 10 and line[:5].lower() == 'model':
                model = int(line.split()[1])
                if model not in models: models[model] = []
        if n >= 54:
            if line[:6].lower() in ['atom  ', 'hetatm']:
                seq = line[6:11].strip()
                atomname = line[12:16].strip()
                altLoc = line[16]
                resname = line[17:20].strip()
                chainid = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = line[54:60].strip() if n >= 60 else ''
                tempfactor = line[60:66].strip() if n >= 66 else ''
                element = line[76:78].strip() if n >= 78 else ''
                charge = line[78:80].strip() if n >= 80 else ''
                models[model].append([
                    seq,atomname,altLoc,resname,chainid,resnum,x,y,z,occ,tempfactor,element,charge,idx
                ])
    return models


class ReadPDBFile:
    def __init__(self,pdbfile=None,*args,**kws):
        if isinstance(pdbfile,str):
            self.pdbfile = pdbfile
            self.pdbfilelines = open(pdbfile).readlines()
            self.pdbmodels = read_pdb(self.pdbfilelines)
        elif isinstance(pdbfile,list):
            self.pdbfile = 'list-inputs-but-not-a-file'
            self.pdbfilelines = pdbfile
            self.pdbmodels = read_pdb(self.pdbfilelines)
        else:
            print('Fatal: invalid input: pdbfile')
            self.pdbfile = ''
            self.pdbfilelines = ''
            self.pdbmodels = {}
        self.models = self.pdbmodels.keys()

    def remove_dup_atoms(self,model):
        """based on resnum, keep largest occupancy of atom with the same names"""
        duplist = self.get_dup_atoms_idx(model)
        if not duplist: return model
        return [a for i,a in enumerate(model) if i not in duplist]

    def get_dup_atoms_idx(self,model):
        """get duplicate atoms index"""
        resnumdict = {}
        for i,a in enumerate(model):
            if a[5] in resnumdict:
                resnumdict.append(i)
            else:
                resnumdict[a[5]] = [i, ]
        duplist = []
        for nl in resnumdict.values():
            if len(nl) < 2: continue
            occdict = {}
            for i in nl:
                atomname = model[i][1]
                if atomname in occdict:
                    occdict[atomname].append(i)
                else:
                    occdict[atomname] = [i, ]
            for ol in occdict.values():
                if len(ol) < 2: continue
                occlist = []
                for j in ol:
                    occ = model[j][9]
                    occlist.append(occ if occ else 0.0)
                keepndx = occlist.index(max(occlist))
                duplist.extend([v for t,v in enumerate(occlist) if t != keepndx])
        return sorted(duplist)


class Kit:
    def __init__(self,*args,**kws):
        pass

    def calc_rmsd(self,a,b):
        ga, ta = self.centroid(a)
        gb, tb = self.centroid(b)
        c = [[vb[i]*ga[n][i] for i in range(3)] for n,vb in enumerate(gb)]
        # singular value decomposition
        r1, s, r2 = np.linalg.svd(c)
        # compute sign, remove mirroring
        if np.linalg.det(c) < 0:
            r2[2,:] *= -1.0
        u = np.dot(r1, r2)
        return u, ta, tb

    def centroid(self,atoms,inplace=None):
        ax = sum([i[0] for i in atoms]) / len(atoms)
        ay = sum([i[1] for i in atoms]) / len(atoms)
        az = sum([i[2] for i in atoms]) / len(atoms)
        if inplace:
            for a in atoms:
                a[0] -= ax
                a[1] -= ay
                a[2] -= az
            return atoms, [ax,ay,az]
        return [[a[0]-ax,a[1]-ay,a[2]-az] for a in atoms], [ax,ay,az]

    def calc_atoms_distance(self,a,b):
        vl = [a[i]-b[i] for i in range(3)]
        dd = sum([i*i for i in vl])
        return pow(dd,0.5)

    def calc_atoms_angle(self,a,b,c,tol=None,degree=None):
        tol = tol if tol else pow(10,-10)
        ba = [a[i]-b[i] for i in range(3)]
        bc = [c[i]-b[i] for i in range(3)]
        xxba = sum([i*i for i in ba])
        xxbc = sum([i*i for i in bc])
        if xxba < tol or xxbc < tol:
            if degree:
                return 90.0
            return math.pi/2
        v = sum([ba[i]*bc[i] for i in range(3)])
        costheta = v / pow(xxba,0.5) / pow(xxbc,0.5)
        theta = math.acos(costheta)
        if degree:
            return theta/math.pi*180.
        return theta

    def calc_atoms_dihedral(self,a,b,c,d,tol=None,degree=None):
        """
        dihedral:   a-bc-d:  costheta = norm(ABxAC) * norm(DBxDC)

        sign is then defined as, if point d is in point a's left in righthanded
        cartesian system, no changes; otherwise, negate it.
        """
        tol = tol if tol else pow(10,-10)
        ab = [b[i]-a[i] for i in range(3)]
        ac = [c[i]-a[i] for i in range(3)]
        db = [b[i]-d[i] for i in range(3)]
        dc = [c[i]-d[i] for i in range(3)]
        nabc = [ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]]
        ndbc = [db[1]*dc[2]-db[2]*dc[1], db[2]*dc[0]-db[0]*dc[2], db[0]*dc[1]-db[1]*dc[0]]
        aa = sum([i*i for i in nabc])
        dd = sum([i*i for i in ndbc])
        if aa < tol or dd < tol:
            return 0.0
        costheta = sum([i*j for i,j in zip(nabc,ndbc)]) / pow(aa*dd,0.5)
        theta = math.acos(costheta)
        v = [nabc[1]*ndbc[2]-nabc[2]*ndbc[1], nabc[2]*ndbc[0]-nabc[0]*ndbc[2], nabc[0]*ndbc[1]-nabc[1]*ndbc[0]]
        cb = [c[i]-b[i] for i in range(3)]
        if sum([v[i]*cb[i] for i in range(3)]) < 0:
            theta = -theta
        if degree:
            return theta/math.pi*180.
        return theta







