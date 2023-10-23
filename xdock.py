#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from openbabel import pybel
from rdkit.Chem import AllChem

from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4LigandPreparation, AD4ReceptorPreparation

import os
import io
import sys
import shutil
import argparse


FEATURES = [
    'version 0.1.0  : XAutoDock',
    'version 0.2.0  : add argparse',
    'version 0.3.0  : add option for modification of maps',
    'version 0.4.0  : add generation of ligand shape file',
    'version 0.5.0  : fix ligand shape problems',
    'version 0.5.1  : add more output info',
    'version 0.6.0  : add option `--centroid-at`',
    'version 0.7.0  : add option `--maps-dir`',
    'version 0.8.0  : add option `--gen-ligand-gridshape`',
    'version 0.9.0  : adjust parameters in `dpf`',
    'version 0.10.0 : make option `--maps-dir` more versatile',
    'version 0.11.0 : add option `--scale`',
    'version 0.12.0 : interactive for `ROOT` selection',
    'version 0.13.0 : add selection mode for `ROOT`',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
__version__ = VERSION

# maximum only 14 atom types can be defined
#AUTODOCK_SUPPORTED_TYPES = 'H,HS,HD,  C,A,  N,NA,NS,  OA,OS,  S,SA,  P,F,Cl,Br,I'
AUTODOCK_SUPPORTED_TYPES = 'H,HS,HD,  C,A,  N,NA,NS,  OA,OS,  S,SA,  F,Cl'


TEMPLATE_GPF = """# XAutoDock GPF
npts {npts}
gridfld {gridfld}
spacing 0.375                        # grid spacing
receptor_types {rectypes}
ligand_types {ligtypes}
receptor {receptor}
gridcenter {gridcenter}
smooth 0.5                           # store minimum energy w/in rad(A)
{maps}
dielectric -0.1465                   #
"""

TEMPLATE_DPF = """# XAutoDock DPF
#autodock_parameter_version 4.2       # used by autodock to validate parameter set
outlev 1                             # diagnostic output level
intelec                              # calculate internal electrostatics
seed pid time                        # seeds for random generator
ligand_types {ligtypes}
fld {gridfld}
{maps}
move {move}
about {about}
tran0 0.0 0.0 0.0                    # initial coordinates/A or random
quaternion0 0.0 0.0 0.0 1.0          # initial orientation or random
dihe0 random                         # initial dihedrals (relative) or random
#torsdof 8                            # torsional degrees of freedom
rmstol 2.0                           # cluster_tolerance/A
rmsmode heavy_atoms_only
extnrg 1000.0                        # external grid energy
e0max 0.0 10000                      # max initial energy; max number of retries
ga_pop_size 150                      # number of individuals in population
ga_num_evals 250000                  # maximum number of energy evaluations
ga_num_generations 27000             # maximum number of generations
ga_elitism 10                        # number of top individuals to survive to next generation
ga_mutation_rate 0.02                # rate of gene mutation
ga_crossover_rate 0.8                # rate of crossover
ga_window_size 10                    #
ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution
ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution
set_ga                               # set the above parameters for GA or LGA
sw_max_its 300                       # iterations of Solis & Wets local search
sw_max_succ 4                        # consecutive successes before changing rho
sw_max_fail 4                        # consecutive failures before changing rho
sw_rho 1.0                           # size of local search space to sample
sw_lb_rho 0.01                       # lower bound on rho
ls_search_freq 0.06                  # probability of performing local search on individual
set_psw1                             # set the above pseudo-Solis & Wets parameters
unbound_model bound                  # state of unbound ligand
ga_run 5                             # do this many hybrid GA-LS runs
analysis                             # perform a ranked cluster analysis
"""


def xnewfile(file,head=None,includedir=None):
    """calculate new file name to avoid overwriting
    
    if includedir is turned on, operation will be performed on the `file` path;
    otherwise, operation will be performed on current dir
    """
    if not head and not os.path.isfile(file): return file
    dire = os.path.dirname(file)
    if not dire: dire = '.'
    base = os.path.basename(file)
    if not head: head = 'new-'
    i = 1
    while True:
        if includedir:
            new = dire + os.path.sep + head + str(i) + '-' + base
        else:
            new = head + str(i) + '-' + base
        if not os.path.isfile(new): break
        i += 1
    return new


class XADTPrep:
    def __init__(
        self,recfile=None,ligfile=None,centroid_at=None,
        npts=None,gridcenter=None,ligand_types=None,receptor_types=None,
        gen_lig_mode=None,
        verbose=None,*args,**kws
    ):
        self.verbose = True if verbose is True else False
        self.recfile = recfile if recfile else 'receptor.pdb'
        # be aware, basename should be used!
        self.recfile_s = os.path.basename(os.path.splitext(self.recfile)[0])
        self.recfile_s = xnewfile(self.recfile_s)

        self.ligfile = ligfile if ligfile else 'ligand.pdb'
        self.ligfile_s = os.path.basename(os.path.splitext(self.ligfile)[0])
        self.ligfile_s = xnewfile(self.ligfile_s)
        self.gridcenter = [0.0, 0.0, 0.0]
        self.centroid_at = centroid_at
        self.ligxyz = None
        self._flp = None
        self._gen_pdbqt_for_ligand()
        self.gen_lig_mode = gen_lig_mode

        self.npts_s = self._norm(npts,'60,60,60','npts')
        if not gridcenter: gridcenter = self.gridcenter     # calc'ed from above
        self.gc_s = self._norm(gridcenter,'auto','gridcenter')
        self.types_lig_s = self._norm(ligand_types,AUTODOCK_SUPPORTED_TYPES,'ligand_types')
        self.types_rec_s = self._norm(receptor_types,AUTODOCK_SUPPORTED_TYPES,'receptor_types')

    def _norm(self,l,default_s,info=None):
        bo = False
        if l is None:
            r = ' '.join(default_s.replace(',', ' ').replace(';', ' ').split())
        elif isinstance(l,str):
            r = ' '.join(l.replace(',', ' ').replace(';', ' ').split())
            if len(r) == 0:
                bo = True
        elif isinstance(l,(list,tuple)):
            r = ' '.join([str(round(i,3)) for i in l])
        else:
            r = ''
            bo = True
        if bo:
            if info:
                print(f'Warning: wrong defined: {info}: {l}')
            else:
                print(f'Warning: wrong defined: {l}')
            print('  >  defaults will be used')
            r = ' '.join(default_s.replace(',', ' ').replace(';', ' ').split())
        return r

    def gen_gpf(self):
        gpf = TEMPLATE_GPF.format(
            npts=self.npts_s,
            gridfld=self.recfile_s+'.maps.fld',
            rectypes=self.types_rec_s,
            ligtypes=self.types_lig_s,
            receptor=self.recfile_s+'.pdbqt',
            gridcenter=self.gc_s,
            maps='\n'.join(
                ['map '+self.recfile_s+'.'+t+'.map' for t in self.types_lig_s.split()] +
                ['elecmap '+self.recfile_s+'.e.map', 'dsolvmap '+self.recfile_s+'.d.map']
            )
        )
        self.write(gpf,self.recfile_s+'.gpf')

    def gen_dpf(self):
        dpf = TEMPLATE_DPF.format(
            ligtypes=self.types_lig_s,
            gridfld=self.recfile_s+'.maps.fld',
            maps='\n'.join(
                ['map '+self.recfile_s+'.'+t+'.map' for t in self.types_lig_s.split()] +
                ['elecmap '+self.recfile_s+'.e.map', 'dsolvmap '+self.recfile_s+'.d.map']
            ),
            move=self.ligfile_s+'.pdbqt',
            about=self.gc_s,
        )
        self.write(dpf,self.recfile_s+'.dpf')

    def write(self,contents,outfile=None):
        if not outfile: outfile = 'xautodock.file'
        print(f'Note: Writing to file: {outfile}')
        with open(outfile,'wt') as f:
            f.write(contents)

    def _gen_pdbqt_for_ligand(self,file=None,verbose=None):
        if not file: file = self.ligfile
        if not os.path.isfile(file): return
        if verbose is None: verbose = self.verbose
        mols = Read(file)
        if verbose:
            print(f'Note: Read ligand: {file}')
        mol = mols[0]
        coord_dict = {a:a.coords for a in mol.allAtoms}
        mol.buildBondsByDistance()
        self._flp = AD4LigandPreparation(
            mol,mode='',charges_to_add='gasteiger',root='auto',cleanup='nphs_lps',
            allowed_bonds='backbone',check_for_fragments=False,attach_singletons=False,
            attach_nonbonded_fragments=False,inactivate_all_torsions=False,
        )

        keys = list(coord_dict.keys())
        bad_list = [
            a for a in mol.allAtoms if a in keys and a.coords!=coord_dict[a]
        ]
        if len(bad_list):
            print(len(bad_list), ' atom coordinates changed!')
            for a in bad_list:
                print(a.name, ":", coord_dict[a], ' -> ', a.coords)
        elif verbose:
            print("No change in atomic coordinates")

        # if mol.returnCode != 0:
        #     print(f'ligand return msg: {mol.returnMsg}')
        self.ligxyz = mol.allAtoms.coords
        # calc gridcenter
        sx = sum([i[0] for i in self.ligxyz])
        sy = sum([i[1] for i in self.ligxyz])
        sz = sum([i[2] for i in self.ligxyz])
        n = len(mol.allAtoms)
        self.gridcenter = [round(sx/n,3), round(sy/n,3), round(sz/n,3)]

        if self.centroid_at:
            dx = self.gridcenter[0] - self.centroid_at[0]
            dy = self.gridcenter[1] - self.centroid_at[1]
            dz = self.gridcenter[2] - self.centroid_at[2]
            self.gridcenter = self.centroid_at
            self.ligxyz = [[v[0]-dx,v[1]-dy,v[2]-dz] for v in self.ligxyz]
            self._flp.molecule.allAtoms.coords = self.ligxyz

    def gen_pdbqt_for_ligand(self,file=None,verbose=None):
        if not file:
            file = self.ligfile
            if not self._flp:
                print(f'Warning: no valid ligand is detected: {self.ligfile}')
                return

        ip = InteractiveRoot(file)
        mols = Read(file)
        mol = mols[0]
        mol.buildBondsByDistance()
        kws = dict(
            mode='',charges_to_add='gasteiger',root='auto',cleanup='nphs_lps',
            allowed_bonds='backbone',check_for_fragments=False,attach_singletons=False,
            attach_nonbonded_fragments=False,inactivate_all_torsions=False,
        )

        if self.gen_lig_mode in [None,0]:       # interactive mode
            prev = None
            while True:
                if prev:
                    kws['root'] = prev - 1      # ROOT starts from zero
                flp = AD4LigandPreparation(mol,**kws)
                fp = io.StringIO()
                flp.writer.xzfpout = fp
                flp.write('placeholder.pdbqt')
                fp.seek(0)
                fout = fp.read()
                fp.close()
                rootcenter = ip.draw(fout)
                bo = False
                while True:
                    c = input('Input your choose (to accept, just press Enter): ')
                    if len(c.strip()) == 0:
                        bo = True
                        break
                    try:
                        c = int(c)
                        if c < 0 or c > len(mol.allAtoms):
                            raise IndexError
                    except:
                        print('Warning: wrong input')
                    else:
                        prev = c
                        break
                if bo:
                    self._flp = flp
                    flp.writer.xzfpout = None       # stop hacking
                    self.gc_s = self._norm(rootcenter,'auto','gridcenter')
                    break
        elif self.gen_lig_mode in [1,2]:            # number of atoms
            reflist = []
            atomndx = set(list(range(1,len(self.ligxyz)+1)))
            root = 'auto'
            while True:
                flp = AD4LigandPreparation(mol,**kws)
                fp = io.StringIO()
                flp.writer.xzfpout = fp
                flp.write('placeholder.pdbqt')
                fp.seek(0)
                fout = fp.read()
                fp.close()
                rootatomndx,rootcenter = ip.draw(fout,bool_draw=False)
                reflist.append((root,len(rootatomndx)))
                atomndx = atomndx.difference(rootatomndx)
                if len(atomndx):
                    root = atomndx.pop()
                    kws['root'] = root - 1          # ROOT starts from zero
                else:
                    break
            ndxes = sorted(reflist,key=lambda k:k[1])   # from smallest to bigest
            if self.gen_lig_mode == 2:
                print('Note: using the least number of atoms to set root')
                kws['root'] = ndxes[0][0] if ndxes[0][0] == 'auto' else ndxes[0][0]-1
            else:
                print('Note: using the most number of atoms to set root')
                kws['root'] = ndxes[-1][0] if ndxes[-1][0] == 'auto' else ndxes[-1][0]-1
            self._flp = AD4LigandPreparation(mol,**kws)
            self.gc_s = self._norm(rootcenter,'auto','gridcenter')

        print(f'Note: Writing ligand pdbqt to: {self.ligfile_s}.pdbqt')
        self._flp.write(self.ligfile_s+'.pdbqt')
        atomtypes = self.get_autodock_atomtypes(self.ligfile_s+'.pdbqt')
        atomtypes = list(set(atomtypes))
        self.types_lig_s = ' '.join(atomtypes)
        print('Note: processed ligand valid types: '+' '.join(atomtypes))
        print('Note: processed ligand overall atoms center: '+' '.join([str(i) for i in self.gridcenter]))

    def gen_pdbqt_for_protein(self,file=None,verbose=None):
        if not file: file = self.recfile
        if not os.path.isfile(file): return
        if verbose is None: verbose = self.verbose
        mols = Read(file)
        if verbose:
            print(f'Note: Read protein: {file}')
        mol = mols[0]
        mol.buildBondsByDistance()
        frp = AD4ReceptorPreparation(
            mol,mode='',charges_to_add='gasteiger',cleanup='nphs_lps_waters_nonstdres'
        )
        print(f'Note: Writing receptor pdbqt to: {self.recfile_s}.pdbqt')
        frp.write(self.recfile_s+'.pdbqt')
        atomtypes = self.get_autodock_atomtypes(self.recfile_s+'.pdbqt')
        atomtypes = list(set(atomtypes))
        print('Note: processed receptor types: '+' '.join(atomtypes))

    def gen_gridbox_pdb(self,outfile=None,gridcenter=None,npts=None,spacing=None):
        if not outfile: outfile = 'box.pdb'
        if not gridcenter: gridcenter = self.gridcenter
        t = self._norm(gridcenter,'0.0,0.0,0.0','gridcenter')
        r = list(map(float,t.split()))
        n = self._norm(npts,'60,60,60','npts')
        npts = list(map(int,n.split()))
        spacing = float(self._norm(spacing,'0.375','spacing'))
        #    5-------8
        #   /|      /|
        #  6-------7 |
        #   |/1----|/4
        #  2|-----|/3
        box = """HEADER    XAutoDock GridBox
            REMARK    CENTER (X Y Z)  {:}
            REMARK    NPTS   (X Y Z)  {:}
            REMARK    SPACING         {:}
            ATOM      1  DUA BOX     1    {:}
            ATOM      2  DUB BOX     1    {:}
            ATOM      3  DUC BOX     1    {:}
            ATOM      4  DUD BOX     1    {:}
            ATOM      5  DUE BOX     1    {:}
            ATOM      6  DUF BOX     1    {:}
            ATOM      7  DUG BOX     1    {:}
            ATOM      8  DUH BOX     1    {:}
            CONECT    1    2    4    5
            CONECT    2    1    3    6
            CONECT    3    2    4    7
            CONECT    4    1    3    8
            CONECT    5    1    6    8
            CONECT    6    2    5    7
            CONECT    7    3    6    8
            CONECT    8    4    5    7
            END
        """
        p = [i/2*spacing for i in npts]
        a = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]-p[0],r[1]-p[1],r[2]-p[2])
        b = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]+p[0],r[1]-p[1],r[2]-p[2])
        c = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]+p[0],r[1]+p[1],r[2]-p[2])
        d = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]-p[0],r[1]+p[1],r[2]-p[2])
        e = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]-p[0],r[1]-p[1],r[2]+p[2])
        f = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]+p[0],r[1]-p[1],r[2]+p[2])
        g = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]+p[0],r[1]+p[1],r[2]+p[2])
        h = '{:8.3f}{:8.3f}{:8.3f}'.format(r[0]-p[0],r[1]+p[1],r[2]+p[2])
        contents = box.format(t,n,spacing,a,b,c,d,e,f,g,h) # left-strip the blank spaces
        self.write('\n'.join([i.lstrip() for i in contents.split('\n')]),outfile=outfile)

    def get_autodock_atomtypes(self,file):
        if not os.path.isfile(file): return []
        raws = []
        with open(file,'rt') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    raws.append(line[77:79])
        if not raws:
            print(f'Warning: no valid ATOM/HETATM entry was found: {file}')
            return []
        pros = [i.strip() for i in raws if i.strip()]
        if len(pros) != len(raws):
            print(f'Warning: some ATOM/HETATM entry were missing: {file}')
            return []
        return pros


class ReadMap:
    """Read Map Files
    
    data format: [z[y[x[float]]]]
    """
    def __init__(
        self,mapfile=None,npts=None,*args,**kws
    ):
        self.mapfile = mapfile
        self.npts = npts if npts else [60, 60, 60]
        self._pars = []
        self.center = [0.0, 0.0, 0.0]
        self.space = 0.375
        self.r = 0.0
        self.d = 2.0

    def read_map_from_autodock(self,mapfile=None,npts=None):
        if not mapfile: mapfile = self.mapfile
        if not npts: npts = self.npts
        if not os.path.isfile(mapfile):
            print(f'Error: not a map file: {mapfile}')
            return []
        # try to find npts if not defined
        if not npts:
            with open(mapfile,'rt') as f:
                for line in f:
                    ltmp = line.split()
                    if len(ltmp) == 4 and ltmp[0].lower() == 'nelements':
                        try:
                            x = int(ltmp[1])
                            y = int(ltmp[2])
                            z = int(ltmp[3])
                        except ValueError:
                            pass
                        else:
                            npts = [x,y,z]
                            break
            if not npts:
                print('Error: npts: NELEMENTS not defined')
                return []
        self.npts = npts
        data = []
        with open(mapfile,'rt') as f:
            # dump first 6 lines
            self._pars = [f.readline().rstrip() for i in range(6)]
            try:
                data = [float(l.strip()) for l in f]
            except ValueError:
                data = []
        if not data:
            print('Error: not all values are the number, double check it')
            return []
        if len(data) != (npts[0]+1)*(npts[1]+1)*(npts[2]+1):
            print('Error: number of data inside map is not equal to defined points')
            return []
        print(f'Note: read total number of {len(data)} points')
        print(f'Note: minimum energy: {min(data)}')
        print(f'Note: maximum energy: {max(data)}')
        dx, dy, dz = npts[0]+1, npts[1]+1, npts[2]+1
        oz = dx * dy
        prodata = [
            [[data[x+y*dx+z*oz] for x in range(dx)] for y in range(dy)]
            for z in range(dz)
        ]

        self._get_center_grids()
        return prodata

    def save_autodock_map_from_data(self,data,filename=None,overwrite=None):
        if not filename: filename = self.mapfile
        if overwrite is not True:
            new = xnewfile(filename,'bak-',includedir=True)
            if filename != new:
                print(f'Note: Backup old map file to: {new}')
                shutil.move(filename,new)
        with open(filename,'wt') as f:
            for s in self._pars[:6]:
                f.write(f'{s}\n')
            for k in range(self.npts[2]+1):
                for j in range(self.npts[1]+1):
                    for i in range(self.npts[0]+1):
                        f.write(f'{data[k][j][i]}\n')
        print(f'Note: file saved to {filename}')

    def _get_center_grids(self):
        if not self._pars: return
        for l in self._pars[:6]:
            if l.startswith('SPACING'):
                try:
                    self.space = float(l.split()[1])
                except (ValueError,IndexError):
                    print(f'Warning: wrong line: {l}')
                    return
            elif l.startswith('CENTER'):
                ls = l.split()
                bo = False
                if len(ls) == 4:
                    try:
                        x = float(ls[1])
                        y = float(ls[2])
                        z = float(ls[3])
                    except ValueError:
                        bo = True
                    else:
                        self.center = [x,y,z]
                else:
                    bo = True
                if bo:
                    print(f'Warning: wrong line: {l}')
                    return
        self.xmin = self.center[0] - self.npts[0]/2*self.space
        self.ymin = self.center[1] - self.npts[1]/2*self.space
        self.zmin = self.center[2] - self.npts[2]/2*self.space

    def calc_indexs(self,ligxyz,r=None,d=None):
        if not r: r = self.r
        if not d: d = self.d
        rd = r + d
        overlaps = set()
        for c in ligxyz:
            ox = (int((c[0]-rd-self.xmin)/self.space),int((c[0]+rd-self.xmin)/self.space))
            oy = (int((c[1]-rd-self.ymin)/self.space),int((c[1]+rd-self.ymin)/self.space))
            oz = (int((c[2]-rd-self.zmin)/self.space),int((c[2]+rd-self.zmin)/self.space))
            for i in range(ox[0],min(ox[1]+1,self.npts[0]+1)):
                for j in range(oy[0],min(oy[1]+1,self.npts[1]+1)):
                    for k in range(oz[0],min(oz[1]+1,self.npts[2]+1)):
                        overlaps.add((i,j,k))
        alls = set()
        for i in range(self.npts[0]+1):
            for j in range(self.npts[1]+1):
                for k in range(self.npts[2]+1):
                    alls.add((i,j,k))
        left = alls.difference(overlaps)
        return overlaps,left

    def fill_values(self,data,indexes,fillval=None,scale=None):
        if scale:
            print(f'Note: scaleing value on grid points by: {scale}')
            for ndx in indexes:
                data[ndx[0]][ndx[1]][ndx[2]] *= scale
            return
        if not fillval: fillval = 0.0
        print(f'Note: unifying value on grid points by: {fillval}')
        for ndx in indexes:
            data[ndx[0]][ndx[1]][ndx[2]] = fillval


class LigandShape(XADTPrep):
    def __init__(
        self,mapfiles=None,maps_dir=None,gen_mapshape=None,r=None,d=None,scale=None,
        *args,**kws
    ):
        super().__init__(*args,**kws)

        if maps_dir is None:
            pass
        elif isinstance(maps_dir,list):
            maps_dir = ' '.join(maps_dir).replace(',',' ').replace(';',' ').split()
        elif isinstance(maps_dir,str):
            maps_dir = maps_dir.replace(',',' ').replace(';',' ').split()
        else:
            print(f'Warning: wrong defined maps_dir: {maps_dir}: ignored')
            maps_dir = []
        if maps_dir:
            d = os.listdir(maps_dir[0]) if os.path.isdir(maps_dir[0]) else []
            d = [i for i in d if i.endswith('.map')]
            if len(maps_dir) > 1:
                vs = [os.path.splitext(os.path.splitext(i)[0])[1] for i in d]
                vs = [i.replace('.', ' ').strip() for i in vs]
                if len(vs) != len(set(vs)):
                    print('Warning: similar maps exist in `maps_dir`, unexpected results will happen')
                p = []
                for k in maps_dir[1:]:
                    for v in zip(d,vs):
                        if k in v:
                            p.append(v[0])      # use full name
                            break
                d = p   # alias
            self.mapfiles_dir = [(os.path.join(maps_dir[0],i),os.path.basename(i)) for i in d]
            if self.mapfiles_dir:
                print(f'Note: processing on maps_dir: {maps_dir[0]}')
                print('   >> map files: '+' '.join(d))
        else:
            self.mapfiles_dir = []

        self.mapfiles = self._norm_mapfiles(mapfiles)
        if self.mapfiles and self.mapfiles_dir:
            print('Warning: conflicts: self.mapfile & maps_dir')
        self.gen_mapshape = False if gen_mapshape is False else True
        self.r = 0.0    # atomic radii
        self.d = 1.0    # distance
        self.scale = scale

    def _norm_mapfiles(self,mapfiles):
        if isinstance(mapfiles,(list,tuple)):
            mf = [i for i in mapfiles if os.path.isfile(i)]
            if mf:
                mapfiles = mf
            else:
                mapfiles = [self.recfile_s+'.'+i+'.map' for i in mapfiles]
                mapfiles = [i for i in mapfiles if os.path.isfile(i)]
        elif isinstance(mapfiles,str):
            if os.path.isfile(mapfiles):
                mapfiles = [mapfiles, ]
            else:
                mapfiles = [self.recfile_s+'.'+mapfiles+'.map', ]
                mapfiles = [i for i in mapfiles if os.path.isfile(i)]
        else:
            mapfiles = []

    def run(self,mapfiles=None,fillval=None):
        if not self.ligxyz:
            print('Fatal: ligand `-l` is not specified')
            return
        mapfiles = self._norm_mapfiles(mapfiles)
        if not mapfiles:
            if self.mapfiles_dir:
                mapfiles = self.mapfiles_dir
            else:
                mapfiles = self.mapfiles
        if not mapfiles:
            print('Warning: no valid map files is defined')
            return
        fn = ReadMap()
        for m in mapfiles:
            if isinstance(m,str):
                o = False
                b = m
            else:
                o = True
                m, b = m[0],m[1]
            print(f'Note: processing mapfile: {m}')
            data = fn.read_map_from_autodock(mapfile=m)
            if not data: continue
            overlaps,left = fn.calc_indexs(self.ligxyz)
            fn.fill_values(data,left,fillval,self.scale)
            fn.save_autodock_map_from_data(data,b,overwrite=o)

            if self.gen_mapshape:
                self._gen_shape(fn.center,fn.npts,fn.space,overlaps,m,'mapshape-')

    def _gen_shape(self,center,npts,space,overlaps,filename,head):
        xmin = center[0] - npts[0]/2*space
        ymin = center[1] - npts[1]/2*space
        zmin = center[2] - npts[2]/2*space
        fxyz = []
        for v in overlaps:
            x = v[0]*space + xmin
            y = v[1]*space + ymin
            z = v[2]*space + zmin
            fxyz.append((x,y,z))
        bm = os.path.basename(filename)
        new = xnewfile(bm+'.xyz',head=head)
        print(f'Note: shape file saved to: {new}')
        with open(new,'wt') as f:
            f.write(f'{len(fxyz)}\nGrid\n')
            for i in fxyz:
                f.write('X   {:.3f}   {:.3f}   {:.3f}\n'.format(*i))

    def gen_shape_on_ligfile(self,ligfile=None):
        if ligfile:
            super.__init__(ligfile=ligfile)
        if not self._flp:
            print('Warning: no valid ligand file is defined')
            return
        npts = list(map(int,self.npts_s.split()))
        space = 0.375
        center = self.gridcenter
        rd = self.r + self.d
        overlaps = set()
        xmin = center[0] - npts[0]/2*space
        ymin = center[1] - npts[1]/2*space
        zmin = center[2] - npts[2]/2*space
        for c in self.ligxyz:
            ox = (int((c[0]-rd-xmin)/space),int((c[0]+rd-xmin)/space))
            oy = (int((c[1]-rd-ymin)/space),int((c[1]+rd-ymin)/space))
            oz = (int((c[2]-rd-zmin)/space),int((c[2]+rd-zmin)/space))
            for i in range(ox[0],min(ox[1]+1,npts[0]+1)):
                for j in range(oy[0],min(oy[1]+1,npts[1]+1)):
                    for k in range(oz[0],min(oz[1]+1,npts[2]+1)):
                        overlaps.add((i,j,k))
        self._gen_shape(center,npts,space,overlaps,self.ligfile,'gridshape-')


class InteractiveRoot:
    def __init__(self,file=None,*args,**kws):
        self.nice = True
        # openbabel must be used to read file
        if file and os.path.isfile(file):
            if file.endswith('.pdbqt'):
                net = next(pybel.readfile('pdbqt',file))
                self.mol = Chem.MolFromPDBBlock(net.write('pdb'))
            elif file.endswith('.pdb'):
                net = next(pybel.readfile('pdb',file))
                self.mol = Chem.MolFromPDBBlock(net.write('pdb'))
            else:
                print(f'Warning: InteractiveRoot: not support: {file}')
                self.nice = False
        else:
            print(f'Warning: InteractiveRoot: not a file {file}')
            self.nice = False
        if self.mol:
            self.newmol = Chem.Mol(self.mol)
            AllChem.Compute2DCoords(self.newmol)
        else:
            self.nice = False

    def draw(self,fout=None,bool_draw=True):
        if not self.nice or not fout: return

        # [
        #   ['XZMARK', 'ROOT',   element="ad_atomtype", name="pdb_atomname", coords="x_y_z"], ...
        #   ['XZMARK', 'BRANCH', element="ad_atomtype", name="pdb_atomname", coords="x_y_z"], ...
        # ]
        marklines = []
        for i in fout.split('\n'):
            if i.startswith('XZMARK:'):
                marklines.append(i.replace(':',' ').replace(',',' ').split())

        rootndx = []
        branchndx = []
        for u in marklines:
            if len(u) == 5:
                elem = u[2].split('=')[1].strip()
                if elem.lower() != 'h':
                    name = u[3].split('=')[1].strip()
                    coor = list(map(float,u[4].split('=')[1].split('_')))
                    rootndx.append((elem,name,*coor))
            elif len(u) == 8:
                elem1 = u[2].split('=')[1].strip()
                name1 = u[3].split('=')[1].strip()
                coor1 = list(map(float,u[4].split('=')[1].split('_')))
                elem2 = u[5].split('=')[1].strip()
                name2 = u[6].split('=')[1].strip()
                coor2 = list(map(float,u[7].split('=')[1].split('_')))
                branchndx.append((elem1,name1,*coor1))
                branchndx.append((elem2,name2,*coor2))

        hitatoms = []
        conf = self.mol.GetConformer()
        for g in rootndx+branchndx:
            for i,c in enumerate(conf.GetPositions()):
                if abs(g[2]-c[0])<0.002 and abs(g[3]-c[1])<0.002 and abs(g[4]-c[2])<0.002:
                    hitatoms.append(i)
                    break

        n = len(rootndx)
        rootatoms = hitatoms[:n]
        rcx = round(sum([i[2] for i in rootndx])/n, 3)
        rcy = round(sum([i[3] for i in rootndx])/n, 3)
        rcz = round(sum([i[4] for i in rootndx])/n, 3)
        print(f'\nNote: processed ligand root atoms center: {rcx} {rcy} {rcz}')

        if not bool_draw:
            return rootatoms, (rcx,rcy,rcz)

        hitbonds = []
        for i in range(0,len(branchndx),2):
            a1 = hitatoms[n+i]
            a2 = hitatoms[n+i+1]
            hitbonds.append(self.mol.GetBondBetweenAtoms(a1,a2).GetIdx())

        for i,a in enumerate(self.newmol.GetAtoms()):
            a.SetProp("atomNote", str(a.GetIdx()+1))

        rootcolor = dict([(i,(0.0,1.0,0.25)) for i in rootatoms])

        p = rdMolDraw2D.MolDraw2DCairo(1000,1000)
        p.drawOptions().addAtomIndices = False
        rdMolDraw2D.PrepareAndDrawMolecule(
            p, self.newmol, highlightAtoms=rootatoms, highlightAtomColors=rootcolor,
            highlightBonds=hitbonds,
        )
        p.FinishDrawing()
        i = 1
        while True:
            file = 'root-picker-'+str(i)+'.png'
            if not os.path.isfile(file):
                break
            i += 1
        p.WriteDrawingText(file)
        print(f'Note: check file for inspection: {file}')
        print('Note: atoms in green means ROOT, bonds in red are Rotatable')
        print('Note: use the number to indicate New ROOT, otherwise ENTER for continuing')
        return [rcx,rcy,rcz]


def parsecmd():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""
XAutoDock {VERSION}. General Usage:

1. generate all needed files
>>> xdock.py -l ligand.pdb -r receptor.pdb -gb -gl -gr -gg -gd -gs

2. meaningful when `-l` specified
>>> xdock.py -l ligand.pdb -gl -gb [-gg] [-gs] [-c 0.0 0.0 0.0]

3. meaningful when `-r` specified
>>> xdock.py -r receptor.pdb -gr [-gd]

4. special cases, default values will be used for parameter files even no inputs
>>> xdock.py -gg -gd -gb

5. to modify exist map files
# use receptor's basename locate maps
>>> xdock.py -l ligand.pdb -r receptor.pdb -ms A C [--gen-mapshape]
# use full mapfiles name
>>> xdock.py -l ligand.pdb -ms receptor.A.map receptor.C.map [--gen-mapshape]
# use a folder
>>> xdock.py -l ligand.pdb -md /path/to/dir-maps [--gen-mapshape]
## suppose /path/to/dir-maps have map files: `rep.A.map`, `rep.C.map`, `rep.N.map`
#### modify all maps
>>> xdock.py -l ligand.pdb -md /path/to/dir-maps [--gen-mapshape]
#### modify `A` and `N` maps, use full name
>>> xdock.py -l ligand.pdb -md /path/to/dir-maps rep.A.map rep.N.map [--gen-mapshape]
#### or, alternatively, prefix and suffix will be auto added
>>> xdock.py -l ligand.pdb -md /path/to/dir-maps A N [--gen-mapshape]

6. centroid ligand first, then generate its shape
>>> xdock.py -l ligand.pdb -ms receptor.A.map -c 0.0 0.0 0.0

7. generate ligand shape
>>> xdock.py -l ligand.pdb -gs
""",
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-l','--ligfile',
        dest='ligfile',
        type=str,
        help='input ligand file'
    )
    parser.add_argument(
        '-r','--recfile',
        dest='recfile',
        type=str,
        help='input receptor file'
    )
    parser.add_argument(
        '-c','--centroid-at',
        dest='centroid_at',
        metavar='v',
        nargs=3,
        type=float,
        help='valid when `-l` is used, centroid `ligfile` at defined center before any generations'
    )
    parser.add_argument(
        '--npts',
        dest='npts',
        metavar='v',
        nargs=3,
        type=int,
        help='number of grid points'
    )
    parser.add_argument(
        '-gb','--gen-gridbox',
        dest='gen_box',
        action='store_true',
        help='generate grid box based on input ligand'
    )
    parser.add_argument(
        '-gl','--gen-ligand-pdbqt',
        dest='gen_lig',
        action='store_true',
        help='parse input ligand file and then generate pdbqt file'
    )
    parser.add_argument(
        '-glm','--gl-mode',
        dest='gen_lig_mode',
        default=0,
        type=int,
        choices={0,1,2},
        help='valid when `-gl`, mode: 0(interactive), 1(max atoms), 2(least atoms)'
    )
    parser.add_argument(
        '-gr','--gen-receptor-pdbqt',
        dest='gen_rec',
        action='store_true',
        help='parse input receptor file and then generate pdbqt file'
    )
    parser.add_argument(
        '-gg','--gen-gpf',
        dest='gen_gpf',
        action='store_true',
        help='generate gpf file'
    )
    parser.add_argument(
        '-gd','--gen-dpf',
        dest='gen_dpf',
        action='store_true',
        help='generate dpf file'
    )
    parser.add_argument(
        '-gs','--gen-ligand-gridshape',
        dest='gen_ligand_gridshape',
        action='store_true',
        help='generate ligand grid shape file'
    )
    parser.add_argument(
        '-ms','--modify-maps-types',
        dest='maps_types',
        metavar='m',
        nargs='+',
        help='ongoing modified maps'
    )
    parser.add_argument(
        '-md','--modify-maps-dir',
        metavar='v',
        nargs='+',
        dest='maps_dir',
        help=(
            'process maps only in this folder, conflict with `-ms`, the first value should '
            'be the path of folder, the following values are maps will be modified'
        )
    )
    parser.add_argument(
        '--gen-mapshape',
        dest='gen_mapshape',
        action='store_true',
        help='valid when `-ms` or `-md` is used, generate grid shape file',
    )
    parser.add_argument(
        '-i','--fill-value',
        dest='fillval',
        metavar='v',
        type=float,
        help='valid when `-ms` or `-md` is used, fill values'
    )
    parser.add_argument(
        '-s','--scale',
        dest='scale',
        metavar='v',
        type=float,
        help='valid when `-ms` or `-md` is used, scale value on grid point'
    )
    parser.add_argument(
        '--verbose',
        dest='verbose',
        action='store_true',
        help='show more info'
    )
    parser.add_argument(
        '--features',
        dest='feature',
        action='store_true',
        help='show development features'
    )
    if len(sys.argv) == 1:
        parser.print_help()
        return
    w = parser.parse_args(sys.argv[1:])
    if w.feature:
        for i in FEATURES:
            print(i)
        return

    if w.maps_dir and w.maps_types:
        print('Fatal: conflict argument: --modify-maps-dir & --modify-maps-types')
        return

    if w.npts:
        new = [i//2*2 for i in w.npts]
        if w.npts != new:
            print(f'Warning: npts: only even integers are accepted: {w.npts}')
            print(f'     ->  converted to: {new}')
            w.npts = new
    xp = LigandShape(**vars(w))

    if w.gen_lig:       # has to be put in first
        xp.gen_pdbqt_for_ligand()
    if w.gen_rec:
        xp.gen_pdbqt_for_protein()
    if w.gen_box:
        xp.gen_gridbox_pdb()
    if w.gen_gpf:
        xp.gen_gpf()
    if w.gen_dpf:
        xp.gen_dpf()
    if w.gen_ligand_gridshape:
        xp.gen_shape_on_ligfile()

    if w.maps_types:
        alls = [i.lower() for i in xp.types_lig_s.split()]
        news = [i for i in w.maps_types if i.lower() in alls]
        if news != w.maps_types:
            print(f'Warning: modify maps types: not correspondent: {w.maps_types}')
            keys = xp.types_lig_s.split()
            news = [keys[alls.index(i)] for i in news]
            print(f'     ->  converted to: {news}')
            w.maps_types = news
        xp.run(mapfiles=w.maps_types,fillval=w.fillval)
    elif w.maps_dir:
        if xp.mapfiles_dir:
            xp.run(fillval=w.fillval)


if __name__ == '__main__':
    ANCHOR = None
    parsecmd()



