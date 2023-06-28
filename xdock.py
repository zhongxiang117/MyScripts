#!/usr/bin/env python3


from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4LigandPreparation, AD4ReceptorPreparation

import os
import sys
import shutil
import argparse

FEATURES = [
    'version 0.1.0  : XAutoDock',
    'version 0.2.0  : add argparse',
    'version 0.3.0  : add option for modification of maps',
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
gridfld {gridfld}
{maps}
move {move}
about {about}
tran0 random                         # initial coordinates/A or random
quaternion0 random                   # initial orientation
dihe0 random                         # initial dihedrals (relative) or random
torsdof 8                            # torsional degrees of freedom
rmstol 2.0                           # cluster_tolerance/A
extnrg 1000.0                        # external grid energy
e0max 0.0 10000                      # max initial energy; max number of retries
ga_pop_size 150                      # number of individuals in population
ga_num_evals 25000                   # maximum number of energy evaluations
ga_num_generations 2700              # maximum number of generations
ga_elitism 1                         # number of top individuals to survive to next generation
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


def xnewfile(file,head=None):
    if not os.path.isfile(file): return file
    if not head: head = 'new-'
    i = 1
    while True:
        new = head + str(i) + '-' + file
        if not os.path.isfile(new): break
        i += 1
    return new


class XADTPrep:
    def __init__(
        self,recfile=None,ligfile=None,
        npts=None,gridcenter=None,ligand_types=None,receptor_types=None,
        verbose=None,*args,**kws
    ):
        self.verbose = True if verbose is True else False
        self.recfile = recfile if recfile else 'receptor'
        # be aware, basename should be used!
        self.recfile_s = os.path.basename(os.path.splitext(self.recfile)[0])
        self.recfile_s = xnewfile(self.recfile_s)

        self.ligfile = ligfile if ligfile else 'ligand'
        self.ligfile_s = os.path.basename(os.path.splitext(self.ligfile)[0])
        self.ligfile_s = xnewfile(self.ligfile_s)
        self.gridcenter = [0.0, 0.0, 0.0]
        self.ligxyz = None
        self._flp = None
        self._gen_pdbqt_for_ligand()

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
            bo = True
        if bo:
            if info:
                print(f'Warning: wrong defined: {info}: {l}')
            else:
                print(f'Warning: wrong defined: {l}')
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
        elif verbose: print("No change in atomic coordinates")

        # if mol.returnCode != 0:
        #     print(f'ligand return msg: {mol.returnMsg}')
        self.ligxyz = mol.allAtoms.coords
        # calc gridcenter
        sx = sum([i[0] for i in self.ligxyz])
        sy = sum([i[1] for i in self.ligxyz])
        sz = sum([i[2] for i in self.ligxyz])
        n = len(mol.allAtoms)
        self.gridcenter = [round(sx/n,3), round(sy/n,3), round(sz/n,3)]

    def gen_pdbqt_for_ligand(self,file=None,verbose=None):
        if file:
            self._gen_pdbqt_for_ligand(self,file,verbose)
        else:
            if not self._flp: return
        print(f'Note: Writing ligand pdbqt to: {self.ligfile_s}.pdbqt')
        self._flp.write(self.ligfile_s+'.pdbqt')

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


class ReadMap:
    """Read Map Files
    
    data format: [z[y[x[float]]]]
    """
    def __init__(
        self,mapfile=None,npts=None,*args,**kws
    ):
        self.mapfile = mapfile
        self.npts = npts
        self._pars = []
        self.center = [0.0, 0.0, 0.0]
        self.space = 0.375

    def read_map_from_autodock(self,mapfile=None,npts=None):
        if not mapfile: mapfile = self.mapfile
        self.mapfile = mapfile      # reset current mapfile name
        if not npts: npts = self.npts
        if not os.path.isfile(mapfile):
            print(f'Error: map: not a file: {mapfile}')
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
            new = xnewfile(filename,'bak-')
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
        print(f'Note: map file saved to {filename}')

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

    def fill_values(self,data,indexes,fillval=None):
        if not fillval: fillval = 0.0
        for x in indexes[0]:
            for y in indexes[1]:
                for z in indexes[2]:
                    data[x][y][z] = fillval

    def calc_indexs(self,span_ranges):
        bo = False
        if not span_ranges:
            bo = True
        elif len(span_ranges) == 3:
            if not span_ranges[0] or not span_ranges[1] or not span_ranges[2]:
                bo = True
        else:
            bo = True
        if bo:
            print('Warning: index: wrong inputs: span_range')
            return []

        """Save for future use
        xmin = self.center[0] - self.npts[0]/2*self.space
        xmax = self.center[0] + self.npts[0]/2*self.space
        ymin = self.center[1] - self.npts[1]/2*self.space
        ymax = self.center[1] + self.npts[1]/2*self.space
        zmin = self.center[2] - self.npts[2]/2*self.space
        zmax = self.center[2] + self.npts[2]/2*self.space
        if span_ranges[0][0][0]<xmin or span_ranges[0][-1][1]>xmax:
            print('Warning: x axis: not correspondent: grid/ligand')
            return []
        if span_ranges[1][0][0]<ymin or span_ranges[1][-1][1]>ymax:
            print('Warning: y axis: not correspondent: grid/ligand')
            return []
        if span_ranges[2][0][0]<zmin or span_ranges[2][-1][1]>zmax:
            print('Warning: z axis: not correspondent: grid/ligand')
            return []
        """
        
        xi = self._calc_index(span_ranges[0],'x')
        yi = self._calc_index(span_ranges[1],'y')
        zi = self._calc_index(span_ranges[2],'z')
        return [xi,yi,zi]


    def _calc_index(self,spans,type_='x'):
        if type_ == 'x':
            t = 0
        elif type_ == 'y':
            t = 1
        else:
            t = 2
        smin = a = self.center[t] - self.npts[t]/2*self.space
        b = self.center[t] + self.npts[t]/2*self.space

        voids = []
        for v in spans:
            if v[0] > a:
                voids.append((a,v[0]))
            a = v[1]
        if a < b:
            voids.append((a,b))

        if not voids:
            return []

        ndxlist = []
        c = 0
        for i in range(self.npts[t]+1):
            v = smin + i * self.space
            while c < len(voids):
                if v >= voids[c][0] and v <= voids[c][1]:
                    ndxlist.append(i)
                    break
                c += 1
        return ndxlist


class LigandShape(XADTPrep):
    def __init__(self,mapfiles=None,r=None,d=None,*args,**kws):
        super().__init__(*args,**kws)
        self.mapfiles = self._norm_mapfiles(mapfiles)
        self.r = 0.0    # atomic radii
        self.d = 0.0    # distance
    
    def _norm_mapfiles(self,mapfiles):
        if isinstance(mapfiles,(list,tuple)):
            mapfiles = [self.recfile_s+'.'+i+'.map' for i in mapfiles]
        elif isinstance(mapfiles,str):
            mapfiles = [self.recfile_s+'.'+mapfiles+'.map', ]
        else:
            mapfiles = []
        return mapfiles

    def run(self,mapfiles=None,fillval=None):
        mapfiles = self._norm_mapfiles(mapfiles)
        if not mapfiles: mapfiles = self.mapfiles
        if not mapfiles: return
        fn = ReadMap()
        for m in mapfiles:
            data = fn.read_map_from_autodock(mapfile=m)
            if not data: continue
            span_ranges = self.calc_span_ranges()
            indexes = fn.calc_indexs(span_ranges)
            fn.fill_values(data,indexes,fillval)
            fn.save_autodock_map_from_data(data,m)

    def _span_ranges(self,l,d=None):
        if not d: d = self.d
        if not l:
            return []
        elif len(l) == 1:
            return [(l[0]-self.r-self.d, l[0]+self.r+self.d), ]
        rd = self.r + self.d
        ls = sorted(l)
        a = ls[0] - rd
        b = ls[0] + rd
        rst = []
        for v in ls[1:]:
            if b < v-rd:
                rst.append((a,b))
                a = v-rd
            b = v+rd
        rst.append((a,b))
        return rst

    def calc_span_ranges(self,d=None):
        lx = [i[0] for i in self.ligxyz]
        xs = self._span_ranges(lx,d)
        ly = [i[1] for i in self.ligxyz]
        ys = self._span_ranges(ly,d)
        lz = [i[2] for i in self.ligxyz]
        zs = self._span_ranges(lz,d)
        return [xs, ys, zs]


def parsecmd():
    parser = argparse.ArgumentParser(
        description='XAutoDock '+VERSION,
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
        '--npts',
        dest='npts',
        nargs=3,
        type=int,
        help='number of grid points'
    )
    parser.add_argument(
        '-gb','--gen-gridbox',
        dest='g_box',
        action='store_true',
        help='generate grid box based on input ligand'
    )
    parser.add_argument(
        '-gl','--gen-ligand-pdbqt',
        dest='g_lig',
        action='store_true',
        help='parse input ligand file and then generate pdbqt file'
    )
    parser.add_argument(
        '-gr','--gen-receptor-pdbqt',
        dest='g_rec',
        action='store_true',
        help='parse input receptor file and then generate pdbqt file'
    )
    parser.add_argument(
        '-gg','--gen-gpf',
        dest='g_gpf',
        action='store_true',
        help='generate gpf file'
    )
    parser.add_argument(
        '-gd','--gen-dpf',
        dest='g_dpf',
        action='store_true',
        help='generate dpf file'
    )
    parser.add_argument(
        '-ms','--modify-maps',
        dest='ms',
        nargs='+',
        help='generate ligand-shaped maps for inputs'
    )
    parser.add_argument(
        '-i','--fill-value',
        dest='fillval',
        type=float,
        help='valid when `-ms` is used, fill values'
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
    if w.npts:
        new = [i//2*2 for i in w.npts]
        if w.npts != new:
            print(f'Warning: npts: only even integers are accepted: {w.npts}')
            print(f'     ->  converted to: {new}')
            w.npts = new
    xp = LigandShape(recfile=w.recfile,ligfile=w.ligfile,npts=w.npts,verbose=w.verbose)

    if w.g_lig:
        xp.gen_pdbqt_for_ligand()
    if w.g_rec:
        xp.gen_pdbqt_for_protein()
    if w.g_box:
        xp.gen_gridbox_pdb()
    if w.g_gpf:
        xp.gen_gpf()
    if w.g_dpf:
        xp.gen_dpf()

    if w.ms:
        alls = [i.lower() for i in xp.types_lig_s.split()]
        news = [i for i in w.ms if i.lower() in alls]
        if news != w.ms:
            print(f'Warning: modify maps: not correspondent: {w.ms}')
            keys = xp.types_lig_s.split()
            news = [keys[alls.index(i)] for i in news]
            print(f'     ->  converted to: {news}')
            w.ms = news
        xp.run(mapfiles=w.ms,fillval=w.fillval)



if __name__ == '__main__':
    parsecmd()



