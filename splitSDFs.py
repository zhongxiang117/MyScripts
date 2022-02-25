#!/usr/bin/env python3

try:
    # version 3
    from openbabel import openbabel
except ImportError:
    # version 2
    import openbabel
import numpy as np

import os
import sys
import time
import argparse

FEATURES = [
    'version 0.1.0  : Split Database File: Feb 18, 2022',
    'version 0.2.0  : add generate conformers'
    'version 0.3.0  : add nonbonded fragment, usually, salts, filtration',
    'version 0.4.0  : add bool_nofrags_check',
    'version 0.5.0  : add more controls on RotorSearch conformer generations method',
    'version 0.6.0  : add RMSD calculation',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
_ffs = openbabel.vectorString()
openbabel.OBPlugin.ListAsVector('forcefields',None,_ffs)
SUPPORTED_FFS = [i.split()[0].lower() for i in _ffs]
OBConv = openbabel.OBConversion()
_fms = OBConv.GetSupportedInputFormat()
SUPPORTED_FORMATS = [i.split()[0].lower() for i in _fms]


class RMSDClosedFormOBMols:
    """Calculate RMSD for obmols by using closedForm algorithm

    Args:
        obmols(List[obmol, ]): the first obmol will be used for reference
        vl : 2D n*3f : List[[float, float, float]] : Left input matrix
        vr : 2D n*3f : List[[float, float, float]] : Right input matrix
        vl = BestFit * vr
        Left set matrix is used for reference

    Methods:
        centroid
        calc_N
        calc_left_rotM
        calc_bestfit

    Reference:
        Closed-form solution of absolute orientation using unit quaternions
        Berthold K. P. Horn
        J. Opt. Soc. Am. A 4, 629-642 (1987)
    """
    def __init__(self,obmols=None,vl=None,vr=None,*args,**kwargs):
        self.obmols = obmols if obmols else []
        self.vl = vl
        self.vr = vr

    def run(self,obmols=None):
        if obmols is None: obmols = self.obmols
        if len(obmols) < 2: return []

        # reference mol
        vl = [[a.GetX(),a.GetY(),a.GetZ()] for a in openbabel.OBMolAtomIter(obmols[0])]
        cvl,tl = self.centroid(vl)

        rmsdlist = [0.0, ]
        for mol in obmols[1:]:
            vr = [[a.GetX(),a.GetY(),a.GetZ()] for a in openbabel.OBMolAtomIter(mol)]
            fit = self.calc_bestfit(cvl,vr,centroid_vl=False)
            rmsdlist.append(self.calc_rmsd(cvl, fit))
        return rmsdlist

    def calc_rmsd(self,vl,vr):
        """calculate overall RMSD for input vl and vr"""
        if len(vl) != len(vr): return 0.0
        rmsd = 0.0
        for i in range(len(vl)):
            t = [vl[i][j]-vr[i][j] for j in range(3)]
            rmsd += sum([k*k for k in t])
        return pow(rmsd,0.5)

    def calc_bestfit(self,vl=None,vr=None,centroid=True,centroid_vl=True,centroid_vr=True):
        """calc best fit structure
        both vl(reference) and vr(target) can be reused
        centroid: whether centroid results or not
        """
        if vl is None: vl = self.vl
        if vr is None: vr = self.vr
        if centroid_vl:
            cvl,tl = self.centroid(vl)
        else:
            cvl = vl
            tl = [0.0, 0.0, 0.0]
        if centroid_vr:
            cvr,tr = self.centroid(vr)
        else:
            cvr = vr
            tr = [0.0, 0.0, 0.0]
        M = self.calc_left_rotM(cvl,cvr)
        fit = []
        for v in cvr:
            rx = v[0]*M[0][0] + v[1]*M[1][0] + v[2]*M[2][0]
            ry = v[0]*M[0][1] + v[1]*M[1][1] + v[2]*M[2][1]
            rz = v[0]*M[0][2] + v[1]*M[1][2] + v[2]*M[2][2]
            fit.append([rx,ry,rz])
        if not centroid:
            for i in range(len(fit)):
                fit[i][0] += tr[0]
                fit[i][1] += tr[1]
                fit[i][2] += tr[2]
        return fit

    def centroid(self,v):
        """calc centroid vector and translation vector"""
        x = [i[0] for i in v]
        y = [i[1] for i in v]
        z = [i[2] for i in v]
        t = len(v)
        ax = sum(x) / t
        ay = sum(y) / t
        az = sum(z) / t
        # care, do not use in-place operation
        cv = [[0.0 for i in range(3)] for j in range(len(v))]
        for i in range(t):
            cv[i][0] = v[i][0] - ax
            cv[i][1] = v[i][1] - ay
            cv[i][2] = v[i][2] - az
        return cv,(ax,ay,az)

    def calc_N(self,vl,vr):
        """calc 4x4 real symmetric N matrix"""
        XxYx = 0.0
        XxYy = 0.0
        XxYz = 0.0
        XyYx = 0.0
        XyYy = 0.0
        XyYz = 0.0
        XzYx = 0.0
        XzYy = 0.0
        XzYz = 0.0
        # careful of the sequence: X-r, Y-l
        # for rotation l = Rr
        for i,p in enumerate(vl):
            XxYx += p[0] * vr[i][0]
            XxYy += p[0] * vr[i][1]
            XxYz += p[0] * vr[i][2]

            XyYx += p[1] * vr[i][0]
            XyYy += p[1] * vr[i][1]
            XyYz += p[1] * vr[i][2]

            XzYx += p[2] * vr[i][0]
            XzYy += p[2] * vr[i][1]
            XzYz += p[2] * vr[i][2]

        N = [[0.0, 0.0, 0.0, 0.0] for i in range(4)]

        N[0][0] = XxYx + XyYy + XzYz
        N[0][1] = XyYz - XzYy
        N[0][2] = XzYx - XxYz
        N[0][3] = XxYy - XyYx

        N[1][0] = N[0][1]
        N[1][1] = XxYx - XyYy - XzYz
        N[1][2] = XxYy + XyYx
        N[1][3] = XzYx + XxYz

        N[2][0] = N[0][2]
        N[2][1] = N[1][2]
        N[2][2] = -XxYx + XyYy - XzYz
        N[2][3] = XyYz + XzYy

        N[3][0] = N[0][3]
        N[3][1] = N[1][3]
        N[3][2] = N[2][3]
        N[3][3] = -XxYx - XyYy + XzYz

        return N

    def calc_left_rotM(self,vl,vr):
        """calc rotation matrix for vr*M,
        quaternion is got from the vector which is corresponding to
        largest positive eigenvalues

        M  : 2D 3*3f : rotation matrix for vr, vr*M, in element-wise operation
        """
        N = self.calc_N(vl,vr)
        values,vectors = np.linalg.eig(N)
        ndx = np.where(values == max(values))
        ndx = ndx[0][0]
        # For numpy, eigenvectors are correspondingly put in column
        # note, this vector has already been normalized
        V = vectors[:,ndx]
        M = [[0.0, 0.0, 0.0] for i in range(3)]

        M[0][0] = 1 - 2 * (V[2]*V[2] + V[3]*V[3])
        M[0][1] = 2 * (V[1]*V[2] - V[3]*V[0])
        M[0][2] = 2 * (V[0]*V[3] + V[2]*V[0])

        M[1][0] = 2 * (V[1]*V[2] + V[3]*V[0])
        M[1][1] = 1 - 2 * (V[1]*V[1] + V[3]*V[3])
        M[1][2] = 2 * (V[2]*V[3] - V[1]*V[0])

        M[2][0] = 2 * (V[1]*V[3] - V[2]*V[0])
        M[2][1] = 2 * (V[2]*V[3] + V[1]*V[0])
        M[2][2] = 1 - 2 * (V[1]*V[1] + V[2]*V[2])

        return M


class ReadRawFile:
    """read data from input file

    Attributes:
        filename(str): input file name
        fileobmols(List[openbabel.OBMol, ]):
        fileobfrags(List[tuple(openbabel.OBMol,SMILE), ]):
        filerftimes(List[float]):
        runtime(float):
    """
    def __init__(self,
                filename=None,
                ignore_h=None,
                add_h=None,
                gen3d=None,
                force_field=None,
                bool_nofrags_check=None,
                rs_confs_totnum=None,
                rs_confs_selnum=None,
                bool_confs_select_all=None,
                bool_calc_rmsd=None,
                *args,
                **kwargs):
        self.filename = filename
        # ignore_h is in precedence with add_h
        self.ignore_h = True if ignore_h is True else False
        self.add_h = True if add_h is True else False
        self.gen3d = True if gen3d is True else False
        self.bool_nofrags_check = True if bool_nofrags_check is True else False
        self.rs_confs_totnum = rs_confs_totnum if rs_confs_totnum else 20
        self.rs_confs_selnum = rs_confs_selnum if rs_confs_selnum else 3
        self.bool_confs_select_all = True if bool_confs_select_all is True else False
        self.bool_calc_rmsd = True if bool_calc_rmsd is True else False
        if self.bool_confs_select_all:
            self.rs_confs_selnum = self.rs_confs_totnum
        if force_field is not None:
            self.force_field = force_field.lower()
            if self.force_field not in SUPPORTED_FFS:
                print('Warning: not support: force_field: ',self.force_field)
                print('Warning: revert force_field to default value')
                self.force_field = None
        else:
            self.force_field = None
        # openbabel miscellany, universally set to 'smi'
        OBConv.SetOutFormat('smi')
        # hidden parameter for conformer generations
        # filter out any energy difference less than this value, kcal/mol
        self._energy_tol = 1.0

    def run(self,filename=None):
        curtime = time.time()
        if filename is None: filename = self.filename
        self.fileobmols, self.fileobfrags, self.filerftimes = self.readfile(filename)
        self.runtime = time.time() - curtime

    def readfile(self,file):
        """read data from file based on extension"""
        bo = False
        if '.' in file:
            ext = file.split('.')[-1].lower()
            if ext not in SUPPORTED_FORMATS: bo = True
        else:
            bo = True
        if bo:
            print('Warning: file not support: ',file)
            return [], [], []

        OBConv.SetInFormat(ext)

        if self.force_field:
            setff = openbabel.OBForceField.FindForceField(self.force_field)
        else:
            setff = openbabel.OBForceField.FindForceField('mmff94')

        obmollist = []
        obfraglist = []
        timelist = []
        obmol = openbabel.OBMol()
        now = time.time()
        bo = OBConv.ReadFile(obmol,file)
        while bo:
            if not self.bool_nofrags_check:
                # first, check whether is fragment
                tmpbo, smi = self.is_fragment(obmol)
                if tmpbo:
                    obfraglist.append((obmol,smi))
                    # continue reading
                    now = time.time()
                    obmol = openbabel.OBMol()
                    bo = OBConv.Read(obmol)
                    continue

            # all work well due to the python implicit index
            myobmols = [obmol, ]
            # optimization should be executed in first
            if self.gen3d:
                obmol.AddHydrogens()
                if setff.Setup(obmol):
                    myobmols = self._f_gen3d_many_rs(setff,obmol)
            if self.ignore_h:
                for i in myobmols: i.DeleteHydrogens()
            elif self.add_h:
                if not self.gen3d:
                    # means no myobmols at all, care on the logical
                    if setff.Setup(obmol):
                        myobmols = self._f_gen3d_many_rs(setff,obmol)
            obmollist.extend(myobmols)
            moment = time.time()
            s = (moment-now) / len(obmollist)
            timelist.extend([s for i in range(len(obmollist))])
            # continue reading
            now = moment
            obmol = openbabel.OBMol()
            bo = OBConv.Read(obmol)

        return obmollist, obfraglist, timelist

    def is_fragment(self,obmol):
        """determine whether input is salt"""
        # format: SMILES\tTitle\n
        smi = OBConv.WriteString(obmol).split()[0]
        if '.' in smi:
            return True, smi
        return False, smi

    def _f_gen3d_many_rs(self,setff,obmol):
        """generate conformers based on torsion rotation search"""
        self.add_h = True
        setff.EnableCutOff(True)
        setff.SetVDWCutOff(10.0)
        setff.SetElectrostaticCutOff(20.0)
        setff.SetUpdateFrequency(10)
        #setff.ConjugateGradients(250,10**(-4))
        setff.SteepestDescent(500)

        #setff.WeightedRotorSearch(250,10)
        #setff.ConjugateGradients(500,10**(-6))

        obmollist = []
        enelist = []
        setff.RandomRotorSearchInitialize(self.rs_confs_totnum)
        while setff.RandomRotorSearchNextConformer(500):
            #setff.ConjugateGradients(500,10**(-6))
            ene = setff.Energy()
            enelist.append(ene)
            obtmp = openbabel.OBMol(obmol)
            data = openbabel.OBPairData()
            data.SetAttribute('OBEnergy')
            data.SetValue('{:.4f}'.format(ene))
            obtmp.CloneData(data)
            setff.UpdateCoordinates(obtmp)
            obmollist.append(obtmp)

        if len(obmollist) <= 1:
            # means done on initialization or before first N steps
            setff.UpdateCoordinates(obmol)
            return [obmol, ]

        # sort obmollist from energy smallest to largest
        ndxlist = sorted(range(len(enelist)),key=lambda k: enelist[k])
        # do energy tolerance filtration
        filterlist = [ndxlist[0], ]
        for i in ndxlist[1:]:
            if abs(enelist[i]-enelist[filterlist[-1]]) > self._energy_tol:
                filterlist.append(i)
        # update
        obmollist = [obmollist[i] for i in filterlist]

        if len(obmollist) <= self.rs_confs_selnum:
            # consideration when bool_confs_select_all
            if len(obmollist) and self.bool_calc_rmsd:
                fn = RMSDClosedFormOBMols(obmols=obmollist)
                rmsdlist = fn.run()
                # implicitly linking
                for i,obmol in enumerate(obmollist):
                    data = openbabel.OBPairData()
                    data.SetAttribute('OBRMSD')
                    data.SetValue('{:.4f}'.format(rmsdlist[i]))
                    obmol.CloneData(data)
            return obmollist

        if self.rs_confs_selnum == 1:
            return [obmollist[0], ]

        if self.rs_confs_selnum == 2:
            return [obmollist[0], obmollist[-1]]

        dt = (len(obmollist) - 1) // (self.rs_confs_selnum - 1)
        reflist = list(range(0,len(obmollist),dt))
        if len(reflist) > self.rs_confs_selnum:
            reflist = reflist[:self.rs_confs_selnum]

        # update
        obmollist = [obmollist[i] for i in reflist]

        if self.bool_calc_rmsd:
            fn = RMSDClosedFormOBMols(obmols=obmollist)
            rmsdlist = fn.run()
            # implicitly linking
            for i,obmol in enumerate(obmollist):
                data = openbabel.OBPairData()
                data.SetAttribute('OBRMSD')
                data.SetValue('{:.4f}'.format(rmsdlist[i]))
                obmol.CloneData(data)

        return obmollist

    def _f_gen3d_many_cg(self,setff,obmol):
        """generate conformers based on conjugate gradients search"""
        self.add_h = True
        setff.EnableCutOff(True)
        setff.SetVDWCutOff(10.0)
        setff.SetElectrostaticCutOff(20.0)
        setff.SetUpdateFrequency(10)
        #setff.ConjugateGradients(250,10**(-4))
        setff.SteepestDescent(500)

        #setff.WeightedRotorSearch(250,10)
        #setff.ConjugateGradients(500,10**(-6))

        obmollist = []
        enelist = []
        setff.ConjugateGradientsInitialize(500,10**(-6))
        while setff.ConjugateGradientsTakeNSteps(5):
            ene = setff.Energy()
            enelist.append(ene)
            obtmp = openbabel.OBMol(obmol)
            data = openbabel.OBPairData()
            data.SetAttribute('OBEnergy')
            data.SetValue(str(ene))
            obtmp.CloneData(data)
            setff.UpdateCoordinates(obtmp)
            obmollist.append(obtmp)

        if not len(obmollist):
            # means done on initialization or before first N steps
            setff.UpdateCoordinates(obmol)
            return [obmol, ]

        if len(obmollist) <= 3:
            return obmollist

        ndxlist = sorted(range(len(enelist)),key=lambda k: enelist[k])
        mid = len(enelist) // 2
        return [obmollist[ndxlist[0]], obmollist[mid], obmollist[ndxlist[-1]]]


class SaveRawFile:
    """save final results"""
    def __init__(self,
                ftype=None,
                data=None,
                fname=None,
                bool_saveproinfo=None,
                bool_savesmiles=None,
                *args,
                **kwargs):
        self.ftype = 'sdf'
        if ftype:
            self.ftype = ftype.lower()
            if self.ftype not in SUPPORTED_FORMATS:
                print('Fatal: currently file format not support: ',self.ftype)
                self.ftype = 'sdf'
        if isinstance(data,list):
            self.data = data
        elif isinstance(data,dict):
            self.data = [data, ]
        else:
            print('Fatal: not data input')
            self.data = None
        if fname:
            self.fname = fname
        else:
            self.fname = 'artrocs'
        self.bool_saveproinfo = False if bool_saveproinfo is False else True
        self.bool_savesmiles = True if bool_savesmiles is True else False
        if self.data:
            self._myinfo = ''
            self.saveobmols()
            if self.bool_saveproinfo:
                self.saveproinfo()
            if self.bool_savesmiles:
                self.savesmiles()

    def savesmiles(self):
        filesmile = self._get_fnew(self.fname+'-smiles', 'txt')
        print('Note: ArtROCS salts info is saved to: '+filesmile+'.txt')
        with open(filesmile+'.txt','wt') as f:
            for d in self.data:
                f.write(d['smile'])
                f.write('\n')

    def saveproinfo(self):
        info = 'ArtROCS: ' + time.ctime() + '\n' + self._myinfo
        info += 'FileName                            molindex     ReadFileTimeCost(s)\n'
        for d in self.data:
            if not len(d): continue
            info += '{:34}  {:^6}        {:>.6f}\n'.format(d['filename'],d['molindex'],d['rftime'])
        info += '\n\n'
        # check existing files to avoid overwriting
        fileinfo = self._get_fnew(self.fname+'-info', 'txt')
        print('Note: ArtROCS processing info is saved to: '+fileinfo+'.txt')
        with open(fileinfo+'.txt','wt') as f: f.write(info)

    def saveobmols(self):
        # check existing files to avoid overwriting
        fileoverlap = self._get_fnew(self.fname, self.ftype)
        tmp = 'Note: result file is saved to: '+fileoverlap+'.'+self.ftype
        print(tmp)
        self._myinfo += tmp + '\n'
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat(self.ftype)
        obconv.WriteFile(self.data[0]['obmol'],fileoverlap+'.'+self.ftype)
        if len(self.data) > 1:
            for d in self.data[1:]:
                if d['obmol']:
                    obconv.Write(d['obmol'])
        obconv.CloseOutFile()

    def _get_fnew(self,fname,ftype):
        if not os.path.isfile(fname+'.'+ftype):
            return fname
        i = 0
        while True:
            i += 1
            if not os.path.isfile(fname+'-'+str(i)+'.'+ftype):
                break
        return fname+'-'+str(i)


class SplitDatabase:
    def __init__(self,
                files=None,
                splitnum=None,
                fname=None,
                *args,**kwargs):
        self.files = files
        self.splitnum = splitnum if splitnum else None

        totdata = []
        totfrags = []
        rf = ReadRawFile(*args,**kwargs)
        for f in self.files:
            rf.filename = f
            rf.run()
            for i,d in enumerate(rf.fileobmols):
                # iterate all data, to find the first valid one
                v = None
                for m in d.GetData().iterator():
                    k = m.GetAttribute()
                    if len(k.split()) == 1:
                        v = m.GetValue()
                        break
                if not v:
                    # not found, use title
                    v = obmol.GetTitle()
                totdata.append({
                    'filename'  : rf.filename,
                    'obmol'     : d,
                    'molindex'  : v,
                    'rftime'    : rf.filerftimes[i]
                })

            for i,d in enumerate(rf.fileobfrags):
                totfrags.append({
                    'filename'  : rf.filename,
                    'obmol'     : d[0],
                    'smile'     : d[1],
                })

        if self.splitnum is None:
            datalist = [totdata,]
        else:
            nmlist = list(range(0,len(totdata)+1,self.splitnum))
            if nmlist[-1] != len(totdata): nmlist.append(len(totdata))
            # care in the trick, index actually starts from 1
            datalist = [totdata[nmlist[i]:v] for i,v in enumerate(nmlist[1:])]
        for data in datalist:
            SaveRawFile(data=data,fname=fname,*args,**kwargs)

        # save fragment
        fname = fname +'-salts' if fname else 'artrocs-salts'
        SaveRawFile(data=totfrags,bool_saveproinfo=False,fname=fname,bool_savesmiles=True,*args,**kwargs)


def parsecmd():
    parser = argparse.ArgumentParser(
        description='ArtROCS',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--files',
        nargs='+',
        help='input files'
    )
    parser.add_argument(
        '--ignore_h',
        action='store_true',
        help='ignore hydrogens, highest priority'
    )
    parser.add_argument(
        '--add_h',
        action='store_true',
        help='add hydrogens for input files if needed'
    )
    parser.add_argument(
        '--gen3d',
        action='store_true',
        help='turn on generating 3D coordinates'
    )
    parser.add_argument(
        '-s','--splitnum',
        type=int,
        help='number to be split file',
    )
    parser.add_argument(
        '--force_field',
        help='valid on gen3d or add_h, supported force fields: '+'; '.join(SUPPORTED_FFS)
    )
    parser.add_argument(
        '--ftype',
        help='output file type, check by option: --show_supported_types',
    )
    parser.add_argument(
        '--show_supported_types',
        action='store_true',
        help='show supported file output types',
    )
    parser.add_argument(
        '--fname',
        help='output file basename'
    )
    parser.add_argument(
        '--bool_nofrags_check',
        action='store_true',
        help='advanced setting, turn off fragments check'
    )
    parser.add_argument(
        '--rs_confs_totnum',
        type=int,
        help=('advanced setting, exclusively for RotorSearch method, '
              'setup total number of conformers generation')
    )
    parser.add_argument(
        '--rs_confs_selnum',
        type=int,
        help=('advanced setting, exclusively for RotorSearch method, '
              'setup number of selections of generated conformers')
    )
    parser.add_argument(
        '--bool_confs_select_all',
        action='store_true',
        help='select all generated conformers, it will mute --rs_confs_selnum'
    )
    parser.add_argument(
        '--bool_calc_rmsd',
        action='store_true',
        help='turn on calculate RMSD'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return {}

    args = parser.parse_args(sys.argv[1:])
    if args.show_supported_types:
        print('Supported output file types:')
        print('; '.join(SUPPORTED_FORMATS))
        return {}
    if args.features:
        for i in FEATURES: print(i)
        return {}
    return vars(args)


def main():
    args = parsecmd()
    if not len(args): return
    SplitDatabase(**args)


if __name__ == '__main__':
    main()


