#!/usr/bin/env python3

try:
    # version 3
    from openbabel import openbabel
except ImportError:
    # version 2
    import openbabel

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
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
_ffs = openbabel.vectorString()
openbabel.OBPlugin.ListAsVector('forcefields',None,_ffs)
SUPPORTED_FFS = [i.split()[0].lower() for i in _ffs]
OBConv = openbabel.OBConversion()
_fms = OBConv.GetSupportedInputFormat()
SUPPORTED_FORMATS = [i.split()[0].lower() for i in _fms]


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
            data.SetValue(str(ene))
            obtmp.CloneData(data)
            setff.UpdateCoordinates(obtmp)
            obmollist.append(obtmp)

        if not len(obmollist):
            # means done on initialization or before first N steps
            setff.UpdateCoordinates(obmol)
            return [obmol, ]

        if len(obmollist) <= self.rs_confs_selnum:
            return obmollist

        # sort obmollist from energy smallest to largest
        ndxlist = sorted(range(len(enelist)),key=lambda k: enelist[k])
        obmollist = [obmollist[i] for i in ndxlist]

        if self.rs_confs_selnum == 1:
            return [obmollist[0], ]
        
        if self.rs_confs_selnum == 2:
            return [obmollist[0], obmollist[-1]]
        
        dt = (len(obmollist) - 1) // (self.rs_confs_selnum - 1)
        reflist = list(range(0,len(obmollist),dt))
        if len(reflist) > self.rs_confs_selnum:
            reflist = reflist[:self.rs_confs_selnum]
        return [obmollist[i] for i in reflist]

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
        '--show_supported_ftypes',
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
        '--features',
        action='store_true',
        help='show development features'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return {}

    args = parser.parse_args(sys.argv[1:])
    if args.show_supported_ftypes:
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


