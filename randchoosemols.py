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
import random
import argparse

FEATURES = [
    'version 0.1.0  : Randomly choose molecules based on given database',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
_fms = openbabel.OBConversion().GetSupportedInputFormat()
SUPPORTED_FORMATS = [i.split()[0].lower() for i in _fms]


class ReadFileOBMols:
    """read data from input file

    Attributes:
        filename(list): input file name
        fileobmols(List[openbabel.OBMol, ]):
    """
    def __init__(self,filename=None,*args,**kwargs):
        self.filename = filename

    def run(self,filename=None):
        curtime = time.time()
        if filename is None: filename = self.filename
        self.fileobmols = self.readfile(filename)
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
            return []

        obconv = openbabel.OBConversion()
        obconv.SetInFormat(ext)

        obmollist = []
        obmol = openbabel.OBMol()
        bo = obconv.ReadFile(obmol,file)
        while bo:
            obmollist.append(obmol)
            # continue reading
            obmol = openbabel.OBMol()
            bo = obconv.Read(obmol)
        return obmollist


class RandChooseMols(ReadFileOBMols):
    def __init__(self,n=None,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.n = n

    def run(self):
        super().run()
        self.obmols = []
        if not self.fileobmols: return
        if self.n > len(self.fileobmols):
            print(f'Warning: too many choices: ({self.n} > {len(self.fileobmols)})')
            return
        self.obmols = random.sample(self.fileobmols,k=self.n)


class SaveFile:
    """save final results"""
    def __init__(self,obmols=None,ftype=None,fname=None,*args,**kwargs):
        if ftype:
            self.ftype = ftype.lower()
            if self.ftype not in SUPPORTED_FORMATS:
                print('Fatal: currently file format not support: ',self.ftype)
                self.ftype = 'sdf'
        else:
            self.ftype = 'sdf'
        if isinstance(obmols,list):
            self.obmols = obmols
        else:
            print('Fatal: not data input')
            self.obmols = None
        self.fname = fname if fname else 'artchoices'
        if self.obmols: self.saveobmols()

    def saveobmols(self):
        # check existing files to avoid overwriting
        fnew = self._get_fnew(self.fname, self.ftype)
        fnew += '.' + self.ftype
        print('Note: result file is saved to: '+fnew)
        obconv = openbabel.OBConversion()
        obconv.SetOutFormat(self.ftype)
        obconv.WriteFile(self.obmols[0],fnew)
        for obmol in self.obmols[1:]:
            obconv.Write(obmol)
        obconv.CloseOutFile()

    def _get_fnew(self,fname,ftype):
        i = 0
        while True:
            i += 1
            if not os.path.isfile(fname+'-'+str(i)+'.'+ftype):
                break
        return fname+'-'+str(i)


def parsecmd():
    parser = argparse.ArgumentParser(
        description='Random Choose Number of Molecules on Given Database',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--filename',
        help='input filename, database'
    )
    parser.add_argument(
        '-n',
        type=int,
        help='number of choose'
    )
    parser.add_argument(
        '--fname',
        help='output file basename'
    )
    parser.add_argument(
        '--ftype',
        help='output file type, check by option: --show_supported_ftypes',
    )
    parser.add_argument(
        '--show_supported_ftypes',
        action='store_true',
        help='show supported file output types',
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
    t = RandChooseMols(**args)
    t.run()
    args['obmols'] = t.obmols
    SaveFile(**args)


if __name__ == '__main__':
    main()


