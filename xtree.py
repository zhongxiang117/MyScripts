#!/usr/bin/env python3

import os
import sys
import argparse
from colorama import Fore, Style


FEATURES = [
    'version 0.1.0  : Jan 19, 2024',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


class Tree:
    def __init__(
        self, cwd=None,
        show_max_n_files=None, show_all_files=None,
        sort_by_file_size=None, sort_by_file_mtime=None, sort_by_file_ctime=None, sort_by_file_name=None,
        show_max_n_dirs=None, show_all_dirs=None,
        sort_by_dir_size=None, sort_by_dir_mtime=None, sort_by_dir_ctime=None, sort_by_dir_name=None,
        reverse_sort_dir=None, reverse_sort_file=None,
        show_size=None,
        *args,**kws,
    ):
        self.cwd = cwd
        self.show_max_n_files = show_max_n_files if show_max_n_files else 3
        self.show_max_n_dirs = show_max_n_dirs if show_max_n_dirs else 3
        self.marks = ['├──', '│', '└──']
        self.sort_by_file_size = sort_by_file_size
        self.sort_by_file_ctime = sort_by_file_ctime
        self.sort_by_file_mtime = sort_by_file_mtime
        self.sort_by_file_name = sort_by_file_name
        self.sort_by_dir_size = sort_by_dir_size
        self.sort_by_dir_ctime = sort_by_dir_ctime
        self.sort_by_dir_mtime = sort_by_dir_mtime
        self.sort_by_dir_name = sort_by_dir_name
        self.reverse_sort_dir = reverse_sort_dir if reverse_sort_dir is True else False
        self.reverse_sort_file = reverse_sort_file if reverse_sort_file is True else False
        self.show_size = show_size
        self.show_all_files = show_all_files
        self.show_all_dirs = show_all_dirs

    def run(self,cwd=None):
        if not cwd: cwd = self.cwd
        if not cwd: cwd = '.'
        if not os.path.isdir(cwd): return
        now = os.getcwd()
        os.chdir(cwd)
        fdict = self.get_fdict()
        self.accumulate_size(fdict)
        print()
        self.xprint(fdict)
        os.chdir(now)

    def xprint(self,fdict,tab=None,nested=None):
        if not tab: tab = ''
        beg = self.marks[0]
        mid = self.marks[1]
        end = self.marks[2]
        arch = fdict.pop('/\\arch')
        dirs = [k for k,v in fdict.items() if isinstance(v,dict)]
        files = [k for k in fdict.keys() if k not in dirs]
        if files:
            if self.sort_by_file_name:
                files = sorted(files, reverse=self.reverse_sort_file)
            elif self.sort_by_file_size:
                files = sorted(files, key=lambda x: fdict[x][0], reverse=self.reverse_sort_file)
            elif self.sort_by_file_mtime:
                files = sorted(files, key=lambda x: fdict[x][1], reverse=self.reverse_sort_file)
            elif self.sort_by_file_ctime:
                files = sorted(files, key=lambda x: fdict[x][2], reverse=self.reverse_sort_file)

            n = 0
            if not self.show_all_files:
                if len(files) >= self.show_max_n_files:
                    n = len(files) - self.show_max_n_files
                    files = [files[i] for i in range(self.show_max_n_files)]

            if nested:
                if len(files) == 1:
                    print(self._fout(tab+beg, files[0], fdict[files[0]][0]))
                elif len(files) == 2:
                    print(self._fout(tab+beg, files[0], fdict[files[0]][0]))
                    print(self._fout(tab+end, files[1], fdict[files[1]][0]))
                else:
                    for i in range(len(files)-1):
                        print(self._fout(tab+beg, files[i], fdict[files[i]][0]))
                    print(self._fout(tab+end, files[-1], fdict[files[-1]][0]))
                if n:
                    print(tab+'    ... ({:} files)'.format(n))
            else:
                last = files.pop(-1)
                for f in files:
                    print(self._fout(tab+beg, f, fdict[f][0]))
                if dirs:
                    print(self._fout(tab+beg, last, fdict[last][0]))
                else:
                    print(self._fout(tab+end, last, fdict[last][0]))
                if n:
                    print(tab+mid+'  ... ({:} files)'.format(n))

        n = 0
        if dirs:
            if self.sort_by_dir_name:
                dirs = sorted(dirs, reverse=self.reverse_sort_dir)
            elif self.sort_by_dir_size:
                dirs = sorted(dirs, key=lambda x: fdict[x]['/\\arch'][0], reverse=self.reverse_sort_dir)
            elif self.sort_by_dir_mtime:
                dirs = sorted(dirs, key=lambda x: fdict[x]['/\\arch'][1], reverse=self.reverse_sort_dir)
            elif self.sort_by_dir_ctime:
                dirs = sorted(dirs, key=lambda x: fdict[x]['/\\arch'][2], reverse=self.reverse_sort_dir)

            if not self.show_all_dirs:
                if len(dirs) > self.show_max_n_dirs:
                    n = len(dirs) - self.show_max_n_dirs
                    dirs = [dirs[i] for i in range(self.show_max_n_dirs)]
            last = dirs.pop(-1)
            for d in dirs:
                print(self._fout(tab+Fore.BLUE+d+'  '+Style.RESET_ALL, '', fdict[d]['/\\arch'][0]))
                self.xprint(fdict[d],tab=tab+mid+'  ',nested=True)
            print(self._fout(tab+Fore.BLUE+last+'  '+Style.RESET_ALL, '', fdict[last]['/\\arch'][0]))
            self.xprint(fdict[last],tab=tab+'  ',nested=True)
        if nested:
            if n:
                print(mid+'  ... ({:} dirs)'.format(n))
        else:
            if n:
                print('... ({:} dirs)'.format(n))
            print('(total: {:})\n'.format(self.bytes2human(arch[0])))

    def _fout(self,prefix,file,size):
        if self.show_size:
            size = self.bytes2human(size)
            sz = '({:}) {:}'.format(size,file)
            return prefix+sz
        return prefix+file

    def bytes2human(self,size):
        unit = 1024
        # order is important
        if size > 1000*unit*unit:   # 1000 MB
            size = size / unit / unit / unit
            hsize = '{:.3f} GB'.format(size)
        elif size > 1000*unit:      # 1000 KB
            size = size / unit / unit
            hsize = '{:.3f} MB'.format(size)
        elif size > 1000:
            size = size / unit
            hsize = '{:.3f} KB'.format(size)
        else:
            hsize = '{:.3f} B'.format(size)
        return hsize

    def get_fdict(self,dir=None):
        if not dir: dir = '.'
        # fdict = {
        #     'file1.csv': (10293731, 1697176420.0, 1704683570.618744),
        #     'file2.py': (38307, 1704966868.7494898, 1704966868.7494898),
        #     '.vscode': {
        #         'launch.json': (2302, 1705645470.2488985, 1705645470.2488985),
        #         'settings.json': (129, 1704683986.352312, 1704683986.352312)
        #     },
        #     'dir1': {
        #         'dirA': {
        #             'fileA.pdb': (663410, 1670244072.0, 1704683991.1563532),
        #             'fileB.sdf': (2917, 1670244072.0, 1704683991.1523533),
        #           },
        #         'dirB': {
        #             'fileX.mol2': (4646, 1670244168.0, 1704683991.492356),
        #         }
        #     }
        # }
        fdict = {}
        for (dirpath, dirnames, filenames) in os.walk(dir):     # `dirpath`: ./a/b/c/d
            g = fdict       # alias
            for k in dirpath.split(os.path.sep):
                g = g.setdefault(k,{})
            total = 0.0
            for f in filenames:
                file = os.path.join(dirpath,f)
                p = os.stat(file)
                g[f] = (p.st_size,p.st_mtime,p.st_ctime)
                total += p.st_size
            d = os.stat(dirpath)
            g['/\\arch'] = (total,d.st_mtime,d.st_ctime)
        fdict = fdict[dir]
        return fdict

    def accumulate_size(self,fdict):
        get = fdict['/\\arch']
        now = get[0]
        for v in fdict.values():
            if isinstance(v,dict):
                now += self.accumulate_size(v)
        fdict['/\\arch'] = (now,get[1],get[2])
        return now


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=f'xtree {VERSION}',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        'cwd',
        nargs='?',
        help='Current work directory'
    )

    gf = parser.add_argument_group('options for file')
    gf.add_argument(
        '-sf', '--show-all-files',
        action='store_true',
        help='Force to show all files, highest priority',
    )
    gf.add_argument(
        '-nf', '--show-max-n-files',
        type=int,
        metavar='v',
        help='Show maxmimum this number of files (default:3)',
    )
    gf.add_argument(
        '-fs',  '--sort-by-file-size',
        action='store_true',
        help='Sort by file size',
    )
    gf.add_argument(
        '-fn',  '--sort-by-file-name',
        action='store_true',
        help='Sort by file name',
    )
    gf.add_argument(
        '-fm', '--sort-by-file-mtime',
        action='store_true',
        help='Sort by file modification time',
    )
    gf.add_argument(
        '-fc', '--sort-by-file-ctime',
        action='store_true',
        help='Sort by file creation time',
    )

    gd = parser.add_argument_group('options for dir')
    gd.add_argument(
        '-sd', '--show-all-dirs',
        action='store_true',
        help='Force to show all dirs, highest priority',
    )
    gd.add_argument(
        '-nd', '--show-max-n-dirs',
        type=int,
        metavar='v',
        help='Show maxmimum this number of dirs (default:3)',
    )
    gd.add_argument(
        '-ds',  '--sort-by-dir-size',
        action='store_true',
        help='Sort by dir size',
    )
    gd.add_argument(
        '-dn',  '--sort-by-dir-name',
        action='store_true',
        help='Sort by dir name',
    )
    gd.add_argument(
        '-dm', '--sort-by-dir-mtime',
        action='store_true',
        help='Sort by file modification time',
    )
    gd.add_argument(
        '-dc', '--sort-by-dir-ctime',
        action='store_true',
        help='Sort by dir creation time',
    )

    parser.add_argument(
        '-ss', '--show-size',
        action='store_true',
        help='show size'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )

    w = parser.parse_args()
    if w.features:
        for i in FEATURES: print(i)
        return

    T = Tree(**vars(w))
    T.run()


if __name__ == '__main__':
    main()



