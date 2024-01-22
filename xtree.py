#!/usr/bin/env python3

import os
import argparse
from colorama import Fore, Style


FEATURES = [
    'version 0.1.0  : Jan 19, 2024',
    'version 0.2.0  : add option `--not-include-folder-size`',
    'version 0.3.0  : add info of number of dirs and files',
    'version 0.4.0  : add more useful options',
    'version 0.5.0  : indicate whether file executable',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


class Tree:
    def __init__(
        self, cwd=None, not_include_folder_size=None,
        show_max_n_files=None, show_all_files=None,
        sort_by_file_size=None, sort_by_file_mtime=None, sort_by_file_ctime=None, sort_by_file_name=None,
        show_max_n_dirs=None, show_all_dirs=None,
        sort_by_dir_size=None, sort_by_dir_mtime=None, sort_by_dir_ctime=None, sort_by_dir_name=None,
        reverse_sort_dir=None, reverse_sort_file=None,
        sort_dir_by_dirnum=None, sort_dir_by_filenum=None,
        show_size=None,
        *args,**kws,
    ):
        self.cwd = cwd
        self.not_include_folder_size = not_include_folder_size
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
        self.sort_dir_by_dirnum = sort_dir_by_dirnum
        self.sort_dir_by_filenum = sort_dir_by_filenum

    def run(self,cwd=None):
        if not cwd: cwd = self.cwd
        if not cwd: cwd = '.'
        if not os.path.isdir(cwd): return
        now = os.getcwd()
        os.chdir(cwd)
        fdict = self.get_fdict(not_include_folder_size=self.not_include_folder_size)
        self.accumulate_size(fdict)
        print('\n'+Fore.BLUE+now+Style.RESET_ALL)
        self.xprint(fdict)
        os.chdir(now)

    def xprint(self,fdict,tab=None,nested=None):
        if not tab: tab = ''
        beg = self.marks[0]
        mid = self.marks[1]
        end = self.marks[2]
        arch = fdict.pop('/\\arch')         # be aware, additional key
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
                if dirs:
                    for i in range(len(files)):
                        print(self._fout(tab+beg, files[i], fdict[files[i]]))
                else:
                    if len(files) == 1:
                        print(self._fout(tab+end, files[0], fdict[files[0]]))
                    elif len(files) == 2:
                        print(self._fout(tab+beg, files[0], fdict[files[0]]))
                        print(self._fout(tab+end, files[1], fdict[files[1]]))
                    else:
                        for i in range(len(files)-1):
                            print(self._fout(tab+beg, files[i], fdict[files[i]]))
                        print(self._fout(tab+end, files[-1], fdict[files[-1]]))
                    if n:
                        print(tab+'   ... ({:} files)'.format(n))
            else:
                last = files.pop(-1)
                for f in files:
                    print(self._fout(tab+beg, f, fdict[f]))
                if dirs:
                    print(self._fout(tab+beg, last, fdict[last]))
                else:
                    print(self._fout(tab+end, last, fdict[last]))
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
            elif self.sort_dir_by_dirnum:
                dirs = sorted(dirs, key=lambda x: fdict[x]['/\\arch'][3], reverse=self.reverse_sort_dir)
            elif self.sort_dir_by_filenum:
                dirs = sorted(dirs, key=lambda x: fdict[x]['/\\arch'][4], reverse=self.reverse_sort_dir)

            if not self.show_all_dirs:
                if len(dirs) > self.show_max_n_dirs:
                    n = len(dirs) - self.show_max_n_dirs
                    dirs = [dirs[i] for i in range(self.show_max_n_dirs)]

            last = dirs.pop(-1)
            for d in dirs:
                print(self._dout(tab+beg+Fore.BLUE+d+Style.RESET_ALL, '', fdict[d]['/\\arch'][0]))
                self.xprint(fdict[d],tab=tab+mid+'  ',nested=True)

            print(self._dout(tab+end+Fore.BLUE+last+Style.RESET_ALL, '', arch[0]))
            self.xprint(fdict[last],tab=tab+'   ',nested=True)      # no `mid`, additional space

            if n:
                print(tab+'   ... ({:} dirs)'.format(n))

        if not nested:
            print(
                '(total: {:}, dirnum: {:}, filenum: {:})\n'.format(
                    self.bytes2human(arch[0]), arch[3], arch[4]
                )
            )

    def _fout(self,prefix,file,stat):
        if self.show_size:
            size = self.bytes2human(stat[0])
            new = '({:}) {:}'.format(size,file)
        else:
            new = file
        if stat[3]:
            return prefix+Fore.GREEN+new+Style.RESET_ALL
        return prefix+new

    def _dout(self,prefix,file,size):
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

    def get_fdict(self,dir=None,not_include_folder_size=None):
        if not dir: dir = '.'
        ## `VALUES` is a tuple
        # fdict = {
        #     'file1.csv': VALUES,
        #     'file2.py': VALUES,
        #     '.vscode': {
        #         'launch.json': VALUES,
        #         'settings.json': VALUES,
        #     },
        #     'dir1': {
        #         'dirA': {
        #             'fileA.pdb': VALUES,
        #             'fileB.sdf': VALUES,
        #           },
        #         'dirB': {
        #             'fileX.mol2': VALUES,
        #         }
        #     }
        # }
        fdict = {}
        for (dirpath, dirnames, filenames) in os.walk(dir):     # `dirpath`: ./a/b/c/d
            g = fdict       # alias
            for k in dirpath.split(os.path.sep):
                g = g.setdefault(k,{})
            total = 0.0 if not_include_folder_size else 4*1024.0
            for f in filenames:
                file = os.path.join(dirpath,f)
                p = os.stat(file)
                g[f] = (p.st_size,p.st_mtime,p.st_ctime,os.access(file,os.X_OK))
                total += p.st_size
            d = os.stat(dirpath)
            g['/\\arch'] = (total,d.st_mtime,d.st_ctime,len(dirnames),len(filenames))
        fdict = fdict[dir]
        return fdict

    def accumulate_size(self,fdict):
        get = fdict['/\\arch']
        now = [get[0],get[3],get[4]]
        for v in fdict.values():
            if isinstance(v,dict):
                new = self.accumulate_size(v)
                now[0] += new[0]
                now[1] += new[1]
                now[2] += new[2]
        fdict['/\\arch'] = (now[0],get[1],get[2],now[1],now[2])
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
    parser.add_argument(
        '-ss', '--show-size',
        action='store_true',
        help='show size'
    )
    parser.add_argument(
        '-I', '--not-include-folder-size',
        action='store_true',
        help='an empty folder takes up 4096 Bytes (4KB) space, this option will exclude it'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
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
    gf.add_argument(
        '-rf', '--reverse-sort-file',
        action='store_true',
        help='Reverse result of file sort',
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
    gd.add_argument(
        '-dd', '--sort-dir-by-dirnum',
        action='store_true',
        help='Sort by dir creation time',
    )
    gd.add_argument(
        '-df', '--sort-dir-by-filenum',
        action='store_true',
        help='Sort by dir creation time',
    )
    gd.add_argument(
        '-rd', '--reverse-sort-dir',
        action='store_true',
        help='Reverse result of dir sort',
    )

    w = parser.parse_args()
    if w.features:
        for i in FEATURES: print(i)
        return

    T = Tree(**vars(w))
    T.run()


if __name__ == '__main__':
    main()



