#!/usr/bin/env python3

import os
import sys
import json
import time
import shutil

# version 0.1.0     :    Nov 2nd, 2022

USAGE = """
    [python3] backup.py         :  backup `parse-*.drawio' files
"""

if '-h' in sys.argv or '--help' in sys.argv:
    print(USAGE)
    sys.exit(0)

# path to recursively backuping
start = os.path.expanduser('~')

# name of parameter file (do not change)
parfile = 'par.json'

# exclude folders
rmlist = [
    'Applications',
    'bin',
    'lib',
    'opt',
    'snap',
    'intel',
    'System Volume Information',
    'Temp',
    'Templates',
    'VirtualBox',
    'zoteroFiles',
    'Zotero',
    'zhongxiang117.github.io',

    '__pycache__',
    '__MACOSX',
    'node_modules',
    'csds_data',
    'Discovery Studio',
    'ArtGUI',
    'ArtGUI-versions',
    'DONEBackup',
    'MEGA X',
    'nvvp_workspace',
]


rstlist = []
for (dirpath, dirnames, filenames) in os.walk(start):
    # cleanup dirnames, in-place
    ignores = [i for i in dirnames if i.startswith('.') or i.startswith('master-')]
    for i in ignores: dirnames.remove(i)
    for k in rmlist:
        if k in dirnames:
            dirnames.remove(k)
    for f in filenames:
        if f.startswith('parse-') or f.startswith('parser-'):
            rstlist.append(os.path.join(dirpath,f))
objnew = {}
for file in rstlist:
    s = os.stat(file)
    objnew.update({
        file : {
            'size'  : s.st_size,
            'mtime' : s.st_mtime
        }
    })


curdirs = [d for d in sorted(os.listdir()) if os.path.isdir(d)]
objold = {}
if os.path.isfile(parfile):
    parold = parfile
else:
    for i in range(len(curdirs)-1,-1,-1):
        file = os.path.join(curdirs[i],parfile)
        if os.path.isfile(file):
            parold = file
            break
if os.path.isfile(parfile):
    with open(parold,'rt') as f:
        o = json.load(f)
        for k,v in o.items():
            objold.update({k:v})


newfiles = []
repfiles = []
for k in objnew:
    if k in objold:
        if objold[k]['size'] == objnew[k]['size'] and objold[k]['mtime'] == objnew[k]['mtime']:
            pass
        else:
            repfiles.append(k)
            newfiles.append(k)
    else:
        newfiles.append(k)


# backup first
if repfiles:
    dirname = time.strftime('%Y-%m-%d-%H-%M-%S')
    os.makedirs(dirname)
    if os.path.isfile(parfile):
        print(f'Note: moving old parfile: {parfile} -> {dirname}')
        shutil.move(parfile,dirname)
    for k in repfiles:
        base = os.path.basename(k)
        print(f'Note: moving old file: {base} -> {dirname}')
        shutil.move(base,dirname)
if newfiles:
    cwd = os.getcwd()
    for k in newfiles:
        print(f'Note: backing up file: {k}')
        shutil.copy(k,cwd)


with open(parfile,'w') as f:
    json.dump(objnew,f,indent=4)



