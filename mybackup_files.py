#!/usr/bin/env python3

import os
import sys
import json
import time
import shutil

# version 0.1.0     : Nov 2nd, 2022
# version 0.2.0     : deal with repeats
# version 0.3.0     : more powerful
#       Logical:
#           recursively find files starting with `#identifier` under `start`,
#           then based on the parameters stored in `par.json`, either got from
#           `start` path or most-recently time-stamped folder, calculate lists
#           `replist` (files need to be overwritten) and `newfiles` (files need
#           to be firstly backed up), then create time-stamped folder for old
#           files in `replist` and update the new ones, and copy files in `newfiles`.
#           Specially, files defined in `parfile` will always exist under `cwd`;
#           folders defined in `rmlist` will be excluded.

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

print(f'Note: start: backing up folder: {start}')
rstlist = []
for (dirpath, dirnames, filenames) in os.walk(start):
    # cleanup dirnames, in-place
    ignores = [i for i in dirnames if i.startswith('.') or i.startswith('master-')]
    for i in ignores: dirnames.remove(i)
    for k in rmlist:
        if k in dirnames:
            dirnames.remove(k)
    for f in filenames:
        #TODO identifier
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
parold = parfile
if not os.path.isfile(parold):
    for i in range(len(curdirs)-1,-1,-1):
        file = os.path.join(curdirs[i],parfile)
        if os.path.isfile(file):
            parold = file
            break
if os.path.isfile(parold):
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
            if 'fbak' in objold[k]:
                repfiles.append(objold[k]['fbak'])
            else:
                repfiles.append(os.path.basename(k))
            newfiles.append(k)
    else:
        newfiles.append(k)


if repfiles:
    dirname = time.strftime('%Y-%m-%d-%H-%M-%S')
    os.makedirs(dirname)
    if os.path.isfile(parfile):
        print(f'Note: moving old parfile: {parfile} -> {dirname}')
        shutil.move(parfile,dirname)
    for k in repfiles:
        if os.path.isfile(k):
            print(f'Note: moving old file: {k} -> {dirname}')
            shutil.move(k,dirname)
if newfiles:
    cwd = os.getcwd()
    for k in newfiles:
        base = os.path.basename(k)
        if os.path.isfile(base):
            i = 1
            while True:
                newbase = base + '-' + str(i)
                if not os.path.isfile(newbase):
                    break
                i += 1
            objnew[k]['fbak'] = newbase
            print(f'Note: backing up file: {k} ++ {i}')
            shutil.copy(k,os.path.join(cwd,newbase))
        else:
            print(f'Note: backing up file: {k}')
            shutil.copy(k,cwd)

    # double check
    for k in objnew:
        cfp = None
        cfn = None
        if 'fbak' in objnew[k]:
            if not os.path.isfile(objnew[k]['fbak']):
                cfp = k
                cfn = objnew[k]['fbak']
        else:
            base = os.path.basename(k)
            if not os.path.isfile(base):
                cfp = k
                cfn = base
        if cfp:
            print(f'Warning: backup missing file: {cfp}  ++ {cfn}')
            shutil.copy(cfp,os.path.join(cwd,cfn))

    with open(parfile,'w') as f:
        json.dump(objnew,f,indent=4)




