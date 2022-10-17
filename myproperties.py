#!/usr/bin/env python3

import os
import sys
import time
import collections
import tkinter as tk

# version 0.2.0     :   avoid filetype is type
# version 0.3.0     :   add fsize for each type of files
# version 0.4.0     :   correctly deal with the hidden file

USAGE = """
myproperties.py  [path]   :   show detail properties for file/dir
"""

if '-h' in sys.argv or '--help' in sys.argv:
    print(USAGE)
    sys.exit(0)


cwp = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~')
seplist = []


if not (os.path.isfile(cwp) or os.path.isdir(cwp)):
    print('Fatal: currently only file and dir are supproted')
    sys.exit(0)


def bytes2human(size):
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
        hsize = '{:.3f} Bytes'.format(size)
    return hsize


data = collections.OrderedDict()

data['@current_work_path'] = os.path.abspath(os.path.split(cwp)[0])
data['@name'] = os.path.basename(cwp)
if os.path.isdir(cwp):
    data['@type'] = 'dir'
elif os.path.isfile(cwp):
    data['@type'] = 'file'
else:
    data['@type'] = 'unknown'
data['==='] = ''       # work as separator
seplist.append('===')

data['created_time'] = time.ctime(os.path.getctime(cwp))
data['modified_time'] = time.ctime(os.path.getmtime(cwp))
data['=*='] = ''
seplist.append('=*=')


# specifically deal with file
if os.path.isfile(cwp):
    data['file_size'] = os.path.getsize(cwp)
else:
    data['view_file_num'] = 0
    data['view_file_total_size'] = 0.0
    data['view_dir_num'] = 0
    data['hidden_file_num'] = 0
    data['hidden_file_total_size'] = 0.0
    data['hidden_dir_num'] = 0
    data['unknown_file_num'] = 0
    data['unknown_file_total_size'] = 0.0

    adddict = {}
    for (dirpath, dirnames, filenames) in os.walk(cwp):
        dbase = os.path.basename(dirpath)       # care when ./.good
        if dbase.startswith('.'):
            data['hidden_dir_num'] += 1
            bo_hidden = True
        else:
            data['view_dir_num'] += 1
            bo_hidden = False
        for f in filenames:
            file = os.path.join(dirpath,f)
            if not os.path.isfile(file): continue       # os.walk does not follow symlink
            base,ext = os.path.splitext(f)
            fsize = os.path.getsize(file)
            if ext:
                pext = ext[1:].lower()
                if pext in adddict:
                    adddict[pext][0] += 1
                    adddict[pext][1] += fsize
                else:
                    adddict[pext] = [1,fsize]
            else:
                data['unknown_file_num'] += 1
                data['unknown_file_total_size'] += fsize
            if bo_hidden or base.startswith('.'):
                data['hidden_file_num'] += 1
                data['hidden_file_total_size'] += fsize
            else:
                data['view_file_num'] += 1
                data['view_file_total_size'] += fsize
    data['total_dir_num'] = data['view_dir_num'] + data['hidden_dir_num']
    data['total_file_num'] = data['view_file_num'] + data['hidden_file_num']
    data['total_file_size'] = data['view_file_total_size'] + data['hidden_file_total_size']

    for k,v in data.items():
        if 'size' in k:
            data[k] = bytes2human(v)

    data['=&='] = ''
    seplist.append('=&=')
    for k in sorted(adddict.keys()):
        s = bytes2human(adddict[k][1])
        data[k] = '{:} ({:})'.format(adddict[k][0],s)


g_width = 800
g_height = 600


root = tk.Tk()
root.geometry('{:}x{:}'.format(g_width,g_height))
root.title('XZProperties')

scrollbar = tk.Scrollbar(root)
scrollbar.pack(side='right', fill='y')

mylist = tk.Listbox(root, yscrollcommand=scrollbar.set, font=('courier',12))
for k,v in data.items():
    if k in seplist:
        info = '='*30
    else:
        info = '    {:30}: {:}'.format(k,v)
    mylist.insert(tk.END, info)
mylist.pack(side='left', fill='both', expand=True)

root.mainloop()


