#!/usr/bin/env python3

import os
import sys
import argparse

FEATURES = [
     'version 0.1.0      : May 18, 2022',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()

parser = argparse.ArgumentParser('Rename file by using given separator')
parser.add_argument('old',help='File/folder old')
parser.add_argument('new',nargs='?',help='File/folder new')
parser.add_argument('-s',help='separator',default='-')
args = parser.parse_args()

sep = '_'.join(args.s.split())
if args.new is None:
    new = sep.join(args.old.split())
else:
    new = args.new

if os.path.exists(new):
    print('Warning: already exits:',new)
    print('Continue? y/yes, else quit',end='  > ')
    if input().lower() not in ['y', 'yes']:
        sys.exit()

# double quote on old is very important
os.system('mv "{:}" {:}'.format(args.old,new))


