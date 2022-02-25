#!/usr/bin/env python3

import os
import sys
import argparse

FEATURES = [
    'version 0.1.0     : Feb 23, 2022',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()

parser = argparse.ArgumentParser(
    description='resourcesum.py',
    allow_abbrev=False,
)
parser.add_argument(
    '-v','--version',
    action='version',
    version=VERSION
)
parser.add_argument(
    'pattern',
    nargs='?',
    help='pattern to be searched'
)
parser.add_argument(
    '-a','--all',
    action='store_true',
    help='bool, show all, default False'
)
parser.add_argument(
    '-u','--user',
    help='specify which user to lookup'
)
parser.add_argument(
    '--features',
    action='store_true',
    help='show development features'
)

if len(sys.argv) <= 1:
    parser.print_help()
    exit()

args = parser.parse_args(sys.argv[1:])
if args.features:
    for i in FEATURES: print(i)
    exit()

# note: info is a str type, delimited by new line '\n'
info = os.popen('ps aux | grep -i {:}'.format(args.pattern)).read()

usage = []
# format: USER  PID  %CPU  %MEM  TIME  COMMAND
for line in info.split('\n'):
    its = line.split()
    if len(its) < 10: continue
    usage.append([its[0],its[1],its[2],its[3],its[9],' '.join(its[10:])])

if args.user:
    usage = [i for i in usage if i[0].lower() == args.user.lower()]

totcpu = sum([float(i[2]) for i in usage])
totmem = sum([float(i[3]) for i in usage])
times = [i[4].split(':') for i in usage]
tottime = sum([float(i[0])*60+float(i[1]) for i in times])
if args.all:
    print('USER      PID      %CPU   %MEM   TIME    COMMAND')
    for i in usage:
        print('{:6}  {:10} {:5}  {:5}  {:6}  {:}'.format(*i))
print('Sum: totcpu:{:.2f}%  totmem:{:.2f}%  tottime:{:}s'.format(totcpu,totmem,tottime))


