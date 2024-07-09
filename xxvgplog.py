#!/usr/bin/env python3

import matplotlib.pyplot as plt

import os
import sys
import argparse


FEATURES = [
    'version 0.1.0  : plot for GMX-XVG file',
    'version 0.2.0  : powerful for continuously parse',
    'version 0.3.0  : support plumed data',
    'version 0.4.0  : make more versatile and calculate min/max/avg',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


def read_gmx_xvg(file):
    if not os.path.isfile(file):
        print('Fatal: not a file {:}'.format(file))
        return
    data = []
    legends = []
    bad = []
    fatal = []
    bo_plumed = False
    with open(file,'rt') as f:
        for line in f:
            l = line.strip()
            if not len(l): continue
            if l.startswith('#'):
                if l[:2] == '#!':   # plumed data
                    g = l[2:].split()
                    if g[0] == 'FIELDS':
                        bo_plumed = True
                        legends.extend(g[1:])
            elif l.startswith('@'):
                l = l[1:].strip()
                if l.startswith('legend'): continue
                if 'legend' in l:
                    g = l[1:].split()
                    if len(g) >= 3 and g[1] == 'legend':
                        k = ' '.join(g[2:]).replace('"',' ').replace("'",' ')
                        legends.append('-'.join(k.split()))
            else:
                try:
                    d = list(map(float,l.split()))
                except (ValueError,IndexError):
                    bad.append(l)
                else:
                    data.append(d)
    if not bo_plumed and legends:
        legends.insert(0,'time')
    # check length
    if data:
        n = len(legends) if legends else len(data[0])
        get = []
        for i,d in enumerate(data):
            if len(d) != n:
                fatal.append(' '.join([str(i) for i in d]))
                get.append(i)
        if len(get):
            for i in sorted(get,reverse=True):
                data.pop(i)
    if not data:
        return
    full = {
        'x'    : [d[0] for d in data],
        'bad'  : bad,
        'fatal': fatal,
        'size' : n,
    }
    keys = {}
    for i,e in enumerate(legends):
        full[e] = t = [d[i] for d in data]
        vmin = min(t)
        vmax = max(t)
        vavg = sum(t) / len(t)
        keys[str(i)] = [e,vmin,vmax,vavg]
    full['keys'] = keys
    return full


def plot_xvg(x,ys,titles=None):
    if len(ys) == 1:
        fig, ax = plt.subplots()
        ax.plot(x,ys[0])
        if titles:
            ax.set_title(titles[0])
        plt.show()
    else:
        if titles is None:
            titles = [None for i in len(ys)]
        fig, ax = plt.subplots()
        for y,t in zip(ys,titles):
            ax.plot(x,y,label=t)
        ax.legend()
        plt.show()


def main():
    parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         usage=f"""
plot_gmx_xvg {VERSION}
""",
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--file',
        dest='file',
        help='input file'
    )
    parser.add_argument(
        '-g','--group',
        dest='group',
        metavar='v',
        help=(
            'column indexes to plot groups, first column will be default to x-axis, '
            'thus, 1 means second column, 2 means third column. Names,case-insensitive, '
            'can also be used, abbreviation is allowed, thus `temperature` and `temp` or `te` '
            'can be the same column if unambiguously identified. Multiple groups '
            'are separated by ":", e.g., "1,2,3 : 1,temp : vir-xx, temp"'
        )
    )
    parser.add_argument(
        '-i','--interactive',
        action='store_true',
        help='interactive parse file'
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return
    w = parser.parse_args(sys.argv[1:])
    if w.features:
        for i in FEATURES: print(i)
        return

    full = read_gmx_xvg(w.file)
    if not full:
        print('Fatal: no valid data: file: {:}'.format(w.file))
        return

    bad = False
    if full['bad']:
        bad = True
        print('Lines in wrong format: not a number:')
        for l in full['bad']:
            print('  -> ',l)
    if full['fatal']:
        bad = True
        print('Lines in the same length:')
        for l in full['fatal']:
            print('  -> ',l)
    if bad:
        t = input('Want to continue? [y]/N: ')
        if not (not len(t.split()) or t.strip().lower() == 'y'):
            print('You have decided to quit, nothing will be ploted')
            return

    if w.interactive:
        print('Note: parsing keys: (index -> key  min  max  avg)')
        keys = full['keys']
        n = full['size']
        ask = None
        xaxis = 0
        xdata = full['x']
        outputs = []
        for i in range(n):
            l = '{:2}  ->  {:15}  {:15.3f}  {:15.3f}  {:15.3f}'.format(i, *keys[str(i)])
            outputs.append(l)
        while True:
            while True:
                print('Note: to reset xaxis, input: xaxis=i')
                if ask:
                    g = ask.split()
                    ask = None    # reset
                else:
                    for i in range(n):
                        if xaxis == i:
                            print(outputs[i] + '  (x-axis)')
                        else:
                            print(outputs[i])
                    t = input('Input group to be parsed (use number & space): ')
                    g = t.replace(',',' ').split()
                    for v in g:
                        if v.startswith('xaxis='):
                            try:
                                k = int(v.replace('xaxis=',' '))
                                if k < 0 or k > n: raise ValueError
                            except:
                                print('  -> wrong setting xaxis, no space is allowed')
                            else:
                                xaxis = k
                                xdata = full[keys[str(k)][0]]
                                continue
                if g:
                    try:
                        g = list(map(int,g))
                    except ValueError:
                        print('  -> wrong input: not a number: re-input')
                        continue
                    else:
                        if any([i<0 or i>n for i in g]):
                            print('  -> wrong input: not in range: re-input')
                            continue
                        select = [keys[str(v)][0] for v in g]
                        break
                else:
                    continue
            plot_xvg(x=xdata,ys=[full[i] for i in select],titles=select)

            ask = input('Continue? [y|numbers]/N : ')
            ask = ask.replace(',',' ')
            if not len(ask.split()) or ask.strip().lower() == 'y': continue
            bo = False
            for i in ask.split():
                try:
                    v = int(i)
                except ValueError:
                    bo = True
                    print('You have decided to quit, done')
                    break
            if bo:
                break
        return

    if w.group:
        use = list(full.keys())
        for k in ['bad','fatal','size','keys']: use.remove(k)
        select = []
        for g in w.group.split(':'):
            p = []
            for v in g.replace(',',' ').split():
                if v in use or v.isdigit():
                    if v.isdigit():
                        if v in full['keys'].keys():
                            p.append(full['keys'][v][0])
                        else:
                            print('Warning: group: unknown key <1>: {:}'.format(v))
                    else:
                        p.append(v)
                else:
                    found = ok = False
                    key = None
                    for k in use:
                        if k.startswith(v):
                            if found:
                                print('Warning: group: duplicate key: {:}'.format(v))
                                ok = False
                            else:
                                found = True
                                ok = True
                                key = k
                    if found:
                        if ok:
                            p.append(key)
                    else:
                        print('Warning: group: unknown key <2>: {:}'.format(v))
            if p:
                select.append(p)
    else:
        select = [[full['keys']['1'][0]], ]

    for g in select:
        plot_xvg(x=full['x'], ys=[full[i] for i in g], titles=g)


if __name__ == '__main__':
    main()



