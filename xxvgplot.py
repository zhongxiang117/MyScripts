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
    'version 0.5.0  : deal with when `xaxis` is used',
    'version 0.6.0  : fix potential issue for plumed file of multiple FIELDS',
    'version 0.7.0  : refactor, support GMX log file',
    'version 0.8.0  : add `begin`',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


class ReadFile:
    """
    Attributes:
        full (dict):
            'x'    : List[float],   for x-axis label
            'raw'  : List[str],     wrong parsed line
            'bad'  : List[float],   inconsistent data
            'size' : int,           number of kinds data parsed
            'keys' : List[str],     keys index
    """
    GMX_PROPERTIES = [
        'Bond', 'Angle', 'Proper Dih.', 'Per. Imp. Dih.', 'LJ-14', 'Coulomb-14',
        'LJ (SR)', 'Disper. corr.', 'Coulomb (SR)', 'Coul. recip.', 'Potential',
        'Kinetic En.', 'Total Energy', 'Conserved En.', 'Temperature',
        'Pres. DC (bar)', 'Pressure (bar)', 'Constr. rmsd'
    ]
    def __init__(self,file,*args,**kws):
        self.nice = False
        ext = ''
        self.full = []
        if os.path.isfile(file):
            base,ext = os.path.splitext(file)
            self.nice = True
            with open(file,'rt') as f:
                prolines = f.readlines()
        elif isinstance(file,list):
            if all([isinstance(t,str) for t in file]):
                self.nice = True
                prolines = file
            elif all([all([isinstance(i,(float,int)) for i in t]) for t in file]):
                #TODO
                self.full = []
        if self.nice:
            if ext in ['.xvg', '']:
                self.read_xvg(prolines)
                if not self.full:
                    self.read_gmx_log(prolines)
            else:
                self.read_gmx_log(prolines)

    def read_xvg(self,prolines):
        """
        Example:

        # GMX comments
        @    title "RMSD"
        @    xaxis  label "Time (ps)"
        @    yaxis  label "RMSD (nm)"
        @TYPE xy
        @ subtitle "C-alpha after lsq fit to Protein"
        0.0000000     0.1231007
        10.0000000    0.1196074

        # Plumed comments
        #! FIELDS time phi psi metad.bias
        #! SET min_phi -pi
        #! SET max_phi pi
        #! SET min_psi -pi
        #! SET max_psi pi
        0.000000 -2.879362 2.919290 0.000000
        0.020000 -2.788458 2.883132 0.000000
        """
        if not prolines: return
        data = []
        legends = []
        axes = []
        title = []
        raw = []
        bo_plumed = False
        for line in prolines:
            l = line.strip()
            if not len(l): continue
            if l.startswith('#'):
                if l[:2] == '#!':   # plumed data
                    g = l[2:].split()
                    if g[0] == 'FIELDS':
                        if bo_plumed:
                            print('Warning: multiple FIELDS')
                            print('  -> ',l)
                        bo_plumed = True
                        legends = g[1:]
            elif l.startswith('@'):
                l = l[1:].lstrip()
                if 'legend' in l:
                    if l.startswith('legend'): continue
                    g = l.replace('"',' ').replace("'",' ').split()
                    if len(g) >= 3 and g[1] == 'legend':
                        legends.append(' '.join(g[2:]))
                elif l.startswith('xaxis') or l.startswith('yaxis'):
                    g = l.replace('"',' ').replace("'",' ').split()
                    if len(g) >= 3:
                        axes.append(' '.join(g[2:]))
                elif l.startswith('title') or l.startswith('subtitle'):
                    g = l.replace('"',' ').replace("'",' ').split()
                    if len(g) >= 2:
                        title.append(' '.join(g[1:]))
            else:
                try:
                    d = list(map(float,l.split()))
                except (ValueError,IndexError):
                    raw.append(l)
                else:
                    data.append(d)
        if axes:
            legends = axes
        elif not bo_plumed and legends:
            legends.insert(0,'time')
        if title:
            title = ' / '.join(title)
        else:
            title = 'Plumed Data' if bo_plumed else 'GMX Data'
        
        # check length
        bad = []
        if data:
            n = len(legends) if legends else len(data[0])
            get = []
            for i,d in enumerate(data):
                if len(d) != n:
                    bad.append(d)
                    get.append(i)
            if len(get):
                for i in sorted(get,reverse=True):
                    data.pop(i)
        if not data:
            return
        full = {
            'x'    : [d[0] for d in data],
            'bad'  : bad,
            'raw'  : raw,
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
        self.full = full

    def read_gmx_log(self,prolines):
        if not prolines: return
        i = 0
        full = {'x':[], 'bad':[], 'raw':[]}
        n = len(prolines) - 1
        while i < n:
            si = '-'.join(prolines[i].split())
            if si == 'Step-Time':
                try:
                    t = float(prolines[i+1].split()[1])
                except:
                    t = 0.0
                full['x'].append(t)
                j = i + 4
                while j < n:
                    l = prolines[j] + ' '*80
                    sj = '-'.join(l.split())
                    if sj == 'Step-Time':
                        break
                    t = [l[b:b+15].strip() for b in range(0,len(l)-15,15)]
                    t = [v for v in t if v]
                    if not any([v in self.GMX_PROPERTIES for v in t]):
                        j += 1
                        continue
                    j += 1
                    if j > n: break
                    ld = prolines[j] + ' '*80
                    d = [ld[b:b+15].strip() for b in range(0,len(ld)-15,15)]
                    d = [v for v in d if v]
                    if len(d) == len(t):
                        for k,v in zip(t,d):
                            try:
                                v = float(v)
                            except:
                                v = 0.0
                            full.setdefault(k,[]).append(v)
                    else:
                        # full['raw'].append(l.strip())
                        # full['raw'].append(ld.strip())
                        pass
                    j += 1      # next entry
                i = j
            else:
                i += 1
        skeys = list(full.keys())
        for k in ['x','bad','raw']: skeys.remove(k)
        keys = {}
        for i,e in enumerate(['x',*skeys]):
            t = full[e]
            vmin = min(t)
            vmax = max(t)
            vavg = sum(t) / len(t)
            keys[str(i)] = [e,vmin,vmax,vavg]
        full['keys'] = keys
        full['size'] = len(keys)
        self.full = full


def plot_xvg(x,ys,titles=None):
    if len(ys) == 1:
        n = min(len(x),len(ys[0]))
        fig, ax = plt.subplots()
        ax.plot(x[:n],ys[0][:n])
        if titles:
            ax.set_title(titles[0])
        plt.show()
    else:
        if titles is None:
            titles = [None for i in len(ys)]
        fig, ax = plt.subplots()
        n = min(len(x),min([len(g) for g in ys]))
        for y,t in zip(ys,titles):
            ax.plot(x[:n],y[:n],label=t)
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

    full = ReadFile(w.file).full
    if not full:
        print('Fatal: no valid data: file: {:}'.format(w.file))
        return

    bad = False
    if full['raw']:
        bad = True
        print('Lines in wrong format: not a number:')
        for l in full['raw']:
            print('  -> ',l)
    if full['bad']:
        bad = True
        print('Lines not in the same length:')
        for l in full['bad']:
            print('  -> ',l)
    if bad:
        t = input('Want to continue? [y]/N: ')
        if not (not len(t.split()) or t.strip().lower() == 'y'):
            print('You have decided to quit, nothing will be ploted')
            return

    if w.interactive:
        keys = full['keys']
        n = full['size']
        ask = None
        xaxis = begin = 0
        xdata = full['x']
        outputs = []
        for i in range(n):
            l = '{:2}  ->  {:15}  {:15.3f}  {:15.3f}  {:15.3f}'.format(i, *keys[str(i)])
            outputs.append(l)
        while True:
            while True:
                print('\nNote: to reset xaxis, input: xaxis=i')
                print('Note: to exclude begin data, input: begin=n')
                if ask:
                    g = ask.split()
                    ask = None    # reset
                else:
                    print('\nNote: parsing keys: (index -> key  min  max  avg)')
                    print(f'  >>: total number of data points: {len(xdata)}')
                    for i in range(n):
                        if xaxis == i:
                            print(outputs[i] + f'  (x-axis) (begin={begin})')
                        else:
                            print(outputs[i])
                    t = input('Input group to be parsed (use number & space): ')
                    g = t.replace(',',' ').split()
                    if not g: continue

                bo = None
                for v in g:
                    if v.startswith('xaxis='):
                        try:
                            k = int(v.replace('xaxis=',' '))
                            if k < 0 or k > n: raise ValueError
                        except:
                            print('  -> wrong setting xaxis, no space is allowed')
                            bo = False
                        else:
                            xaxis = k
                            xdata = full[keys[str(k)][0]]
                            bo = True
                            break
                if bo is True:
                    g.remove(v)
                elif bo is False:
                    continue
                for v in g:
                    if v.startswith('begin='):
                        try:
                            begin = int(v.replace('begin=',' '))
                            if begin < 0: raise ValueError
                        except:
                            print('  -> wrong setting begin, no space is allowed')
                        else:
                            g.remove(v)
                            break       # hack, must break
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
                        select = g
                        titles = [keys[str(v)][0] for v in select]
                        break
                else:
                    continue
            if begin == 0:
                ys=[full[i] for i in titles]
            else:
                xdata = xdata[begin:]
                ys=[full[i][begin:] for i in titles]
                
            print('\nNote: parsing keys: (index -> key  min  max  avg)')
            for c,q in zip([xaxis,*select], [xdata,*ys]):
                q = q[begin:]
                vmin = min(q)
                vmax = max(q)
                vavg = sum(q) / len(q)
                m = keys[str(c)][0]
                l = '{:2}  ->  {:15}  {:15.3f}  {:15.3f}  {:15.3f}'.format(c,m,vmin,vmax,vavg)
                print(l + f'  (begin={begin})')
            
            begin = 0   # reset
            plot_xvg(x=xdata,ys=ys,titles=titles)

            ask = input('Continue? [y|numbers]/N : ')
            ask = ask.replace(',',' ').strip()
            if not ask or ask.lower() == 'y': continue
            bo = False
            for i in ask.split():
                if i[:5] in ['begin','xaxis']: continue
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
        for k in ['bad','raw','size','keys']: use.remove(k)
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



