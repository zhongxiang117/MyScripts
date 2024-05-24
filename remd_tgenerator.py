#!/bin/env python3

import os
import sys
import math
import argparse
import configparser


FEATURES = [
    'version 0.1.0  : Replica Exchange Molecular Dynamics Temperature Generator',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


parameters = {
    'pdes':     (0.3,   float,  'Probability density'),
    'tlow':     (None,  float,  'Lower temperature'),
    'thigh':    (None,  float,  'Higher temperature'),
    'nw':       (None,  int,    'Number of water molecules'),
    'np':       (None,  int,    'Number of protein atoms'),
    'tol':      (0.0001,float,  'Tolerance'),
    'pc':       (1,     int,    'Protein constraints, 0: flexible, 1: hydrogen bonds only, 2: rigid'),
    'wc':       (1,     int,    'Water constraints, 0: flexible, 2: angles only, 3: rigid'),
    'hff':      (0,     int,    'Whether constraints on all or polar hydrogens, 0: all, 1: polar'),
    'vs':       (0,     int,    'Whether virtual site constraints exist, 0: no, 1: yes'),
}


def myintegral(m12,s12,cc):
    # for equation 7
    umax = m12 + 5*s12
    du = umax / 100
    u = 0
    tot = 0.0
    div = 2 * s12 * s12
    while u < umax:
        ui = u + du/2
        a = -cc*ui - (ui-m12)*(ui-m12)/div
        tot += math.exp(a)
        u += du
    return du*tot/(s12*math.sqrt(2*math.pi))


def get_remd_temperatures(
    pdes=None,tlow=None,thigh=None,nw=None,np=None,tol=None,
    pc=None,wc=None,hff=None,vs=None, **kws
):
    a0 = -59.2194
    a1 = 0.07594
    b0 = -22.8396
    b1 = 0.01347
    d0 = 1.1677
    d1 = 0.002976
    maxiter = 100
    kb = 0.008314

    npp = 0
    nprot = 0
    vc = 0
    if hff == 0:
        nh = round(np*0.5134)
        if vs == 1:
            vc = round(1.91*nh)
        nprot = np
    else:
        npp = round(np/0.65957)
        nh = round(np*0.22)
        if vs == 1:
            vc = round(np+1.91*nh)
        nprot = npp

    if pc == 1:
        nc = nh
    elif pc == 2:
        nc = np
    else:
        nc = 0

    ndf = (9-wc)*nw + 3*np - nc - vc
    flexener = 0.5*kb*(nc+vc+wc*nw)

    info = '\n'.join([
        'probability                    : {:}'.format(pdes),
        'temperature range              : {:} ~ {:} (K)'.format(tlow, thigh),
        'number of water molecules      : {:}'.format(nw),
        'number of protein atoms        : {:}'.format(np),
        'number of hydrogens in proteins: {:}'.format(nc),
        'number of hydrogens            : {:}'.format(npp),
        'whether use virtual site       : {:}'.format('yes' if vc else 'no'),
        'degree of freedom              : {:}'.format(ndf),
        'energy loss due to constraints : {:} (kJ/mol/K)'.format(flexener),
    ])

    fn_mu = lambda temp: (a0+a1*temp)*nw + (b0+b1*temp)*nprot - temp*flexener
    fn_sig = lambda temp: math.sqrt(ndf) * (d0+d1*temp)

    tlist = [(tlow,fn_mu(tlow),fn_sig(tlow),0.0,0.0,pdes), ]
    while tlow < thigh:
        piter = 0
        forward = 1
        iter = 0
        t1 = tlow
        t2 = min(t1+1, thigh)
        low = t1
        high = thigh
        while True:
            iter += 1
            mu12 = (t2-t1) * ((a1*nw)+(b1*nprot)-flexener)
            cc = (1/kb) * (1/t1 - 1/t2)
            var = ndf*(d1*d1*(t1*t1+t2*t2)+2*d1*d0*(t1+t2)+2*d0*d0)
            sig12 = math.sqrt(var)
            erfarg1 = mu12 / (sig12*math.sqrt(2))
            i1 = 0.5 * math.erfc(erfarg1)
            i2 = myintegral(mu12,sig12,cc)
            piter = (i1+i2)
            if piter > pdes:
                if forward == 1:
                    t2 = t2 + 1.0
                elif forward == 0:
                    low = t2
                    t2 = low + (high-low)/2
                if t2 > thigh:
                    t2 = thigh
            elif piter < pdes:
                if forward == 1:
                    forward = 0
                    low = t2 - 1.0
                high = t2
                t2 = low + (high-low)/2
            if not (abs(pdes-piter)>tol and iter<maxiter):
                tlist.append((t2,fn_mu(t2),fn_sig(t2),mu12,sig12,piter))
                break
        tlow = t2
    return info,tlist


def parsecmd():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            usage=f"""
xREMD Temperature Generator {VERSION}:

source: https://virtualchemistry.org/remd-temperature-generator/
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
        help='configuration file'
    )
    parser.add_argument(
        '-t','--template',
        dest='template',
        action='store_true',
        help='get template configuration file'
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

    if w.template:
        txt = ''
        for k,v in parameters.items():
            c = v[-1]
            t = '' if v[0] is None else v[0]
            txt += '# {:}\n{:6} = {:}\n\n'.format(c,k,t)
        print('Note: writing to template: remd_tgenerator.txt')
        with open('remd_tgenerator.txt','wt') as f:
            f.write(txt)
        return

    with open(w.file,'rt') as f:
        txt = f.read()
    txt = '[root]\n' + txt
    parser = configparser.ConfigParser()
    parser.read_string(txt)
    root = parser['root']
    unknown = set(root.keys()).difference(parameters.keys())
    if unknown:
        print('Warning: unknon keys: ',list(unknown))
    miss = set(parameters.keys()).difference(root.keys())
    if miss:
        print('Fatal: missing keys: ',list(miss))
        return
    new = {}
    for k,v in parameters.items():
        g = root[k]
        try:
            g = v[1](g)
        except ValueError:
            print('Fatal: wrong parameter type: {:}'.format(g))
        else:
            new[k] = g
    info,tlist = get_remd_temperatures(**new)
    print()
    print('  n     T           u         sig        D12u    D12sig    D')
    for i,v in enumerate(tlist):
        print('{:3} {:8.2f} {:12.1f}  {:8.2f}  {:8.2f}  {:8.2f}  {:5.3f}'.format(i+1,*v))
    print()
    print(info)
    print('\nTemperature Series:')
    print('  '.join([str(round(i[0],2)) for i in tlist]))
    print()

if __name__ == '__main__':
    parsecmd()


