#!/usr/bin/env python3

import os
import json
import itertools


FEATURES = [
    'version 1.0.0  : Nov 8th, 2024',
    'version 1.1.0  : conclude keep to sell',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


class Hedge:
    """Hedge for bid

    Args:
        inlist: [(value,percent,sell_percent), (percent,sell_percent), (percent) ...]
                -> e.g. [(80, 0.80, 0.76), (15,0.15), 0.3, ...]
        delta_percent: float
        sell_num: int
    """
    def __init__(
            self,inlist=None,
            delta_percent=None,percent_num=None,sell_num=None,
            *args,**kws
    ):
        self.inlist = inlist if inlist else []
        self.delta_percent = delta_percent if delta_percent else 0.05
        self.percent_num = percent_num if percent_num else 3
        self.sell_num = sell_num if sell_num else 50
        # after run
        self.sell = []
        self.keys = {}

        self.full_keys = []
        if os.path.isfile('./hedge-bak.json'):
            with open('./hedge-bak.json','r') as f:
                self.full_keys = json.load(f)
        
        if not self.inlist and self.full_keys:
            print('Saved keys:')
            for i,v in enumerate(self.full_keys):
                print('\n{:3}  {:}'.format(i,v))
            print()
            k = input('Which to use? ')
            k = k.split()
            if k and k[0].isdigit() and 0<=int(k[0])<len(self.full_keys):
                k = int(k[0])
                self.inlist = self.full_keys[k]['inlist']
            print('Warning: no inputs')
            k = input('Which to remove? ')
            if k and k[0].isdigit() and 0<=int(k[0])<len(self.full_keys):
                self.full_keys.pop(int(k[0]))
                self.save_json(force=True)
            print()

    def save_json(self,force=False):
        if not force and not self.keys: return
        if self.keys: self.full_keys.append(self.keys)
        with open('./hedge-bak.json','wt') as f:
            json.dump(self.full_keys,f)

    def sell_on(self,now):
        if not self.keys:
            print('Warning: not preprocessed: run `H.run()` first')
            return
        if not (isinstance(now,(tuple,list)) and len(now) == len(self.inlist)):
            print('Fatal: --wrong input')
            return
        for v in now:
            if not (isinstance(v,(float,int)) and v-1<0.001):
                print('Fatal: wrong input')
                return

        get = sum([self.keys['share'][i]*v for i,v in enumerate(now) if v>0])
        print('\nSell on: ')
        self._print_begin()
        line = '  >>>:'
        for v in now:
            if v>0:
                line += '    x   '
            else:
                line += '        '
        print(line)
        self._print_begin(now,begin='now')
        d = get - self.keys['base']
        print('  Get:  net={:.3f}  return={:.3f}'.format(d,d/self.keys['base']))

    def run(self):
        if not (self.inlist and isinstance(self.inlist,(tuple,list))):
            print(f'Fatal: wrong input inlist: {self.inlist}')
            return
        make = []
        bo = False
        for g in self.inlist:
            if not isinstance(g,(list,tuple)) or len(g)<=0 or \
            not (len(g)!=1 or isinstance(g[0],(float,int)) and 0.05<g[0]<0.95):
                print(f'Fatal: wrong input: --inlist: {self.inlist}: {g}')
                return
            make.append(list(g[:]))
            if len(g) == 1:
                make[-1].insert(0,-1)
                bo = True
        if bo: self.inlist = make
        for g in self.inlist:
            if not isinstance(g,(list,tuple)) or len(g)<2 or len(g)>3 or \
                not (isinstance(g[1],(float,int)) and 0.05<g[1]<0.95) or \
                not (len(g)<3 or isinstance(g[2],(float,int)) and 0.05<g[2]<0.95):
                print(f'Fatal: wrong input: -inlist: {self.inlist}: {g}')
                return
            if len(g) == 2:
                if g[0] > 5 or g[0] < 0: continue
                if 0.05<g[0]<0.95 and 0.05<g[1]<0.95:
                    bo = True
                    g.insert(0,-1)
                else:
                    print(f'Fatal: wrong input: inlist: {self.inlist}: {g}')
                    return
        if bo:
            v = max([t[0]/t[1] for t in self.inlist])
            if v < 10:
                for g in self.inlist:
                    g[0] = 100 * g[1]
            else:
                for g in self.inlist:
                    if g[0]<0.0: g[0] = v * g[1]

        self.inlist = sorted(self.inlist,key=lambda x:x[1],reverse=True)

        buys = []
        base = 0.0
        expect = 0.0
        for g in self.inlist:
            buys.append((g[1],g[0]/g[1]))
            if len(g) == 2:
                base += g[0]
            else:
                base += 2*g[0] - g[0]/g[1]*g[2]
            expect += g[0]/g[1]
        self.suggest(buys,base)
        expect = expect / len(buys) - base

        sell = sorted(self.sell,key=lambda x: x[2],reverse=True)    # format: [(p,idx,e), ...]
        self.sell = sell[:self.sell_num]
        self.print_suggestion_on_sell()

        cost = sum([t[0] for t in self.inlist])
        share = [t[0]/t[1] for t in self.inlist]
        percent = [t[1] for t in self.inlist]
        print(' Hold:',end='')
        use = '  {:6.3f}'*len(self.inlist)
        print(use.format(*[t[0] for t in self.inlist]),end='')
        print('  =>  base={:.2f}'.format(base))
        print('Share:',end='')
        use = '  {:6.2f}'*len(self.inlist)
        print(use.format(*share),end='')
        print('  =>  {:.2f}'.format(sum(share)))
        print('  Net:',end='')
        use = '  {:6.2f}'*len(self.inlist)
        print(use.format(*[t-cost for t in share]))
        print('  Now: expect={:.2f}  min={:.2f}  max={:.2f}'.format(expect,min(share),max(share)))
        print('    +:  cost={:.2f}  return={:.2f}  percent={:.2f}'.format(cost,expect/cost,sum(percent)))

        self.keys = {
            'inlist'        : self.inlist,
            'share'         : share,
            'price'         : [t[0] for t in self.inlist],
            'percent'       : percent,
            'percent_sell'  : [(t[2] if len(t)==3 else t[1]) for t in self.inlist],
            'percent_worn'  : [(t[1]-t[2] if len(t)==3 else 0.0) for t in self.inlist],
            'sell_num'      : self.sell_num,
            'percent_num'   : self.percent_num,
            'delta_percent' : self.delta_percent,
            'cost'          : cost,
            'base'          : base,
            'expect'        : expect,
        }
        self.save_json()

    def print_suggestion_on_sell(self):
        cost = sum([t[0] for t in self.inlist])
        print('    Sell when color, red increase, blue decrease')
        print(' sell:',end='')
        use = '  {:6.2f}'*len(self.inlist)
        print(use.format(*[t[0] for t in self.inlist]),end='')
        print('  >>    net      R')
        for g in self.sell:
            line = '    >>'
            for i,t in enumerate(g[0]):
                if i in g[1]:
                    if t-self.inlist[i][1]>0.001:
                        line += '\x1b[31m' + '  {:6.3f}' + '\033[0m'
                    elif abs(t-self.inlist[i][1])<0.001:
                        line += '  {:6.3f}'
                    else:
                        line += '\x1b[34m' + '  {:6.3f}' + '\033[0m'
                else:
                    line += '  {:6.3f}'
            print(line.format(*g[0]),end='')
            print('  -> {:8.3f}  {:.3f}'.format(g[2],g[2]/cost))
        self._print_begin()

    def _print_begin(self,vl=None,begin=None):
        if not vl: vl = [t[1] for t in self.inlist]
        if begin:
            print('{:>5}:'.format(begin),end='')
        else:
            print('begin:',end='')
        use = '  {:6.3f}'*len(vl)
        print(use.format(*vl))

    def suggest(self, buys, base):
        cnt = 0
        self.sell = []
        delta = [i*self.delta_percent for i in range(1,self.percent_num+1)]
        for k in range(1,len(buys)+1):
            for r in itertools.combinations(range(len(buys)),k):
                for d in delta:
                    for g in itertools.product(*[(-1,1) for i in range(len(r))]):
                        new = [t[0] for t in buys]   # deep copy
                        bo = False
                        for i,t in enumerate(r):
                            new[t] += d*g[i]
                            if new[t]<0.0 or new[t]>1.0:
                                bo = True
                                break
                        if bo or sum(new)<=0.9: continue

                        cnt += 1
                        if cnt > 1000000: return

                        # for sell
                        for v in range(1,len(new)+1):
                            for s in itertools.combinations(range(len(new)),v):
                                e = sum([buys[t][1]*new[t] for t in s])
                                if e-base <= 1.0: continue

                                self.sell.append((new,s,e-base))
                                if len(self.sell) > 10*self.sell_num:
                                    p = sorted(self.sell,key=lambda x: x[2],reverse=True)
                                    self.sell = p[:3*self.sell_num]


h = Hedge()
#h.inlist = [[440,0.20],[0.14,0.13],  [0.17,0.16],[0.08,0.07],[0.062,0.06]]
#h.inlist = [[0.375,0.35],[0.334,0.33],[0.186,0.184]]
#h.inlist = [[500,0.962],[19,0.035]]
h.run()

h.sell_on([-1.0, -1, 0.96])


