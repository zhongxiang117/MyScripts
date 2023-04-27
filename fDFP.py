#!/usr/bin/env python3

import numpy as np


def dfpmin(x,fobj,dfobj,niters=None,rho=None,sigma=None):
    """
    Parameters:
        x(list)     :   initial guesses
        fobj(func)  :   objective function to get f(x)
        dfobj(func) :   first-order derivatives of function fobj
        niters(int) :   maximum number of iterations
        rho(float)  :   
        sigma(float):   
    """
    niters = 100 if niters is None else niters
    rho = 0.55 if rho is None else rho
    sigma = 0.4 if sigma is None else sigma
    n = len(x)
    Hk = [[0. for i in range(n)] for i in range(n)]
    for i in range(n): Hk[i][i] = 1.0
    for i in range(niters):
        gk = dfobj(x)
        dk = [0. for i in range(n)]
        for t in range(n):
            dk[t] = sum([Hk[t][j]*gk[t] for j in range(n)])
        for m in range(20):
            vfold = fobj(x)
            vfnew = fobj([x[j]+(rho**m)*dk[j] for j in range(n)])
            res = sum([gk[j]*dk[j] for j in range(n)])
            if vfnew < vfold + sigma*(rho**m)*res:
                break
        xn = [x[j] + rho**m +dk[j] for j in range(n)]
        sk = [xn[j]-x[j] for j in range(n)]
        gkn = dfobj(xn)
        yk = [gkn[j]-gk[j] for j in range(n)]
        if sum([sk[j]*yk[j] for j in range(n)]) > 0:


if __name__ == '__main__':
    def fobj(x):
        """
        f(x) = 100 * (x1*x1 - x2)**2 + (x1 - 1)**2
        """
        return 100*(x[0]*x[0]-x[1])**2 + (x[0]-1)**2
    def dfobj(x):
        r1 = 400*x[0]*(x[0]*x[0]-x[1]) + 2*(x[0]-1)
        r2 = -200*(x[0]*x[0]-x[1])
        return [r1,r2]

    print('test-1 = ',fNewton(-2, fobj, dfobj, ddfobj))
    print('test-2 = ',fNewton(0, fobj, dfobj, ddfobj))
    print('test-3 = ',fNewton(1, fobj, dfobj, ddfobj))
    print('test-4 = ',fNewton(2, fobj, dfobj, ddfobj))
    print('test-5 = ',fNewton(4, fobj, dfobj, ddfobj))
    print('test-6 = ',fNewton(5, fobj, dfobj, ddfobj))

    print('\ntest-1 = ',fGlobalNewton(-2, fobj, dfobj, ddfobj))
    print('test-2 = ',fGlobalNewton(0, fobj, dfobj, ddfobj))
    print('test-3 = ',fGlobalNewton(1, fobj, dfobj, ddfobj))
    print('test-4 = ',fGlobalNewton(2, fobj, dfobj, ddfobj))
    print('test-5 = ',fGlobalNewton(4, fobj, dfobj, ddfobj))
    print('test-6 = ',fGlobalNewton(5, fobj, dfobj, ddfobj))

