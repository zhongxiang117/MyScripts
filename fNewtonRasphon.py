#!/usr/bin/env python3

def fNewtonRasphon(x, fobj, dfobj, tol=None, prec=None, iter=None):
    """
    Parameters:
        x(float)    :   initial guess
        fobj(func)  :   objective function to get f(x)
        dfobj(func) :   first-order derivatives of function fobj
        tol(float)  :   tolerance for x, positive value
        prec(float) :   precision for y, positive value
        iter(int)   :   iterations, if not defined, controlled by tol
    """
    tol = 0.000001 if tol is None else tol
    prec = 0.000001 if prec is None else prec
    iter = 20 if iter is None or iter <= 0 else iter
    for i in range(iter):
        y = fobj(x)
        if abs(y-prec) < 0:
            return x,y
        xn = x - fobj(x) / dfobj(x)
        if abs(x-xn) < tol:
            return xn, fobj(xn)
        xold = x
        yold = y
        x = xn
    if i >= iter:
        return xold, yold
    return x, y


if __name__ == '__main__':
    def fobj(x):
        """
        three roots in range (-2,0), (1,2), and (4,5)
        """
        return x*x*x - 5*x*x + 2*x + 5
    def dfobj(x):
        return 3*x*x - 10*x + 2
    
    print('test-1 = ',fNewtonRasphon(-2, fobj, dfobj))
    print('test-2 = ',fNewtonRasphon(0, fobj, dfobj))
    print('test-3 = ',fNewtonRasphon(1, fobj, dfobj))
    print('test-4 = ',fNewtonRasphon(2, fobj, dfobj))
    print('test-5 = ',fNewtonRasphon(4, fobj, dfobj))
    print('test-6 = ',fNewtonRasphon(5, fobj, dfobj))



