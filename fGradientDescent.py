#!/usr/bin/env python3

def fGradientDescent(x,fobj,dfobj,step=None,prec=None):
    """
    Parameters:
        x(float)    :   initial guess
        fobj(func)  :   objective function to get f(x)
        dfobj(func) :   first-order derivatives of function fobj
        step(int)   :   stepsize, default 0.001
        prec(float) :   precision for y, default 1e-6
    """
    step = 0.001 if step is None else step
    if prec is None:
        prec = step / 1000
    else:
        if prec < step*1000: prec = step / 1000
    dv = 1.0
    niter = 0
    while dv > prec:
        niter += 1
        v = fobj(x)
        xn = x - dfobj(x)*step
        vn = fobj(xn)
        dv = abs(vn-v)
        x = xn
    return xn, vn, niter


if __name__ == '__main__':
    def fobj(x):
        return x*x
    def dfobj(x):
        return 2*x
    print('test-1: = ',fGradientDescent(1.0, fobj, dfobj))
    print('test-2: = ',fGradientDescent(2.0, fobj, dfobj))
    print('test-3: = ',fGradientDescent(3.0, fobj, dfobj))
    print('test-4: = ',fGradientDescent(4.0, fobj, dfobj,0.00001))