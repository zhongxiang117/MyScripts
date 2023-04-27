#!/usr/bin/env python3

def fNewton(x, fobj, dfobj, ddfobj, tol=None, niters=None):
    """
    Parameters:
        x(float)    :   initial guess
        fobj(func)  :   objective function to get f(x)
        dfobj(func) :   first-order derivatives of function fobj
        ddfobj(func):   second-order derivatives of function fobj
        tol(float)  :   tolerance for x, positive value
        niters(int) :   number of iterations
    """
    tol = 0.000001 if tol is None else tol
    niters = 20 if niters is None else niters
    for i in range(niters):
        dx = dfobj(x)
        if abs(dx) <= tol:
            break
        x = x - dx / ddfobj(x)
    return x, fobj(x)


def fGlobalNewton(x,fobj,dfobj,ddfobj,delta=None,sigma=None,tol=None,niters=None):
    """
    Parameters:
        x(float)    :   initial guess
        fobj(func)  :   objective function to get f(x)
        dfobj(func) :   first-order derivatives of function fobj
        ddfobj(func):   second-order derivatives of function fobj
        delta(float):   in range (0.0, 1.0)
        sigma(float):   in range (0.0, 0.5)
        tol(float)  :   tolerance for x, positive value
        niters(int) :   number of iterations
    """
    delta = 0.5 if delta is None else delta
    sigma = 0.3 if sigma is None else sigma
    tol = 0.000001 if tol is None else tol
    niters = 20 if niters is None else niters
    gm = 20
    for i in range(niters):
        dx = dfobj(x)
        if abs(dx) <= tol:
            break
        dk = -dx / ddfobj(x)
        for m in range(gm):
            left = fobj(x+pow(delta,m)*dk)
            right = fobj(x) + sigma * pow(delta,m) * dfobj(x) * dk
            if left <= right:
                break
        x = x + pow(delta,m) * dk
    return x, fobj(x)


if __name__ == '__main__':
    def fobj(x):
        """
        three roots in range (-2,0), (1,2), and (4,5)
        """
        return x*x*x - 5*x*x + 2*x + 5
    def dfobj(x):
        return 3*x*x - 10*x + 2
    def ddfobj(x):
        return 6*x - 10

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



