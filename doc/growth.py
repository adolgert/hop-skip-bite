# This file makes sure that equations deriving
# the growth rate of bugs are correct.
import math
import logging
import numpy as np
import scipy.integrate
import scipy.optimize

logger=logging.getLogger(__file__)

# This is the modifier on the maximum hazard.
# Want to see if something small breaks it.
p=0.000035

def hazard(r, n0, k, t):
    return r*n0*n0*p*k*(math.exp(r*t)/(k+n0*(math.exp(r*t)-1)))**2

def u_of_t(r, n0, k, t):
    return math.exp(r*t)/((k/n0)-1)

def t_of_u(r, n0, k, u):
    return math.log(u*((k/n0)-1))/r

def inthazard(r, n0, k, t0, t1):
    u0=u_of_t(r, n0, k, t0)
    u1=u_of_t(r, n0, k, t1)
    return p*k*(
        math.log((u1+1)/(u0+1)) -
        (u1-u0)/((u0+1)*(u1+1))
        )

def invhazard(r, n0, k, t0, xa):
    u0=u_of_t(r, n0, k, t0)
    lhs=xa/(p*k) + math.log(u0+1) - u0/(u0+1)
    f=lambda x: math.log(x+1)-x/(x+1)-lhs
    xmax=math.exp(lhs+1)-1
    u1=scipy.optimize.brentq(f, u0, xmax)
    return t_of_u(r, n0, k, u1)


def invhazardvieta(r, n0, k, t0, xa):
    logger.error("This one fails")
    u0=u_of_t(r, n0, k, t0)
    # a=1
    b=u0-(1+u0)*xa/k
    c=-(1+u0)**2*xa/k
    x1=0.5*(-b-math.sqrt(b*b-4*c))
    x2=c/x1
    return t_of_u(r, n0, k, x2)

def invhazardriemann(r, n0, k, t0, xa):
    u0=u_of_t(r, n0, k, t0)
    delta_u=xa*(u0+1)**2/(p*k*u0)
    return t_of_u(r, n0, k, u0+delta_u)


def testuandt():
    r=.1
    n0=2.0
    k=1000.0
    maxerr=0
    for v in np.arange(0, 300, 0.2):
        u=u_of_t(r, n0, k, v)
        t=t_of_u(r, n0, k, u)
        logger.debug("{0} {1} {2}".format(v, u, t))
        maxerr=max(maxerr, abs(v-t))
    print("max error {0}".format(maxerr))


def testint():
    r=.1
    n0=2.0
    k=1000.0
    maxerr=0
    for t0 in np.arange(0, 20, 0.3):
        for t1 in np.arange(0, 100, 12.7):
            a=inthazard(r, n0, k, t0, t0+t1)
            b=scipy.integrate.quad(lambda x: hazard(r, n0, k, x), t0, t0+t1)[0]
            maxerr=max(maxerr, abs(a-b))
            logger.debug("{0}, {1}".format(a, b))
    print(maxerr)

def testinv():
    r=.1
    n0=2.0
    k=1000.0
    maxerr=0
    for t0 in np.arange(0, 20, 0.3):
        for t1 in np.arange(0.01, 100, 12.7):
            xa=inthazard(r, n0, k, t0, t0+t1)
            tf=invhazard(r, n0, k, t0, xa)
            maxerr=max(maxerr, abs(t0+t1 - tf))
            logger.debug("{0} {1} {2}".format(xa, t0+t1, tf))
    print(maxerr)

def testinvv():
    r=.1
    n0=2.0
    k=1000.0
    maxerr=0
    for t0 in np.arange(0, 20, 0.3):
        for t1n in range(4,9):
            t1=10**(-t1n)
            xa=inthazard(r, n0, k, t0, t0+t1)
            tf=invhazardriemann(r, n0, k, t0, xa)
            maxerr=max(maxerr, abs(t0+t1 - tf))
            logger.debug("{0} {1} {2}".format(xa, t0+t1, tf))
    print(maxerr)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    testuandt()
    testint()
    testinv()
    testinvv()
