import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import math as m
import growth


def logistic_growth(r, N0, K, t):
    return N0*K*np.exp(r*t)/(K+N0*(np.exp(r*t)-1))

def exponential_growth(r, N0, K, t):
    t0=np.min(t)
    start=logistic_growth(r, N0, K, t0)
    return start*np.exp(r*(t-t0))

def excess_bug_count():
    fig = plt.figure(1,figsize=(8,5))
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.5,3.5), ylim=(-5,1400))

    N0=1
    K=1000
    r=4.5
    t = np.arange(0.0, 3, 0.01)
    s = logistic_growth(r, N0, K, t)
    line, = ax.plot(t, s, lw=3, color='purple')

    growth_delta=np.max(t)/30
    for t0 in np.arange(0.0, 3, 0.2):
        tdelta=np.arange(t0, t0+growth_delta, 0.01)
        sdelta=exponential_growth(r, N0, K, tdelta)
        line, =ax.plot(tdelta, sdelta, lw=1, color="black")

    ax.annotate('no migration', xy=(0, 1),  xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )


    ax.annotate('lots of migration', xy=(1.5, 1100),  xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )

    ax.set_xlabel("Time [years]")
    ax.set_ylabel("Bugs [count]")
    ax.set_title("Exponential Growth from Logistic")

excess_bug_count()
plt.savefig("excess.pdf", format="pdf")
