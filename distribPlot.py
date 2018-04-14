import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_data(fname):
    t = []
    binCenters = []
    binVals = []
    with open(fname, 'r') as f:
        cents = []
        vals = []
        for line in f:
            line = line.strip().split();
            if line[0] == "TIME:":
                t.append(float(line[1]))
                binCenters.append(cents)
                binVals.append(vals)
                cents = []
                vals = []
            else:
                cents.append(.5*(float(line[0]) + float(line[1])))
                vals.append(float(line[2]))
    return (t, binVals, binCenters)

def plot_data(t, cents, vals, saveName='tmp.png'):
    numBins = len(cents[0])
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i, elem in enumerate(t):
        myX = [elem for _ in range(numBins)]
        myY = cents[i]
        myZ = vals[i]
        ax.plot(myX, myZ, zs=myY)
    ax.view_init(30, 0)
    ax.set_xlabel('time', labelpad=15)
    ax.set_ylabel(r'$\log(\rho)$', labelpad=15)
    ax.set_zlabel('Count', labelpad=15)
    plt.savefig(saveName)
    plt.show()

fname = 'hist_data.txt'
(t, c, v) = load_data(fname)
plot_data(t,c,v, saveName='histo.png')
