import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

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
    myMeans = []
    for i, elem in enumerate(t):
        myX = np.array([elem for _ in range(numBins)])
        myY = np.array(cents[i])
        myZ = np.array(vals[i])
        # myX = myX[myZ > -5.8]
        # myY = myY[myZ > -5.8]
        # myZ = myZ[myZ > -5.8]
        myY /= np.sum(myY)
        myY *= len(myY)
        myMeans.append(10. ** (np.dot(myY,np.array(myZ)) / np.sum(myY)))
        ax.plot(myX, myZ, zs=myY)
    # ax.plot(t, myMeans, [1e2 for _ in range(len(t))], color='black')
    ax.view_init(30, 0)
    ax.set_xlabel('time', labelpad=15)
    # ax.set_ylabel(r'$\log(\rho)$', labelpad=15)
    ax.set_ylabel(r'$\delta$', labelpad=15)
    ax.set_zlabel('Count', labelpad=15)
    plt.savefig(saveName)
    plt.figure()
    tmpBins = np.sum(cents, axis=0)
    plt.plot(vals[0], tmpBins)
    plt.xlabel(r'$\delta$')
    plt.ylabel('Count')
    plt.title('Gap Histogram - All Times')
    plt.figure()
    plt.plot(vals[0], cents[-1])
    plt.title('Gap Histogram - Last Time')
    plt.figure()
    # plt.semilogy(t,myMeans)
    # plt.show()

def plot_hmaps(a=0.1, dt=.25, saveStr='tmp_hmap', logScale=False, Gamma=-.84701, b=-1e-4, kappa=0.001, r=1e-3):
    dataP = np.loadtxt('rho_data.txt')
    dataH = np.loadtxt('h_data.txt')
    dataH = -(Gamma*b)*dataH + (1-Gamma)
    # dataP = -(Gamma*b)*dataP + (1-Gamma)
    for i, elem in enumerate(dataH):
        dataH[i] += dt*i*kappa
        # dataP[i] += dt*i*0.0003
    if logScale:
        dataP = np.log10(dataP + 1e-11)
        dataH = np.log10(dataH + 1e-11)
    plt.imshow(dataP, aspect=a)
    plt.colorbar()
    locs, labels = plt.yticks()
    totTime = dt * len(dataP)
    nDivs = 6
    newLocs = [i*len(dataP) / nDivs for i in range(nDivs+1)]
    newLabs = ['%.1f' % (l*dt) for l in newLocs]
    # nDivs = len(locs) - 2
    # newLabels = ['%.1f' % (totTime*(i-1)/nDivs) for i in range(1,len(locs))]
    # plt.yticks(locs[1:], tuple(newLabels))
    plt.yticks(newLocs, tuple(newLabs))
    plt.xlabel('Lattice Site')
    plt.ylabel('Time')
    # plt.title(r'$\Gamma = $' + str(Gamma) + ', b = ' + str(b) + r', $\kappa$ = ' + str(kappa))
    titleStr = r'$\kappa$: ' + str(kappa) + ', r: ' + str(r)
    plt.title(titleStr)
    plt.savefig(saveStr + '_rho.png')
    plt.figure()
    plt.imshow(dataH, aspect=a)
    plt.colorbar()
    plt.xlabel('Lattice Site')
    plt.ylabel('Time')
    plt.yticks(newLocs, tuple(newLabs))
    # plt.title(r'$\Gamma = $' + str(Gamma) + ', b = ' + str(b) + r', $\kappa$ = ' + str(kappa))
    plt.title(titleStr)
    plt.savefig(saveStr + '_h.png')
    plt.show()

fname = 'hist_data.txt'
# (t, c, v) = load_data(fname)
dat = np.loadtxt('hist_data.txt')
x = 0.5*(dat[:,0] + dat[:,1])
y = dat[:,2]
y /= np.sum(y)
y *= 100
x = x[y > 0]
y = y[y > 0]
plt.plot(x,y)
plt.xlabel(r'$\delta$')
plt.ylabel('Probability density')
titleStr = r'$\kappa$: 5e-4, r: 5e-4'
plt.title(titleStr)
plt.show()
# plot_data(t,c,v, saveName='histo.png')

dat2 = pd.read_csv('mean_data.txt', delimiter='\t', names=['t', 'm'])
t = np.array(dat2['t'])
m = np.array(dat2['m'])
t2 = t[m > 1e-8]
m2 = m[m > 1e-8]
#plt.plot(dat2['t'], dat2['m'])
plt.plot(t2,m2)
# plt.yscale('log')
# plt.xscale('log')
plt.xlabel('Time')
plt.ylabel(r'$\langle\delta\rangle$')
plt.title(titleStr)
plt.show()


# plot_hmaps(a=.03, Gamma=-0.9, b=-5e-3, kappa=5e-4, r=5e-4)
