import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import skimage.segmentation
from scipy.signal import find_peaks
import sys

def integrated_avalanches(dat, interval=5, offset=0, bins=10, prom=.03):
    res = []
    offset = int(offset / interval)
    for i in range(offset,len(dat)-interval,interval):
        res.append(np.sum(dat[i:i+interval]))
    res = np.array(res)
    resLog = np.log(res[res > 0])
    peaks2, _ = find_peaks(resLog,prominence=prom)
    return (peaks2, resLog[peaks2])

def avalanche_stats(fname,cutoff=4):
    print("load data\n")
    dat = np.loadtxt("Results/rho_data_" + fname + ".txt")
    L = np.shape(dat)[1]
    dat = dat[:,int(float(L)/4.):int(3.*float(L)/4.)]
    print("segmenting\n")
    im_labeled, n_labels = skimage.measure.label((dat > 0),background=0,return_num=True)
    im_props = skimage.measure.regionprops(im_labeled)
    res = []
    for elem in im_props:
        a = elem.area
	b = elem.bbox
        if (a < cutoff):
            continue
        w = np.abs(b[0] - b[2])
        h = np.abs(b[1] - b[3])
        w *= .25
        h *= np.sqrt(1.5)
        p = np.sum(dat[b[0]:b[2],b[1]:b[3]])
        t = b[0]
        res.append((t,w,h,a,p))
    intDat = np.sum(dat,axis=1) / float(np.shape(dat)[1])
    print("peaks\n")
    (pkX,pkY) = integrated_avalanches(intDat,interval=50,offset=72000,bins=15,prom=.003)
    with open("segRes_" + fname + ".txt","w") as f:
        for r in res:
            r = list(r)
            r = ','.join([str(elem) for elem in r])
            f.write(r + '\n')
    with open("peakRes_" + fname + ".txt","w") as f:
        for (x,y) in zip(pkX,pkY):
            f.write(str(x) + ',' + str(y) + '\n')

if __name__ == "__main__":
    print("starting")
    # fname = raw_input("File name: ")
    fname = sys.argv[1]
    if fname == '':
        fname = raw_input("File name: ")
    avalanche_stats(fname)

# print("starting")
# avalanche_stats("N3e4")
