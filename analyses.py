from collections import namedtuple
import os.path
import numpy as np
from scipy import ndimage

import Measurements
def gaussian2d(rsquared, sigma):
    return np.exp(- rsquared / (2 * sigma**2))
def addGaussianBlurAdv(data2d, sigma, maxradius):
    #r = np.array([[2,1,2],[1,0,1],[2,1,2]])
    r = np.zeros((maxradius*2+1,maxradius*2+1))
    
    for ind, _ in np.ndenumerate(r):
        r[ind] = sum((i - maxradius)**2 for i in ind)
    gauss = gaussian2d(r, sigma)
    gauss = 1/np.sum(gauss) * gauss
    
    return ndimage.convolve(data2d, gauss, mode = "nearest")
def syncFilter(data, cutoff, maxradius):
    r = np.zeros((maxradius*2+1,maxradius*2+1))
    for ind, _ in np.ndenumerate(r):
        r[ind] = cutoff*np.sqrt(sum((i - maxradius)**2 for i in ind))
    sincfilter = np.sinc(r)
    sincfilter= 1/np.sum(sincfilter) * sincfilter
    return ndimage.convolve(data, sincfilter, mode = "nearest")

def addBoxBlur(data):
    blur = 1/9*np.array([[1,1,1],[1,1,1],[1,1,1]])
    return ndimage.convolve(data, blur, mode = "nearest")
def getColumnLine(data2d):
    ret = np.empty(data2d.shape[1])
    for i, column in enumerate(data2d.T):
        ret[i] = np.mean(column)
    return ret
def getRowLine(data2d):
    ret = np.empty(data2d.shape[0])
    for i, column in enumerate(data2d):
        ret[i] = np.mean(column)
    return ret

def averageTotal(data2d):
    return np.mean(data2d)

def plotData(reduced_measurements):
    colors = ["b", "g", "r", "c", "m", "y"]
    f1 = plt.figure()
    ax1 = f1.add_subplot("111")
    
    d = [(m, analyses.getColumnLine(m.data)) for m in reduced_measurements]
    g = groupMeasurementData(d)
    #d = [(m, analyses.averageTotal(m.data)) for m in reduced_measurements]

    ind = 0
    for n in g:
        i = g[n]
        c = colors[ind]
        added_to_leg = False
        for j in i:
            d = j[0]
            m = j[1]
            xval = np.array(range(len(d),0,-1)) / m.scale
            if added_to_leg:
                plt.plot(xval, d, color=c)
            else:
                added_to_leg = True
                ax1.plot(xval, d, label=os.path.split(n)[1], color=c)
        ind = (ind + 1) % len(colors)
    ax1.set_ylabel("K $[\cdot]$")
    ax1.set_xlabel("Horizontal position $[mm]$")

    ax1.legend()
    
    
    f2 = plt.figure()
    ax2 = f2.add_subplot("111")
    
    d = [(m, analyses.getRowLine(m.data)) for m in reduced_measurements]
    g = groupMeasurementData(d)
    #d = [(m, analyses.averageTotal(m.data)) for m in reduced_measurements]

    ind = 0
    for n in g:
        i = g[n]
        c = colors[ind]
        added_to_leg = False
        for j in i:
            d = j[0]
            m = j[1]
            xval = np.array(range(len(d),0,-1)) / m.scale
            if added_to_leg:
                plt.plot(xval, d, color=c)
            else:
                added_to_leg = True
                ax2.plot(xval, d, label=os.path.split(n)[1], color=c)

        ind = (ind + 1) % len(colors)
    ax2.set_ylabel("K $[\cdot]$")
    ax2.set_xlabel("Vertical position $[mm]$")
    ax2.legend()