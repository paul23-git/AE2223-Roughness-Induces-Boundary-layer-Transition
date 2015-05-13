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

def getColumnLine(data2d):
    ret = np.empty(data2d.shape[1])
    for i, column in enumerate(data2d.T):
        ret[i] = np.mean(column)
    return ret