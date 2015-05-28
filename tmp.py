import numpy as np
from scipy import interpolate
import numpy as np
import Data_Reduction
import weakref
from numpy import inf
import Stanton


def matlab_convert(Im, tint, Tcam, eps):
    print(tint, Tcam, eps)
    if(np.round(tint*1.e6) == 400):
        print("--- Converting (Matlab 400 microsec) ---")
        elambda=np.array([2498026.51787585,1784.71556672801,6.67171636375203e-12,3215.59874368947,493144.349437419])
    elif(math.round(tint*1.e6) == 200):
        print("--- Converting (Matlab 200 microsec) ---")
        elambda=np.array([436445.608980990,1408.69064523959,0.394428724501193,0.000889315394332730,186268.678717702])
    else:
        raise ValueError("Incorrect time integral")
    
    v1 =elambda[4]/(np.exp(elambda[1]/Tcam)-elambda[2])
    
    Tm=elambda[1]/np.log((eps*elambda[0]/((Im)-elambda[3]-elambda[4]/(np.exp(elambda[1]/Tcam)-elambda[2])))+elambda[2])
    print("--- Converting done ---")
    return Tm

def python_convert( Im): #Convert the raw data to temperatures
    print("--- Converting (python way) ---")
    dl_tab = [5700,6200,6683,7205,8363,9668,11110,12700];
    t_tab = [-5,0,5,10,20,30,40,50];
    tck = interpolate.splrep(dl_tab,t_tab,s=0)
    temp = np.zeros(Im.shape)
    for n in range(0,Im.shape[2]):
        temp_dum = interpolate.splev(Im[:,:,n].flatten(0),tck,der=0)
        temp[:,:,n] = np.reshape(temp_dum,(Im.shape[0],Im.shape[1]))
    temp += 273.15
    #self.calcDeltaT()
    print("--- Converting done ---")
    return temp



tmp = np.array([[[100000]]])
v1 = python_convert(tmp)
v2 = matlab_convert(tmp, 0.00004, 28+273.15, 0.86)
print(v1,v2)


