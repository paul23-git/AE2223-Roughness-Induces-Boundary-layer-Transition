# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:19:24 2015

@author: Martijn
"""
import numpy as np
import math





#function for delta T
def deltaT_numeric(x, stepsize):
    #Input of variables
    T = 590
    P = 28 * 100000
    R = 287
    c = 1005
    k = 1.4
    q = 20

    #compute rho
    rho = P/(R*T)
    val = 0.1    
    return  np.nan_to_num((2*q*np.sqrt(x-val))/(np.sqrt(rho*c*k)*np.sqrt(pi))*np.sqrt(stepsize))


#Generating test data



#COOK AND FELDERMAN


def preCalcTimeDivisors(timelist, period, num):
    sqrt = np.sqrt
    rooted_period = sqrt(period)
    for n in range(len(timelist), num):
        
        divisor = 0;
        divlist = []
        last_sq = sqrt(n)
        for i in range(1, n +1):
            sq = sqrt(n-i)
            divlist.append(rooted_period * (sq+last_sq))
            last_sq = sq
        timelist.append(divlist)
    return timelist
                                                        


def CandF(dT, time, rho, c, k):
    sqrt = np.sqrt
    sumdata = np.zeros(len(time))
    for n in range(len(time)):
        tn = time[n]
        # precompute sqrt(tn-time[0])
        last_sq = sqrt(tn - time[0])
        sum_n = 0
        for j in range(1, n+1):
            sq_j = sqrt(tn - time[j])
            sum_n += dT[j] / (sq_j + last_sq)
            last_sq = sq_j
        sumdata[n] = sum_n    
    V = (2*sqrt(rho*c*k))/(sqrt(math.pi))
    return sumdata * V

def CandF_PreCalcedTime(dT, time, rho, c, k):
    sqrt = np.sqrt
    sumdata = np.zeros(len(time))
    for n in range(len(time)):
        t = time[-1]
        l = len(t)
        #print(t)
        # precompute sqrt(tn-time[0])
        sum_n = 0
        for j in range(1, n+1):
            #print(j, n, t)
            i = -n+(j-1)
            dtj = dT[j]
            sum_n += dtj / t[i]
        sumdata[n] = sum_n
    V = (2*sqrt(rho*c*k))/(sqrt(math.pi))
    return sumdata * V

def data_reduction_constant_time(dT, timestep, rho, cp, k):
    #print (len(dT), timestep)
    if timestep in data_reduction_constant_time.preCalcedPeriods:
        time = data_reduction_constant_time.preCalcedPeriods[timestep]
    else:
        time = []
        data_reduction_constant_time.preCalcedPeriods[timestep] = time
    
    if len(time) < len(dT):
        preCalcTimeDivisors(time, timestep, len(dT))

    return CandF_PreCalcedTime(dT, time, rho, cp, k )
data_reduction_constant_time.preCalcedPeriods = {}

def data_reduction(dT, time, rho, cp, k):
    return CandF(dT, time, rho, cp, k)

def plotIt(t, q):
    plt.plot(t,q)
    plt.ylabel('delta Temperature')
    plt.xlabel('Seconds [ms]')
    plt.title('Example Data')
    plt.legend()   
    plt.xlim([0,0.8])
    plt.show()



#Plotting example data
#plt.plot(t,dT)



if (__name__ == "__main__"):
    class data:
        def __init__(self):
            #normally lots of file reading
    
            self.integration = 400
            self.cols = 1
            self.rows = 1
            self.depth = 10
            self.readData()
    
        def readData(self):
            self.ml_delta_temp = 3 * np.random.random_sample((self.cols, self.rows, self.depth))
    
        def calcQML(self):
            print("--- Calculating q ---")
            print (self.cols, self.rows, self.depth)
            v = np.zeros([self.cols,self.rows,self.depth])
            T = 590
            P = 28 * 100000
            R = 287
            c = 1005
            k = 1.4
            q = 20
            rho = P/(R*T)  
            print(v)
    
            print("-------------")
            for x in range(self.cols ):
                print ("calculating column: ", x)
                for y in range(self.rows ):
                    v_tmp = data_reduction_constant_time(self.ml_delta_temp[x,y,:], 
                                                                    self.integration,
                                                                    rho, c, k)
                    v[x,y,:] = v_tmp
            print("--- q done ---")
            print(v)
    dat = data()
    dat.calcQML()
    #t = np.linspace(0,300,1000)
    #T = deltaT_numeric(t, 300/1000)
    #plotIt(t, T)