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