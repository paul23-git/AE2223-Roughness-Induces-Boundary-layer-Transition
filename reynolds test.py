# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:49:48 2015
Reynoldstest
@author: Roeland
"""

from Reynoldsnumber import *
import numpy as np

T=np.array([[500, 300], [510, 600]])
P=10000000.0
X = flowproperties(T,P)
print(X.reynolds)
