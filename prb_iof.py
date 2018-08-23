# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd

#fnamei='HMTX/H36.txt'
#fnameo='numso.txt'

def writeHMatrix(fname, matx):
   # df=pd.DataFrame(data=matx.astype(int))
    df=pd.DataFrame(data=matx)
    df.to_csv(fname,sep=',', header=False, index=False)
    
def readHMatrix(fname):
    data=pd.read_csv(fname,sep=',', header=None )
    return data.as_matrix() #.tolist()

#H=readHMatrix(fnamei)
#writeHMatrix(fnameo,H)