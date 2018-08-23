# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 17:30:49 2018

@author: Suksmono
PRB_BACKTRACK
"""
from prb_symbolic0 import vq2s, vs2q

#-----------------------------------------
# s-vector -> hexa 
#-----------------------------------------
def svec2hex(v):
    tVal=0
    NC=len(v)
    for m in range(0,NC):
        tV=v[m]
        tQ=int(vs2q(tV))
        tVal=tVal+tQ*2**(NC-m-1)
        #print(tQ,tVal)
#    return hex(int(tVal)), int(tVal)
    return hex(int(tVal))

#-----------------------------------------
# Hexa -> q-vector
#-----------------------------------------
def hex2qvec(v):
    bchar=bin(int(v, 16))[2:]
    NC=len(bchar)
    tVec=[]
    for m in range(0,NC):
        tV=int(bchar[m])
        tVec.append(tV)
        #print(tVec)
    return tVec

#-----------------------------------------
# Hexa -> s-vector
#-----------------------------------------
def hex2svec(v):
    bchar=bin(int(v, 16))[2:]
    NC=len(bchar)
    tVec=[]
    for m in range(0,NC):
        tV=int(bchar[m])
        tV1=int(vq2s(tV))
        tVec.append(tV1)
        #print(tVec)
    return tVec

#-----------------------------------------
# q-vector -> Hexa
#-----------------------------------------
def qvec2hex(v):
    tVal=0
    NC=len(v)
    for m in range(0,NC):
        tQ=v[m]
        tVal=tVal+tQ*2**(NC-m-1)
        #print(tQ,tVal)
#    return hex(int(tVal)), (tVal)
    return hex(int(tVal))


