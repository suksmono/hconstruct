# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 17:30:49 2018

@author: Suksmono
PRB_graph
"""
from prb_symbolic0 import vq2s, vs2q
import numpy as np
import networkx as nx

###############################################################################        
#------------------------------------
# define functions
#------------------------------------
def cleanupG(G, TRIAL, Nmin):
    # remove low connectivity nodes
    ndsG=set(G.nodes())
    #vH=[]
    tHex=list()
    rmvList=list()
    while (len(ndsG)>0):
        tHex=ndsG.pop()
        if (G.degree(tHex)<Nmin):
            try: 
                G.remove_node(tHex)
                rmvList.add(tHex)
            except:
                print('Node already removed')
    #remove corresponding entries in trial
    tTRIAL=TRIAL.copy()
    newTRIAL=set()
    while(len(tTRIAL)>0):
        lstTRIAL=tTRIAL.pop()
        setLstTRIAL=set(lstTRIAL)
        mm=0
        TF=True
        while (mm<len(rmvList) and TF):
            trmvList=rmvList[mm]
            print('mm=',mm,'->evaluated node:',trmvList)
            if trmvList in setLstTRIAL:
                TF=False
            #-- end if
            mm=mm+1
        #-- end while
        # print('TF:', TF,', vecTRIAL=', tuple(setLstTRIAL))
        if (TF):
            newTRIAL.add(tuple(setLstTRIAL))
        #--end-if
    #-- end-while
    return(G, newTRIAL)
###############################################################################        
    
#------------------------------------
# define functions
#------------------------------------
def updGStruct(G,M):
    # list all nodes into set
    ndsG=set(G.nodes())
    # construct a list of H-vectors
    vH=[]
    vHex=list()
    while (len(ndsG)>0):
        tHex=ndsG.pop()
        vHex.append(tHex)
        vH.append(hex2svec(tHex,M))        
    #create adjacency matrix
    # D=np.abs(np.matmul(sk, np.transpose(sk)))  
    A=np.abs(np.matmul(vH, np.transpose(vH)))
    G=nx.Graph()
    for m in range(0,len(vH)):
        for nn in range(m+1,len(vH)):
            if A[m][nn]==0:
                G.add_edge(vHex[m], vHex[nn])
            #--end-if
        #--end for nn
    #-- end for m
    return(G)
###############################################################################        

##################################################
# get info from list of clique in G
# > a clique of maximum length
# > vector of clique length
##################################################

def maxCLQ(sCLQ):
    if len(sCLQ)<1:
        return set()
    else:
        vIdx=list()
        for m in range(0,len(sCLQ)):
            vIdx.append(len(sCLQ[:][m]))
        # get max 
        idxMax=np.argmax(vIdx)
        myClique=set(sCLQ[:][idxMax])
        return myClique

def listSCLQ(sCLQ):
    if len(sCLQ)<1:
        return len(set())
    else:
        vIdx=list()
        for m in range(0,len(sCLQ)):
            vIdx.append(len(sCLQ[:][m]))
    return vIdx

##################################################
# IS an orthoset?
##################################################

def isOrthoSet(CLQ,M):
    sk=CLQ2svlist(CLQ,M)
    D=np.abs(np.matmul(sk, np.transpose(sk)))  
    N=len(sk)
    TF= (np.sum(D)-M*N) == 0 
    return TF #, D

def matxIndicator(CLQ,M):
    sk=CLQ2svlist(CLQ,M)
    D=np.abs(np.matmul(sk, np.transpose(sk)))  
    return D
       
#    return D
#
##################################################
        #sk=CLQ2sk(CLQ)
##################################################
def CLQ2svlist(CLQ,M):
    sk=[np.int_(np.ones(M)).tolist()]
    #while nn in range(0,len(CLQ)):
    listCLQ=list(CLQ)
    #while len(ttCLQ)>0:
    for mm in range(0,len(listCLQ)):
        #tndHex=ttCLQ.pop()
        tndHex=listCLQ[mm]
        sk.append(hex2svec(tndHex,M))
    return(sk)
##################################################
def sortCLQS(allCLQ):
    # find out order of all cliques in graph
    # input : all cliques [ [1,3], [2,4,5,7], [0] ]
    # output: sorted all cliques [ [2,4,5,7], [1,3], [0] ] and set
    #        [ (2,4,5,7), (1,3), (0) ]

    tLQ=np.shape(allCLQ)
    LQ=tLQ[0]
    vLQ=list()
    for m in range(0,LQ):
        vLQ.append(len(allCLQ[:][m]))
    #print('Cliques orders in graph:\n',vLQ)
    
    # get sorted index
    clqIdx=np.argsort(vLQ)
    # arrange allCLQ according to length
    #print('Descending sorted cliques')
    NQ=len(vLQ)
    sCLQ=[]
    #setCLQ=[]
    for m in range(0,NQ):
        tIdx=NQ-m-1
        #print(tIdx)
        sCLQ.append(allCLQ[:][clqIdx[tIdx]])
        #setCLQ.append(set(allCLQ[:][clqIdx[tIdx]]))
        #print(sCLQ)
    #
    return sCLQ #, setCLQ
##################################################
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
def hex2qvec(v, N):
    #bchar=bin(int(v, 16))[2:]
    bchar=format(int(v, 16), 'b').zfill(N)
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
def hex2svec(v,N):
    #bchar=bin(int(v, 16))[2:]
    bchar=format(int(v, 16), 'b').zfill(N)
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


