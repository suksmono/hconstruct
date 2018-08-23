"""
# -----------------------------------------------------------------------------
# Created: on 13-Aug-2018, updated 15-Aug-2018
# Author: Suksmono@{STEI-ITB, MDR Inc.}
# Summary: By reperesenting binary vectors as node, we search M-orthoset or
#          M-order clique by growing a random graph randomly
# -----------------------------------------------------------------------------
#    ALGORITHM
# -----------------------------------------------------------------------------
#(0) DEFINE 
#        -Problem size: H-Order M
#        -Problem effort: MAXITR, tolerance: eps, NSWEEP, NREADS
#(1). INITialization (from the start or continues from file)
#        -DISCRETE STRUCTURES
#            -Master graph: G<-{}
#            -Record of trial TRIAL<- {}
#            -List of cliques: LS_CLQ<- Clique(G), ordered by |tCLQ|
#                -Init master-Clique: CLQ <-MaxClique(tCLQ)
#        -PARAMS: NV<-|CLQ|, itr<-0
#(2) WHILE (NV< H-ORDER and itr<MAXITR):   
#    -From Largest to Smallest, fetch tCLQ in NLS_CLQ not yet in TRIAL
#        -Convert into matrix of variables: tCLQ->sk
#        -Get next ortho vectors: 
#            -update list of tried tCLQ: TRIAL<-TRIAL U tCLQ
#            -get vSol
#        - ----------------- replaced ------------------
#        -#Generate complete graph from F<-{vSol} U tCLQ
#        -#Combine F to master: G<-F U G
#        - ----------------- --------- ------------------
#        -Insert nodes into graph: G<-{vSol}
#        -Create adjacency matrix: A
#        -Update graph structure based on adjcy mtx: G<-fUpd(A)
#        - ----------------- not yet implemented ------------------
#        #Remove low -degree nodes from grap: G<- G-badNodes(G) 
#        #Remove list in trial incl bad nodes: TRIAL<-TRIAL-[badNodes(G)]
#        - ----------------- --------- ------------------
#    -Update List of Clique: LS_CLQ=[{tCLQ}]<-Clique(G)
#        -Update master clique: CLQ <- Max(tCLQ)
#    ## save to file
#    -save structures into file    
#    -itr=itr+1, NV=|CLQ|
# -----------------------------------------------------------------------------
# this verision found H28 at iteration 210 out of 280
# -----------------------------------------------------------------------------

"""

# -----------------------------------------------------------
# IMPORT PACKAGES AND MODULES
# -----------------------------------------------------------
from sympy import *
import numpy as np
from prb_findOrthoVectors_ORG import * #findOrthoVectorsNEW
from prb_iof import readHMatrix, writeHMatrix 
import networkx as nx
#from prb_graph import svec2hex, hex2svec, CLQ2svlist, sortCLQS, isOrthoSet
#from prb_graph_ORG import *
from prb_graph import *
from networkx.algorithms.approximation import clique as xclq
import matplotlib.pyplot as plt
import pickle as pk

#------------------------------------
# define functions
#------------------------------------
def updGStruct(G):
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
####
def genSHVector(M):
    rd=np.random.permutation(M)
    v=np.ones(M).tolist()
    for mm in range(0,int(M/2)):
        v[rd[mm]]=-1.
    #
    return(v)

# codelettes
#----------------------------------------------------------------------------
def updateMAXREADS(NSWEEPS, NV):
    #NTSWEEPS=NV*int(np.sqrt(NV))*NSWEEPS
    NTSWEEPS=NV*NSWEEPS
    return NTSWEEPS
    
'''
# ============================================================================
# MAIN PROGRAM: RANDOM GROWING OF H-GRAPH, OBTAIN LARGER CLIQUE
# ============================================================================
'''    
if __name__ == '__main__':
    '''
    # -----------------------------------------------------------
    # DEFINE PROBLEM
    # -----------------------------------------------------------
    '''
    # define H-order, initial parameters
    kk=3
    M = kk*4
    MAXITR=10*M
    # define number of sweeps/iteration
    NSWEEPS=1*1*1*M*10*1000
    # number of read or solutions in the NEAL
    NREADS=int(M/2)
    eps= 1e-9 #0.001/M
    
    '''
    # -----------------------------------------------------------
    # PARAMETERS INITIALIZATION
    # -----------------------------------------------------------
    '''  
    fName_CLQ   = r'BKTRACK/CLQ'+str(M)+'.pkl'
    fName_G     = r'BKTRACK/G'+str(M)+'.pkl'
    fName_TRIAL = r'BKTRACK/TRIAL'+str(M)+'.pkl'
    fName_LS_CLQ  = r'BKTRACK/LS_CLQ'+str(M)+'.pkl'
    '''
    # -----------------------------------------------------------
    # INITIALIZATION OF DISCRETE STRUCTURES
    # -----------------------------------------------------------
    '''
    np.random.seed(17)
    NDZ_svec=np.int_(np.ones(M)).tolist()
    ###############################################################
    F_SEED=0# 0: first time running, 1: continu old job
    NV_LIM=0 # simulate failure in finding M-clique, continu next 
    ###############################################################
    
    if F_SEED>0:    # continue previous job
        #
        CLQ     = pk.load(open(fName_CLQ, 'rb'))
        G       = pk.load(open(fName_G, 'rb'))
        TRIAL   = pk.load(open(fName_TRIAL, 'rb'))
        OLS_CLQ    = pk.load(open(fName_LS_CLQ, 'rb'))
    else:
        tvHex=svec2hex(genSHVector(M))
        # Init master Clique
        CLQ={''}
        CLQ.pop()
        #
        CLQ.add(tvHex)
        # Init master Graph
        G=nx.Graph()
        G.add_node(tvHex)
        # Init record of trial TRIAL<- {}
        TRIAL={()} #set()
        OLS_CLQ= [[tvHex]]
    # end-if-else --
    '''
    # ============================================================================
    # MAIN LOOP
    # ============================================================================
    '''
    itr=0
    vgapE=[]    # record energy
    vszG=[]
    vszCLQ=[]
    NV=len(CLQ)
    MAXBRANCH=int(NREADS)
    # init random seed for reproducability
    while ( (NV<M-1-NV_LIM) and (itr<MAXITR) ) :
        print('===================================================')
        print('>>> FINDING -', NV+1,'TH VECTOR OF -',M, 'ITER=', itr, 'of', MAXITR)
        print('===================================================')

        #NTSWEEPS=NV*NSWEEPS
        NTSWEEPS=updateMAXREADS(NSWEEPS, NV)
        
        # current status
        print('Est-CLQ(G):', len(xclq.max_clique(G)), ', |G|=', len(G.nodes()))
        # -------------------------------------------------------------------------
        # fetch a clique not yet in the trialset
        # -------------------------------------------------------------------------
        TF=True
        mm=0
        tCLQ=set()
        while(TF and (mm<len(OLS_CLQ)) ):
            #print('Trying:',mm, 'with |CLQ|',len(), ' of', len(LS_CLQ))
            tCLQ=OLS_CLQ[:][mm].copy()
            #print('tCLQ->', tCLQ)
            tpCLQ=tuple(set(tCLQ))  #force ordering of the nodes ID
            notInTrial=not(tpCLQ in TRIAL)
            enoughLength=len(tCLQ)>1
            TF= not(notInTrial and enoughLength and isOrthoSet(set(tCLQ),M))
            mm=mm+1
        # now we have qualified tCLQ
        TRIAL.add(tuple(set(tCLQ))) # add into list of tested clique
    
        # -------------------------------------------------------------------------
        # prepare for sampling, use sk = [1] U CLQ
        # -------------------------------------------------------------------------
        print('Extending N(CLQ):', len(tCLQ) )
        sk=CLQ2svlist(tCLQ,M)
        baseE, minE, vmSol=findOrthoVectorsNEW(sk,NTSWEEPS, NREADS)
        #baseE, minE, vmSol= findOrthoVectors(sk,NTSWEEPS, NREADS)
        gapE=abs(baseE+minE)
        print('baseE=', baseE, 'Energy gap=', gapE)
        vgapE.append(gapE)    
        # if orthogonal, add to orthoset sk
        if ( (gapE/M)<eps):
            #print('CLQ(G) before:', len(xclq.max_clique(G)), ', |G|=', len(G.nodes()))
            NSOL=min(MAXBRANCH,len(vmSol))
            for mm in range(0,NSOL):
                tndHex=svec2hex(vmSol[:][mm])
                G.add_node(tndHex)
            #--end for mm ---
            # update graph structure in G
            G=updGStruct(G)
            tndHex=svec2hex(vmSol[:][0])
            CLQ.add(tndHex)
            # create complete graph from CLQ
            print('Est-CLQ(G) after:', len(xclq.max_clique(G)), ', |G|=', len(G.nodes()))
            NV=len(sk)
    
        # generate list of cliques from the graph G
        NLS_CLQ=list(nx.clique.find_cliques(G))
        OLS_CLQ=sortCLQS(NLS_CLQ)
        CLQ=maxCLQ(OLS_CLQ).copy()
        NV=len(CLQ)
    
        '''
        # SAVE ALL PARAMS TO FILE: > CLQ, G, TRIAL, LS_CLQ
        '''
        pk.dump(CLQ, open(fName_CLQ, 'wb+'))
        pk.dump(G, open(fName_G, 'wb+'))
        pk.dump(TRIAL, open(fName_TRIAL, 'wb+'))
        pk.dump(OLS_CLQ, open(fName_LS_CLQ, 'wb+'))
        # --
        itr=itr+1
        vszG.append(len(G.nodes()))
        vszCLQ.append(NV)

        # show distribution of cliques in G
        num_bins=M
        lstCLQ=listSCLQ(NLS_CLQ)
        n, bins, patches = plt.hist(lstCLQ, num_bins, facecolor='blue', alpha=0.5)
        # Add title and axis names
        plt.title('Clique Distribution')
        plt.xlabel('Clique Size')
        plt.ylabel('Number of Cliques')
        plt.show()
        
    ############
    if (NV> M-2) and isOrthoSet(CLQ,M):
        print('Hadamard matrix of order',M, 'is found !!!')
    else:
        print('Hadamard matrix has not been found ...')
    
    '''
    # ============================================================================   
    # SAVE AND DISPLAY RESULTS
    # ============================================================================   
    '''
    
    # write h-matrix if the search is successful
    fmatxname='HMTX/H'+str(M)+'.txt'
    if isOrthoSet(CLQ,M):
        HMAT=CLQ2svlist(CLQ,M)
        writeHMatrix(fmatxname,HMAT)
        
    import matplotlib.pyplot as plt
    import pandas as pd
    
    mc=xclq.max_clique(G)
    print('Order of final CLQ:', len(CLQ))
    print('Estimated max-Clique in Graph:', len(mc))
    
    DD=matxIndicator(CLQ,M)
    imgplot = plt.imshow(abs(DD))
    plt.show()
    
    # draw graph ?
    nx.draw(G,with_labels=True, font_weight='bold')
    
    # plot curves
    fig, ax = plt.subplots()
    ax.plot(vszG,'r', label='Graph-size') 
    ax.plot(vszCLQ,'b', label='Max-Clique size')
    #ax.axis('equal')
    leg = ax.legend();
    plt.title('Growth of graph')
    plt.xlabel('Iteration (t)')
    plt.ylabel('Size') 
    plt.show()
    
    # display histogram of cliques
    #
    NBIN=10
    tbins=[]
    dbin=1
    for m in range (0,M,dbin):
        tbins.append(m)
    #
    sdat=listSCLQ(OLS_CLQ)
    n, bins, patches = plt.hist(sdat, num_bins, facecolor='blue', alpha=0.5)
    # Add title and axis names
    plt.title('Clique Distribution')
    plt.xlabel('Clique Size')
    plt.ylabel('Number of Cliques')
    plt.show()

    sdat1=listSCLQ(list(TRIAL))
    n, bins, patches = plt.hist(sdat1, num_bins, facecolor='blue', alpha=0.5)
    # Add title and axis names
    plt.title('Distribution of Attempts')
    plt.xlabel('Clique Size (before attempt)')
    plt.ylabel('Number of Attemps')
    plt.show()
