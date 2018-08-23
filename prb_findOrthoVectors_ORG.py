# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 08:26:04 2018
@author: Asus
"""
from sympy import *
from prb_symbolic0 import gen_csqmatrix_inc, formHks,\
rmvIdSquareNEW, isingCoeffs
import neal
import numpy as np
from prb_graph import sortCLQS

'''
###############################################################################
Convert finding lost vector to function
###############################################################################
'''
def findOrthoVectors(sk,NSWEEPS, NREADS):
    [NC, M]=np.shape(sk)
    NO= 1 # M-NC    # number of unknown vector
    # number of required qubits
    NQ=M*NO + M*int(NO*(NO-1)/2)
    
    '''
    -------------------------------------------------------------------------------
    1. Formulate Hks
    -------------------------------------------------------------------------------
    '''
    
    ss=symbols('s0:%d'%NQ)
    qq=symbols('q0:%d'%NQ)
    
    '''
    -------------------------------------------------------------------------------
    generate s and q matrix
    -------------------------------------------------------------------------------
     |<---problem-->|<------ ancillas ------>|
    -------------------------------------------------------------------------------
      s0  1  1  | * 
      s1 -1  1  | *
      s2  1  1  | *
      s3 -1  1  | * 
    --------------------------------------------------------------------------------------------------------------------------------------------------------------
    '''
    s,q=gen_csqmatrix_inc(M,sk)
    '''
    -------------------------------------------------------------------------------
    calculate Hks
    -------------------------------------------------------------------------------
    '''
    
    print('Calculating Hks ...')  
    Hks=formHks(NC+1,s)
    '''
    ---------------------------------------------------
    simplify by substition of all si**2 terms: si**2->1
    ---------------------------------------------------
    '''
    print('Substitution of si**2<-1 ...')
    Hks=rmvIdSquareNEW(Hks,ss)
    
    H2s=Hks
    #print('H2s:\n', H2s)
    
    '''
    ------------------------------------------------------------
    3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
    ------------------------------------------------------------
    '''
    print('Obtaining Ising coefficients ...')
    b, hi, Jij = isingCoeffs(H2s,NQ)
    print('b:',b, 'hi:', hi, 'Jij', Jij)
   
    # normalize coefficients
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    #
    b=b/maxCoeff
    '''
    -----------------------------------------------------------------------------
    convert the problem into Ising coefficients
    -----------------------------------------------------------------------------
    '''
    #in dictionary format
    h={0:0}
    J={(0,1):1}
    
    for m in range(0,len(hi)):
        h[m]=hi[m]
        for n in range (m+1,len(hi)):
            J[m,n]=Jij[m,n]
        
    '''
    -----------------------------------------------------------------------------
    4. SOLVE THE PROBLEM
    -----------------------------------------------------------------------------
    select a solver
    > dimod: ExaxtSolver
    > neal:  SimulatedAnnealingSampler
    '''
    #
    print('Solving the problem using neal  ...')
    solver=neal.SimulatedAnnealingSampler()
    #NSWEEPS=1*1*10*10*1000
    response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=NREADS)
    #
    vE=response.data_vectors['energy']
    aSol=response.samples_matrix
    #print('All energy',vE)
    #print('All configurations',aSol)
    # report all minimum energy
    minE=min(vE)
    idxVE=[i for i, e in enumerate(vE) if e == minE]
    lSol=aSol.tolist()
    vmSol=list()
    # report all ortho-vectors 
    for mm in range(0,len(idxVE)):
        tVect=lSol[:][mm]
        if (not( (tVect in vmSol) or (np.negative(tVect).tolist() in vmSol) )):
            vmSol.append(tVect)        
        #print('tVect',tVect, 'vmSol:', vmSol) ##DEBUG ##
    #idxMinE=np.argmin(vE)
    #tSol=aSol[idxMinE]
    #vSol=tSol[0]
    '''
    return solution, base energy
    '''
    #################################  
    # groundstate=-b, 
    # achived ground state=minE
    # all eigenvectors in vmSol
    return b, minE, vmSol
    #################################   
'''
###############################################################################
'''

# ---- try new method ---
# from sk,extract Jij by formula
# Jij= SUM vi*vj', hi=[0.0, 0., ..., 0]
# b=M*len(sk)
def IsingCoeffsNEW(sk):
    [NC, M]=np.shape(sk)
    hi=np.zeros(M)
    Jij=np.zeros([M,M])
    A=np.zeros([M,M])
    hi=[]
    for m in range(0,len(sk)):
        vi=sk[:][m]
        A=A+2*np.outer(vi, vi)
   
    #
    b=len(sk)*M
    # -- print result --
    #print(A)
    for m in range(0,M):
        hi.append(0. )
        for nn in range(m+1, M):
          Jij[m][nn]=A[m][nn]    
    # 
    return b, hi, Jij
#

def findOrthoVectorsNEW(sk,NSWEEPS, NREADS):

    '''
    ------------------------------------------------------------
    3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
    ------------------------------------------------------------
    '''
    print('Obtaining Ising coefficients ...')
    #b, hi, Jij = isingCoeffs(H2s,NQ)
    b, hi, Jij = IsingCoeffsNEW(sk)
    #print('b:',b, 'hi:', hi, 'Jij', Jij)

    # normalize coefficients
    aJij=np.abs(Jij)
    ahi=np.abs(hi)
    maxCoeff=np.max([np.max(ahi), np.max(aJij)])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    #
    b=b/maxCoeff
    '''
    -----------------------------------------------------------------------------
    convert the problem into Ising coefficients
    -----------------------------------------------------------------------------
    '''
    #in dictionary format
    h={0:0}
    J={(0,1):1}
    
    for m in range(0,len(hi)):
        h[m]=hi[m]
        for n in range (m+1,len(hi)):
            J[m,n]=Jij[m,n]
        
    '''
    -----------------------------------------------------------------------------
    4. SOLVE THE PROBLEM
    -----------------------------------------------------------------------------
    select a solver
    > dimod: ExaxtSolver
    > neal:  SimulatedAnnealingSampler
    '''
    #
    print('Solving the problem using neal  ...')
    solver=neal.SimulatedAnnealingSampler()
    #NSWEEPS=1*1*10*10*1000
    response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=NREADS)
    #
    vE=response.data_vectors['energy']
    aSol=response.samples_matrix
    #print('All energy',vE)
    #print('All configurations',aSol)
    # report all minimum energy
    minE=min(vE)
    idxVE=[i for i, e in enumerate(vE) if e == minE]
    lSol=aSol.tolist()
    vmSol=list()
    # report all ortho-vectors 
    for mm in range(0,len(idxVE)):
        tVect=lSol[:][mm]
        if (not( (tVect in vmSol) or (np.negative(tVect).tolist() in vmSol) )):
            vmSol.append(tVect)        
        #print('tVect',tVect, 'vmSol:', vmSol) ##DEBUG ##
    #idxMinE=np.argmin(vE)
    #tSol=aSol[idxMinE]
    #vSol=tSol[0]
    '''
    return solution, base energy
    '''
    #################################  
    # groundstate=-b, 
    # achived ground state=minE
    # all eigenvectors in vmSol
    return b, minE, vmSol
    #################################   
'''
###############################################################################
'''
