#Quantum walk with gradient descent techniques with return Renyi entropy

#Quantum walk -- moving forward 06/05/18

import numpy as np
import cmath
from scipy import optimize
import math

def tensor(vectorA,vectorB) :
    m = np.size(vectorA,0)
    n = np.size(vectorB,0)
    tens=np.zeros((m,n))
    for i in range(m) :
        for j in range(n) :
            tens[i][j] = vectorA[i]*vectorB[j]
    return (tens);


def func(z) :    
    n=3 #number of steps
    k=n+1 #number of sites at the final state
    
    initial = np.zeros((2*k,1),dtype=complex)
    initial[0][0]= z[0]
    initial[1][0]= z[1]
    initial/= np.linalg.norm(initial)
    
    Initial = initial
    #print (Initial)
    
    f = open("0.txt","a+")
    f.write("Initial")
    f.close()
    with open('results_niter10_T1_NORM.txt', 'a+') as f:
       print (Initial,file=f)
    f.close()
   
    
    #definition of invS
    invS = np.zeros((2*k,2*k),dtype=complex)
    matrixS = np.zeros((2*k,2*k),dtype=complex)
    for i in range (0,2*k,2) :
        invS[0+i][0+i] =1.
        matrixS[0+i][0+i] =  1.
        if (i+3)< 2*k :
            invS[1+i][3+i] = 1. #S-1
            matrixS[3+i][1+i] = 1.
    
    listSt = []
    listC = []
    
    listSt.append (initial)
    
    #Quantum Walk Forward : results_niter10_T1_NORM in a superpotion over (n+1) number of sites. This state satisfies the conditions of eq.8
    
    for j in range (0,n,+1) : 
        #print (j+1)
        #definition of C
        v1=np.array([[initial[0][0]],[initial[3][0]]])
        #print (v1)
        matrixC = np.zeros((2*k,2*k),dtype=complex)
     
        
        #define arbitary a
        a=cmath.pi/2.
        c = np.array([[v1[0][0],-cmath.exp(-a*1j)*v1[1][0].conjugate()],[v1[1][0],cmath.exp(a*1j)*v1[0][0].conjugate()]])
        c = c/np.linalg.norm(c)
    
        
        for i in range (0,2*k,2):
            matrixC[0+i][0+i] = c[0][0]
            matrixC[1+i][1+i] = c[1][1]
            matrixC[0+i][1+i] = c[0][1]          
            matrixC[1+i][0+i] = c[1][0]   
        #matrixC = matrixC/np.linalg.norm(matrixC)
         
        listC.append (matrixC)    
        
        
        m1 = np.dot(matrixC,initial)
        m2 = np.dot(matrixS,m1)   #previous state
        #print (m2)
        listSt.append (m2)
        initial = m2/np.linalg.norm(m2)
    
    #initial now is the output--the final-- at the paper mentioned as Phi 
    phi=initial
    
    #output final state--superposition over (n+1) sites -- at the paper mentioned as phi 
    final = np.zeros((n+1,1),dtype = complex) 
    
    for i in range (0,n,1):
        final[i][0] = phi[2*i][0]+phi[2*i+1][0]
        
    final /= np.linalg.norm(final)
    #print ("final",final)
    
    
    f = open("results_niter10_T1_NORM.txt","a+")
    f.write("Final-output")
    f.close()
    with open('results_niter10_T1_NORM.txt', 'a+') as f:
        print (final,file=f)
    f.close()
    

    psiA, l, psiB = np.linalg.svd(final,full_matrices=1) #decomposition of initial matrix
    
    NORM=0.0
    p=1.0 # p has to be larger or equal than 1 for the algorithm to work
    
    m = np.size(final,0)  #number of rows of initial matrix
    n = np.size(final,1)  #number of columns of initial matrix
    sum=np.zeros((m,n))
    
    for i in range(k) :
         tens= tensor(psiA[i],psiB[i])
         sum += math.pow(l[i],p-1)*tens
         NORM = NORM + math.pow(l[i],p) 

    NORM = math.pow(NORM,1./p)

    
    #print ("fidelity",Fidelity)
    #print (z)
    f = open("results_niter10_T1_NORM.txt","a+")
    f.write("NORM")
    f.close()
    with open('results_niter10_T1_NORM.txt', 'a+') as f:
        print (NORM,file=f)
    f.close()
    print ()
    return (-NORM+math.sqrt(2))


x0= [0.2,0.5]

minimizer_kwargs = {"method": "BFGS"}
#ret = optimize.basinhopping(func,x0,niter=1 ,T=1.0, stepsize=0.5, minimizer_kwargs=minimizer_kwargs,interval = 5, niter_success = 6 )
ret = optimize.basinhopping(func,x0, minimizer_kwargs=minimizer_kwargs,niter=10, T=1.0, disp = True )
   