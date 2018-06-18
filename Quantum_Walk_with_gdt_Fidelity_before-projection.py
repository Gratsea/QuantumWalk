#Quantum walk with gradient descent techniques with return function Fidelity with the desired state

#Quantum walk -- moving forward 06/05/18

import numpy as np
import cmath
from scipy import optimize

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
    with open('Fidelity_today.txt', 'a+') as f:
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
    listc=[]
    
    listSt.append (initial)
    
    #Quantum Walk Forward : Fidelity_today in a superpotion over (n+1) number of sites. This state satisfies the conditions of eq.8
    
    for j in range (0,n,+1) : 
        #print (j+1)
        #definition of C
        v1=np.array([[initial[0][0]],[initial[3][0]]])
        #print (v1)
        matrixC = np.zeros((2*k,2*k),dtype=complex)
     
        
        #define arbitary a
        a=cmath.pi
        c = np.array([[v1[0][0],-cmath.exp(-a*1j)*v1[1][0].conjugate()],[v1[1][0],cmath.exp(a*1j)*v1[0][0].conjugate()]])
        c = c/np.linalg.norm(c)
        listc.append(c)
        
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
    Phi=initial
    
    Phi_target = np.array([[ 0.35355339 +0.00000000e+00j],  [ 0.00000000 +0.00000000e+00j],        [ 0.70710678 +8.65956056e-17j],        [-0.35355339 +4.32978028e-17j],        [-0.35355339 +0.00000000e+00j],        [ 0.00000000 -8.65956056e-17j],        [ 0.00000000 +0.00000000e+00j],        [-0.35355339 +1.29893408e-16j]])
    
    
    f = open("Fidelity_today.txt","a+")
    f.write("Final-output")
    f.close()
    with open('Fidelity_today.txt', 'a+') as f:
        print (Phi,file=f)
    f.close()
    

    Fidelity = np.dot(Phi.transpose(),Phi_target)*np.dot(Phi.transpose(),Phi_target)

    f = open("Fidelity_today.txt","a+")
    f.write("1-Fidelity")
    f.close()
    with open('Fidelity_today.txt', 'a+') as f:
        print (1-Fidelity.real,file=f)
    f.close()
    print (1-Fidelity.real,z)    
    return (1-Fidelity.real)


x0= [0.2,0.2]

minimizer_kwargs = {"method": "BFGS"}
#ret = optimize.basinhopping(func,x0,niter=1 ,T=1.0, stepsize=0.5, minimizer_kwargs=minimizer_kwargs,interval = 5, niter_success = 6 )
ret = optimize.basinhopping(func,x0, minimizer_kwargs=minimizer_kwargs,niter=100, T=1.0, disp = True )
   