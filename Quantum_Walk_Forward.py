#Quantum walk -- moving forward 06/05/18

import math
import numpy as np
import cmath

#define final state
k=5 #number of sites

initial = np.zeros((2*k,1),dtype=complex)
initial[0][0]=initial[1][0]=1./math.sqrt(2)
#initial[0][0]=1.0


Initial = initial
print (Initial)


#definition of invS
invS = np.zeros((2*k,2*k),dtype=complex)
matrixS = np.zeros((2*k,2*k),dtype=complex)
for i in range (0,2*k,2) :
    invS[0+i][0+i] =1.
    matrixS[0+i][0+i] =  1.
    if (i+3)< 2*k :
        invS[1+i][3+i] = 1. #S-1
        matrixS[3+i][1+i] = 1.
        


n=4 #number os steps
listSt = []
listC = []

listSt.append (initial)

for j in range (0,n,+1) : 
    print (j+1)
    #definition of C
    v1=np.array([[initial[0][0]],[initial[3][0]]])
    #print (v1)
    matrixC = np.zeros((2*k,2*k),dtype=complex)
 
    
    #a=1.0
    a=cmath.pi
    c = np.array([[v1[0][0],-cmath.exp(-a*1j)*v1[1][0].conjugate()],[v1[1][0],cmath.exp(a*1j)*v1[0][0].conjugate()]])
    c = c/np.linalg.norm(c)
    '''if (j==1) :
        print (c)'''
    
    for i in range (0,2*k,2):
        matrixC[0+i][0+i] = c[0][0]
        matrixC[1+i][1+i] = c[1][1]
        matrixC[0+i][1+i] = c[0][1]          
        matrixC[1+i][0+i] = c[1][0]   
    #matrixC = matrixC/np.linalg.norm(matrixC)
     
    listC.append (matrixC)    
    
    
    m1 = np.dot(matrixC,initial)
    #print ("m1")
    #print (m1)
    m2 = np.dot(matrixS,m1)  
    #print ("m2-1")
    print (m2)
    #m2 /= np.linalg.norm(m2)
    #print ("m2-2")
    #print (m2)
    
    #print ("2",  prevstate )
    #print ("after step j=",j," ",prevstate)
    listSt.append (m2)
    initial = m2
    '''
    x1 = np.array([[initial[0][0]],[initial[3][0]]])
    x2 = np.array([[initial[2][0]],[initial[5][0]]])
    print (np.dot(x1.transpose(),x2))'''



#demonastration


