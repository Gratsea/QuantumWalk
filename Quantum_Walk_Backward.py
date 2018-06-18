#Quantum walk -- moving backward 06/05/18

import math
import numpy as np
import cmath

#define final state
n=4 #number of steps
k=n+1 #number of sites


initial = np.zeros((2*k,1),dtype=complex)
initial[0][0]=initial[1][0]=1./math.sqrt(2)

'''
final = np.zeros((2*k,1),dtype=complex) #2k number of total sites (up and down)
#final[0][0] = 1. #spin up 
#final[1][0] = 1. #spin down

u1 = np.ones((2,1),dtype=complex)/math.sqrt(4)
u2 = np.ones((2,1),dtype=complex)/math.sqrt(4)

final[0][0] = u1[0][0]
final[2][0]= u1[1][0]
final[2*k-3][0]= u2[0][0]
final[2*k-1][0]= -u2[1][0]


final = np.ones((8,1),dtype=complex)/math.sqrt(4)
final[2*k-2][0]= final[1][0] = 0.0 
final[2*k-1][0] = -final[2*k-1][0] 
final[3][0] = 0.0
final[4][0] =1.0

final = final/np.linalg.norm(final) '''

final = np.array([[ 0.015625 +0.00000000e+00j], [ 0.000000 +0.00000000e+00j], [ 0.031250 +3.82702125e-18j], [-0.015625 +1.91351062e-18j], [-0.015625 +0.00000000e+00j], [ 0.000000 -3.82702125e-18j], [ 0.000000 +0.00000000e+00j], [-0.015625 +5.74053187e-18j]])
#final = np.array([[ 0.00035438+0.00000000e+00j], [ 0.        +0.00000000e+00j], [ 0.00141753+1.73597209e-19j], [-0.00035438+4.33993023e-20j], [ 0.        +2.60395814e-19j], [-0.00070876-8.67986045e-20j], [ 0.        -1.73597209e-19j], [ 0.00070876-8.67986045e-20j], [-0.00035438+8.67986045e-20j], [-0.00070876+8.67986045e-20j], [ 0.        +0.00000000e+00j], [-0.00035438+2.16996511e-19j],[ 0.000+0.0000j],[ 0.000+0.0000j],[ 0.000+0.0000j],[ 0.000+0.0000j]])

#example of final state with 4 sites - n=3
final = np.array([[ 0.0625+0.00000000e+00j], [ 0.    +0.00000000e+00j], [ 0.125 +1.53080850e-17j], [-0.0625+7.65404249e-18j], [-0.0625+6.35179545e-34j], [ 0.    -1.53080850e-17j], [ 0.    +0.00000000e+00j], [-0.0625+2.29621275e-17j], [ 0.    +0.00000000e+00j], [ 0.    +0.00000000e+00j]])

#example of final state with 4 sites - n=3
final = np.array([[ 0.0625+0.00000000e+00j], [ 0.    +0.00000000e+00j], [ 0.1875+2.29621275e-17j], [-0.0625+7.65404249e-18j], [-0.0625+1.53080850e-17j], [-0.0625-1.53080850e-17j], [ 0.0625-7.65404249e-18j], [ 0.0625+7.65404249e-18j], [ 0.    +0.00000000e+00j], [ 0.0625-3.06161700e-17j]])

Final = final
print ("Final")
print (Final)



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
listc = []
listinvc = []

listSt.append (Final)

for j in range (n,2,-1) : 
    print (j)
    #definition of C
    v1=np.array([[final[0][0]],[final[3][0]]])
    v2=np.array([[final[2*k-4][0]],[final[2*k-1][0]]])    
    print (v1)
    print (v2)
    matrixC = np.zeros((2*k,2*k),dtype=complex)
 
    
    a=cmath.pi
    c = np.array([[v1[0][0],cmath.exp(a*1j)*v2[0][0].conjugate()],[v1[1][0],cmath.exp(a*1j)*v2[1][0].conjugate()]])
    #c = np.array([[v1[0][0],-cmath.exp(-a*1j)*v1[1][0].conjugate()],[v1[1][0],cmath.exp(a*1j)*v1[0][0].conjugate()]])
    '''if (j==2) :
        print (c/np.linalg.norm(c))'''
    c/=np.linalg.norm(c)
    listc.append(c)
    invc = np.linalg.inv(c)
    listinvc.append(invc)
    
    invc = invc / np.linalg.norm(invc)
    for i in range (0,2*k,2):
        matrixC[0+i][0+i] = invc[0][0]
        matrixC[1+i][1+i] = invc[1][1]
        matrixC[0+i][1+i] = invc[0][1]          
        matrixC[1+i][0+i] = invc[1][0]
    matrixC /= np.linalg.norm(matrixC)
   
    listC.append (matrixC)    
    
    
    m1 = np.dot(invS,final)
    #♦print (m1)
    m2 = np.dot(matrixC,m1) 
    print (m2)
    m2 /= np.linalg.norm(m2)
    #print ("2",  prevstate )
    #print ("after step j=",j," ",prevstate)
    listSt.append (m2)
    final = m2
    
    '''
    x1 = np.array([[final[0][0]],[final[3][0]]])
    x2 = np.array([[final[2][0]],[final[5][0]]])
    print (np.dot(x1.transpose(),x2)) '''

    
C1 = np.dot(final,initial.transpose())
listC.append (C1)
listSt.append(initial)


#demonastration



