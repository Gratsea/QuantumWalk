#projection of an arbitary state and creation of target state

from scipy.optimize import fsolve
import numpy as np
import cmath
ite = 0


def system(d):
    m=3 #number of steps
    s=m-1
    #superposition over 4 sites -- at the paper mentioned as phi
    u = np.ones((m+1,1),dtype = complex)
    u /= np.linalg.norm(u)
    #system of equations
    equations = []
    
    x = np.zeros((m),dtype=complex)
    for k in range (0,s,1) :
        x[k] = d[0+2*k]+d[1+2*k]*1j
        print (x[k])
    
    '''
    x = np.zeros((m+1),dtype=complex)
    x[0] = 0.
    x[m] = u[m][0]
    for k in range (0,s,1) :
        x[k+1] = d[0+2*k]+d[1+2*k]*1j
        print (x[k])
    
    #derivation of the system of equations -- general
    for s in range (1,m,1):
        condition = 0.
        for i in range (1,s,1) :
            condition += (u[i-1][0].conjugate()-x[i-1].conjugate())*(  u[m-s+i-1][0] - x[m-s+i-1] ) + x[i].conjugate()*x[m-s+i]
        equations.append( condition.real )
        equations.append( condition.imag ) 
        
        
            if (i == 0) :
                #print (i)
                condition = u[i][0].conjugate()*(  u[i+1][0] - x[i] ) + x[i].conjugate()*x[i+1]
                print ("0",condition.real , condition.imag)
                #print (condition.real,condition.imag)
                equations.append( condition.real)
                equations.append( condition.imag)
            else :            
                #print (i)        
                condition = ( u[i][0].conjugate() - x[i-1].conjugate() )*(  u[i+1][0] - x[i] ) + x[i].conjugate()* u[m][0] 
                print ("1",condition)
                equations.append( condition.real )
                equations.append( condition.imag )    
                #print (condition.real,condition.imag)
    '''
    '''
    #derivation of the system of equations -- for 4 sites     
    condition1 =  u[0][0].conjugate()*(  u[1][0] - x[0] ) + x[0].conjugate()*x[1]
    #without double summation
    #condition2 = ( u[1][0].conjugate() - x[0].conjugate())*(  u[2][0] - x[1] ) + x[1].conjugate()*u[m][0] 
    #with double summation
    condition2 = ( u[0][0].conjugate() )*(  u[1][0] - x[0] ) + x[0].conjugate()*x[1] + ( u[1][0].conjugate() - x[0].conjugate())*(  u[2][0] - x[1] ) + x[1].conjugate()*u[m][0] 
    #print (condition1,condition2)
    '''
    
    condition1 =  u[0][0].conjugate()*(  u[2][0] - x[1] ) + x[0].conjugate()*u[m][0]
    condition2 = ( u[0][0].conjugate() )*(  u[1][0] - x[0] ) + x[0].conjugate()*x[1] + ( u[1][0].conjugate() - x[0].conjugate())*(  u[2][0] - x[1] ) + x[1].conjugate()*u[m][0] 
    
    equations.append( condition1.real )
    equations.append( condition1.imag )  
    equations.append( condition2.real )
    equations.append( condition2.imag )
    print (d)
    global ite
    ite += 1
    
    return equations


d0 = fsolve(system,[0.7,0.5,-1.7,1.5])
print ("final",d0)
#print (ite)
