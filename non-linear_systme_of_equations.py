#projection of an arbitary state and creation of target state

from scipy.optimize import fsolve
import numpy as np


def system(d):
    m=3 #number of steps
    s=m-1
    
    #superposition over 4 sites -- at the paper mentioned as phi
    u = np.ones((m+1,1),dtype = complex)
    u /= np.linalg.norm(u)
    

    #set of parameters to be determined
    x = np.zeros((m+1),dtype=complex)
    x[0] = 0.
    x[m] = u[m][0]
    for k in range (0,s,1) :
        x[k+1] = d[0+2*k]+d[1+2*k]*1j
        print (x[k])
        
    #system of equations
    equations = []
    
    #derivation of the system of equations -- general
    for s in range (1,m,1):
        condition = 0.
        for i in range (1,s+1,1) :
            condition += (u[i-1][0].conjugate()-x[i-1].conjugate())*(  u[m-s+i-1][0] - x[m-s+i-1] ) + x[i].conjugate()*x[m-s+i]
        equations.append( condition.real )
        equations.append( condition.imag ) 
       
    return equations


d0 = fsolve(system,[0.7,0.5,-1.7,1.5])
print ("final",d0)
#print (ite)
