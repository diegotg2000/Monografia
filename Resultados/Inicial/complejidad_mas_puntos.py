import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import entropy

def entropia(A):
    return entropy(A,base=2)

def D(A,B):
    return entropia((A+B)/2)-entropia(A)/2-entropia(B)/2

def H(A):
    return entropia(A)/np.log2(len(A))

def Dstar(n):
    return -0.5*((n+1)/n*np.log2(n+1)+np.log2(n)-2*np.log2(2*n))

def complejidad(P):
    n=len(P)
    U=np.ones(len(P))/len(P)
    return D(P,U)*H(P)/Dstar(n)

NOMBRES=[]
BETAS=np.arange(10,51)/10
for i in range(len(BETAS)):
    NOMBRES.append('random'+str(BETAS[i])+'.BSKIndex')
for i in range(len(BETAS)):
    NOMBRES.append('norandom'+str(BETAS[i])+'.BSKIndex')
    
Srandom=np.zeros(len(BETAS))
Snorandom=np.zeros(len(BETAS))
Crandom=np.zeros(len(BETAS))
Cnorandom=np.zeros(len(BETAS))
Nrandom=np.zeros(len(BETAS))
Nnorandom=np.zeros(len(BETAS))

for i in range(len(NOMBRES)):
    A=np.loadtxt(NOMBRES[i], usecols=0)
    B=np.loadtxt(NOMBRES[i], usecols=1)
    C=np.concatenate((A,B),axis=0)
    
    nodos, apariciones = np.unique(C, return_counts=True)
    
    x,y=np.unique(apariciones, return_counts=True)
    
    R=y
    
    if i <= (len(BETAS)-1):
        puntos=np.loadtxt('random.dat',usecols=0)
        total=len(puntos)
        nonzero=len(nodos)
        p0=np.array([(total-nonzero)])
        P=np.concatenate((p0,R),axis=0)
        P=P/total
        Srandom[i]=entropia(P)
        Crandom[i]=complejidad(P)
        Nrandom[i]=len(P)
    else:
        puntos=np.loadtxt('norandom.dat',usecols=0)
        total=len(puntos)
        nonzero=len(nodos)
        p0=np.array([(total-nonzero)])
        P=np.concatenate((p0,R),axis=0)
        P=P/total
        Snorandom[i-len(BETAS)]=entropia(P)
        Cnorandom[i-len(BETAS)]=complejidad(P)
        Nnorandom[i-len(BETAS)]=len(P)
            
np.savetxt('Betas_mas_puntos.txt', BETAS)
np.savetxt('S_random_mas_puntos.txt',Srandom)
np.savetxt('S_norandom_mas_puntos.txt', Snorandom)
np.savetxt('C_random_mas_puntos.txt',Crandom)
np.savetxt('C_norandom_mas_puntos.txt', Cnorandom)
np.savetxt('numero_random.txt', Nrandom)
np.savetxt('numero_norandom.txt', Nnorandom)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    