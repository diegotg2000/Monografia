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

NOMBRES=['random1.0.BSKIndex','random1.5.BSKIndex','random2.0.BSKIndex','random2.5.BSKIndex','random3.0.BSKIndex','random3.5.BSKIndex','random4.0.BSKIndex','random4.5.BSKIndex','random5.0.BSKIndex','norandom1.0.BSKIndex','norandom1.5.BSKIndex','norandom2.0.BSKIndex','norandom2.5.BSKIndex','norandom3.0.BSKIndex','norandom3.5.BSKIndex','norandom4.0.BSKIndex','norandom4.5.BSKIndex','norandom5.0.BSKIndex']
BETAS=[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]

Srandom=np.zeros(9)
Snorandom=np.zeros(9)
Crandom=np.zeros(9)
Cnorandom=np.zeros(9)

for i in range(len(NOMBRES)):
    A=np.loadtxt(NOMBRES[i], usecols=0)
    B=np.loadtxt(NOMBRES[i], usecols=1)
    C=np.concatenate((A,B),axis=0)
    
    nodos, apariciones = np.unique(C, return_counts=True)
    
    x,y=np.unique(apariciones, return_counts=True)
    
    R=y
    
    if i <=8:
        puntos=np.loadtxt('random.dat',usecols=0)
        total=len(puntos)
        nonzero=len(nodos)
        p0=np.array([(total-nonzero)])
        P=np.concatenate((p0,R),axis=0)
        P=P/total
        Srandom[i]=entropia(P)
        Crandom[i]=complejidad(P)
    else:
        puntos=np.loadtxt('norandom.dat',usecols=0)
        total=len(puntos)
        nonzero=len(nodos)
        p0=np.array([(total-nonzero)])
        P=np.concatenate((p0,R),axis=0)
        P=P/total
        Snorandom[i-9]=entropia(P)
        Cnorandom[i-9]=complejidad(P)
            
    
plt.figure(figsize=(5,4))
plt.plot(BETAS, Srandom, label='Random',marker='o')
plt.plot(BETAS, Snorandom, label='Simulation',marker='o')
plt.xlabel('\u03b2')
plt.ylabel('S')
plt.legend()
plt.savefig('SvBeta.pdf')

plt.figure(figsize=(5,4))
plt.plot(BETAS, Crandom, label='Random',marker='o')
plt.plot(BETAS, Cnorandom, label='Simulation',marker='o')
plt.xlabel('\u03b2')
plt.ylabel('C')
plt.legend()
plt.savefig('CvBeta.pdf')

plt.figure(figsize=(5,4))
plt.plot(Srandom, Crandom, label='Random',marker='o')
plt.plot(Snorandom, Cnorandom, label='Simulation',marker='o')
plt.xlabel('S')
plt.ylabel('C')
plt.legend()
plt.savefig('SvC.pdf')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    