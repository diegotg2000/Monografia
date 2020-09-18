import numpy as np
from scipy.stats import entropy

def entropia(A):
    return entropy(A,base=2)

def D(A,B):
    return entropia((A+B)/2)-entropia(A)/2-entropia(B)/2         

def H(A):
    return entropia(A)/np.log2(len(A))

def Dstar(n):
    return -0.5*((n+1)/n*np.log2(n+1)+np.log2(n)-2*np.log2(2*n))

def complejidad(A,B):
    n=len(A)
    return D(A,B)*H(A)/Dstar(n)


BETAS=np.arange(1.0,5.1,0.5)
np.savetxt('Betas_varianza.txt',BETAS)
numeros=[]
for i in range(20):
    if i<10:
        numeros.append('0'+str(i))
    else:
        numeros.append(str(i))


for j in range(20):
    Com=np.zeros(len(BETAS))
    for i in range(len(BETAS)):
        #se lee los datos random
        A=np.loadtxt('random_sphere'+numeros[j]+'_'+str(BETAS[i])+'.BSKIndex', usecols=0)
        B=np.loadtxt('random_sphere'+numeros[j]+'_'+str(BETAS[i])+'.BSKIndex', usecols=1)
        C=np.concatenate((A,B),axis=0)
    
        nodos_ran, conexiones_ran = np.unique(C, return_counts=True)
    
        #se lee los datos simulados
        A=np.loadtxt('sphere'+numeros[j]+'_'+str(BETAS[i])+'.BSKIndex', usecols=0)
        B=np.loadtxt('sphere'+numeros[j]+'_'+str(BETAS[i])+'.BSKIndex', usecols=1)
        C=np.concatenate((A,B),axis=0)
    
        nodos_dat, conexiones_dat = np.unique(C, return_counts=True)
    
        #se hace el histograma del numero de conexiones
    
        max_con=max([max(conexiones_dat),max(conexiones_ran)])
    
        y_ran,x_ran=np.histogram(conexiones_ran, bins=max_con, range=(1,max_con+1))
    
        y_dat,x_dat=np.histogram(conexiones_dat, bins=max_con, range=(1,max_con+1))
    
        #se analiza random
        puntos_ran=np.loadtxt('random_sphere'+numeros[j]+'.dat',usecols=0)
        total_ran=len(puntos_ran)
        nonzero_ran=len(nodos_ran)
        p0_ran=np.array([(total_ran-nonzero_ran)])
        P=np.concatenate((p0_ran,y_ran),axis=0)
        y_ran=P/total_ran
        #se analiza datos
        puntos_dat=np.loadtxt('sphere'+numeros[j]+'.dat',usecols=0)
        total_dat=len(puntos_dat)
        nonzero_dat=len(nodos_dat)
        p0_dat=np.array([(total_dat-nonzero_dat)])
        P=np.concatenate((p0_dat,y_dat),axis=0)
        y_dat=P/total_dat

    
        if len(x_dat)==len(x_ran):   
            Com[i]=complejidad(y_dat,y_ran)
        else:
            print('ojo!')
            
    np.savetxt('CvB_{}.txt'.format(j+1), Com)
    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    