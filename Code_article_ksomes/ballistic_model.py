#################################################################################
### Program generated ballistic density profiles for several values of R_1(0) ###
#################################################################################

# Author: Carole Chevalier

# Running with python3

# It is important to well determine the name of files for the next programs:
#   - For same parameters but several p_i sequences, change the letter after 'p_iSeq1' (use 'p_iSeq1a', 'p_iSeq1b', 'p_iSeq1c', ...)
#   - For other parameters you can use 'p_iSeq2' for example

# Libraries
import random as rand
from numba import jit
from numba.typed import List
import time
from math import*
from decimal import *
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib import rc
import numpy as np
import scipy.special as spp
import scipy.linalg
import scipy.integrate as integrate
from scipy.optimize import minimize
import sys
import statistics

#----------------------------------------Parameters-------------------------------
L=100 # number of codons
n_experience=1
F1sS1min=1 #R_1(0) min
F1sS1max=3 #R_1(0) max
seed=201034
rand.seed(seed)

dt=1
alpha=0.08*dt
t_half_life=log(2)/0.0001
p_min=0.2*dt
p_max=4*dt
dp=p_max-p_min
jump_r=List()

def ksomes_stochastic_death(L,k,p,alpha,t_half_life,T_cross_RNA):
    #Fk : first part of the density (mRNA in filling)
    #SK : second part of the density (mRNA in its own stationary state)
    density=[]
    Fk=[]
    Sk=[]
    T_i=0
    T_im1=0
    omega=log(2)/t_half_life
    s=omega+alpha
    P_k=omega/s*(alpha/s)**k*spp.gammainc(k+1,s*T_cross_RNA)+(alpha*T_cross_RNA)**k/factorial(k)*exp(-s*T_cross_RNA)
    if k>0:
        for i in range(L):
            T_i+=1/p[i]
            if i>0:
                    T_im1+=1/p[i-1]
            I=(T_i*spp.gammainc(k,s*T_i)-T_im1*spp.gammainc(k,s*T_im1))*factorial(k-1)+1/s*factorial(k)*(spp.gammainc(k+1,s*T_im1)-spp.gammainc(k+1,s*T_i))
            Fk.append(omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*T_cross_RNA)/p[i]-I/factorial(k-1)))
            Sk.append(alpha/p[i]/P_k*(alpha*T_cross_RNA)**(k-1)/factorial(k-1)*exp(-s*T_cross_RNA))
            density.append(omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*T_cross_RNA)/p[i]-I/factorial(k-1))+alpha/p[i]/P_k*(alpha*T_cross_RNA)**(k-1)/factorial(k-1)*exp(-s*T_cross_RNA))
    return density, P_k, Fk, Sk
    
def Density_polysomes_stochastic_death(L,alpha,p,t_half_life,T_cross_RNA):
	density=[]
	omega=log(2)/t_half_life
	T_i=0
	T_im1=0
	T_l=0
	#for i in range(l) :
		#T_l+=1/p[i]
	for i in range(L):
		T_i+=1/p[i]
		if i>0:
			T_im1+=1/p[i-1]
		density.append(alpha/omega*(exp(-omega*T_im1)-exp(-omega*T_i))) # J'avais mis ça :(1+alpha*T_l), pourquoi ? exclusion a l'entrée ?
		#density.append((exp(-omega*T_im1)-exp(-omega*T_i))/(1-exp(-omega*T_cross_RNA))) #normalised to 1
	return density

for i in range(L):
	jump_r.append(rand.random()*dp+p_min)
		
T_cross_RNA=0
for i in range(L):
	T_cross_RNA+=1/jump_r[i]
	
#plt.figure("model-curves",figsize=(7, 5))
for m in range(F1sS1min-1,F1sS1max): #boucle pour tracer plusieurs jeux de courbes avec paramètres différents	
        omega=(m+1)*alpha/(exp(alpha*T_cross_RNA)-1)
        #omega=1/2400
        t_half_life=log(2)/omega
        name='ballistic_p_iSeq1a_F1sS1_'+str(m+1)
        n_ksomes=[0]*(L+1)
        analytic_dens=np.zeros((5,L))
        analytic_prob=np.zeros((5,L))
        norm_tot=[0]*5
        count_tot=0
        mean_k_tot=0
        tab_reads=np.zeros((n_experience,L+1,L))
        tab_norm=np.zeros((n_experience,5))
        tab_0somes=[0]*n_experience
        count=[0]*n_experience
        mean_k=[0]*n_experience
        #------------------------ Calculate profiles ----------------------
        for k in range(1,5):
                analytic_dens[k]=ksomes_stochastic_death(L,k,jump_r,alpha,t_half_life,T_cross_RNA)[0]
                analytic_prob[k]=[i/k for i in analytic_dens[k]]
                #plt.plot(analytic_prob[k],label=r'analytic %s-somes, $\alpha=%s$' % (k,alpha))
        analytic_dens[0]=Density_polysomes_stochastic_death(L,alpha,jump_r,t_half_life,T_cross_RNA) 
        norm_poly=sum(analytic_dens[0])
        analytic_prob[0]=[j/norm_poly for j in analytic_dens[0]]
        #plt.plot(analytic_prob[0],color='k',label='analytic polysomes')

        #------------------------ Record in a file ------------------------
        fichier = open("Ribo-Seq_"+name+".dat", "w")
        fichier.write("#codon\tjump rates\tpolysome\t1some\t2some\t3some\t4some\n")
        for i in range(L):
                fichier.write("%d\t%g\t%g\t%g\t%g\t%g\t%g\n" % (i,jump_r[i]/dt,analytic_dens[0][i],analytic_dens[1][i],analytic_dens[2][i],analytic_dens[3][i],analytic_dens[4][i]))
        fichier.write("#alpha = "+str(alpha/dt)+"\n#p_min = "+str(p_min)+", p_max = "+str(p_max)+", T(L) = "+str(T_cross_RNA)+"\n#half-life time = "+str(t_half_life)+"\n<k> = "+str(alpha/dt*t_half_life*dt/log(2)*(1-exp(-log(2)/t_half_life/dt*T_cross_RNA*dt)))+"\n#generator seed = "+str(seed))
        fichier.close()
plt.show()
