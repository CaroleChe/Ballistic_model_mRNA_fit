#################################################################################################################
### Program plotting the score function using several fit test of a same mRNA but different p_i distributions ###
#################################################################################################################
# Score function: errors made on the fitted parameters (alpha, p_i) in function of the value of (F_1(0)/S_1=R_1(0))

# Author: Carole Chevalier

# Running with python3
# Name of the first sequence to be read at line 41. The other names are generated automatically if they have the proper format but the constant part of the names have to be adapted manually (l. 43).

# Libraries:
import matplotlib.pyplot as plt
import numpy as np
from math import*
from decimal import *
import matplotlib.ticker as mticker
import re
import pandas as pd
import os
import string
#os.nice(20)

# MAKE THE CHOICES !
n_seq=2 # Number of profiles with same parameters except p_i distribution
F1sS1_min=1 #R1(0) min
F1sS1_max=3 #R1(0) max
L=100  # number of codons
n_dot=F1sS1_max-F1sS1_min+1
Mean_over_sequences=True
Folder='' # Folder where the data of the fit are
getcontext().prec = 2
alphaR=0.08 # Alpha which have to be retrieve 

# Extracting p_i that have to be found
RiboSeq_data=[0]*n_seq
filename=[]
minimize_filenames=[]
real_p=[0]*n_seq 
real_TL=[0]*n_seq

filename.append(Folder+"Ribo-Seq_ballistic_p_iSeq1a_F1sS1_"+str(F1sS1_min)+".dat")
for i in range(1,n_seq):
	filename.append(Folder+"Ribo-Seq_ballistic_p_iSeq1"+list(string.ascii_lowercase)[i]+"_F1sS1_"+str(F1sS1_min)+".dat")

for i in range(n_seq):
	RiboSeq_data[i]=pd.read_csv(filename[i],delimiter="\t",usecols=[1,2],engine='python',skipfooter=5)
	real_p[i]=RiboSeq_data[i]["jump rates"]
	real_TL[i]=sum(1/real_p[i])
#print(real_p)	

# Function calculating the corresponding omega axis
def Omega(X):
	V=[0]*len(X)
	for n in range(len(X)):
		for i in range(n_seq):
			V[n]+=X[n]*alphaR/(exp(alphaR*real_TL[i])-1)
	V=[z/n_seq*100 for z in V]
	return ["%.2f" % z for z in V]

def F1sS1_exacte(omega):
	V=[]
	for i in range(len(omega)):
		V.append(round(omega[i]/(omega[i]+alphaR)*(exp((omega[i]+alphaR)*77.85)-1)))
	return V

def relative_error(true_value,fit_value):
	return abs(true_value-fit_value)/true_value

def ReadData(seq_names):
	filename=Folder+seq_names
	fichier=open(filename,'r')#?
	pattern = re.compile("x: array")
	startE=0
	for m in re.finditer(pattern, fichier.read()):
		startE=m.end()+2
	fichier.seek(startE)
	pattern = re.compile("]")
	for m in re.finditer(pattern, fichier.read()):
		stopE=m.start()
	fichier.seek(startE)
	fit_results=fichier.readlines(stopE-30)
	data = []
	for line in fit_results:
		f_list = [float(i) for i in line.split(",") if i.strip()]
		data+=f_list[0:len(f_list)]
	line = fichier.readline()
	data+=[float(i) for i in line.split("])") if i.strip()]
	fichier.close()
	alpha=data[0]
	p_i=data[1:]
	return alpha,p_i

def pi_from_Ti(Seq):
	found_p=[1/Seq[1][0]]
	for i in range(0,len(Seq[1])-1):	
		found_p.append(1/(Seq[1][i+1]-Seq[1][i]))
	return found_p

#Function calculating R1(0)=F_1(0)/S_1
def R1_approx(alphaTL,omegaTL):
	F1sS1=omegaTL/alphaTL*(exp(alphaTL)-1)
	return F1sS1	

# Errors and omega calculation
MeanErr_alpha=[0]*n_dot
MeanErr_mean_p=[0]*n_dot
MeanErr_max_p=[0]*n_dot
found_p=[0]*n_seq
result_seq=[0]*n_seq
SeqMean_p_error=np.zeros((n_seq,L))
Err_alpha=np.zeros((n_seq,n_dot))
Err_mean_p=np.zeros((n_seq,n_dot))
Err_max_p=np.zeros((n_seq,n_dot))

for n in range(F1sS1_min,F1sS1_max+1):
        for i in range(n_seq):
                #print(n)
                result_seq[i]=ReadData('minimize_R-S_ballistic_p_iSeq1'+list(string.ascii_lowercase)[i]+'_F1sS1_'+str(n)+'.dat')
                found_p[i]=pi_from_Ti(result_seq[i])
                Err_alpha[i][n-F1sS1_min]=relative_error(alphaR,result_seq[i][0])
                MeanErr_alpha[n-F1sS1_min]+=Err_alpha[i][n-F1sS1_min]

                for j in range(L):
                        #print('real ',real_p[i][j])
                        #print('found ',found_p[i][j])
                        SeqMean_p_error[i][j]=relative_error(real_p[i][j],found_p[i][j])
                        
                        #print('err ',SeqMean_p_error[i][j])

                Err_mean_p[i][n-F1sS1_min]=sum(SeqMean_p_error[i])/L;Err_max_p[i][n-F1sS1_min]=max(SeqMean_p_error[i])
        
                MeanErr_mean_p[n-F1sS1_min]+=Err_mean_p[i][n-F1sS1_min]

                MeanErr_max_p[n-F1sS1_min]+=Err_max_p[i][n-F1sS1_min]
                
MeanErr_alpha=[i/n_seq for i in MeanErr_alpha]
MeanErr_mean_p=[i/n_seq for i in MeanErr_mean_p]
MeanErr_max_p=[i/n_seq for i in MeanErr_max_p]

F1sS1=np.arange(F1sS1_min,F1sS1_max+1,1)

plt.figure(figsize=[5.8,3.2],constrained_layout=True)
plt.axes(axisbelow=True) # Pour que la grille soit en arrière plan
for i in range(n_seq):
	plt.scatter(F1sS1,Err_max_p[i],s=80,marker='v',color='palegreen')#,color='palegreen'
	plt.scatter(F1sS1,Err_alpha[i],s=110,marker='*',color='peachpuff')#,color='peachpuff'
	plt.plot(F1sS1,Err_mean_p[i],'--o',color='lightblue')#,color='lightblue'
		
if Mean_over_sequences==True :
	plt.scatter(F1sS1,MeanErr_max_p,s=80,marker='v',color='limegreen',label='max $p_i$ error',zorder=8)
	plt.scatter(F1sS1,MeanErr_alpha,s=110,marker='*',color='r',label=r'$\alpha$ error',zorder=9)
	plt.plot(F1sS1,MeanErr_mean_p,'--o',label='mean $p_i$ error',zorder=10)
plt.legend()

plt.axhline(y = 0.1, color = 'k', linestyle = '--',lw=1)

F1sS1_axe=[]
for i in range(0,len(F1sS1),2):
	F1sS1_axe.append(F1sS1[i])
plt.xticks(F1sS1_axe)
plt.xlim(1/2,15.5)
extraticks=[0.1]
plt.yticks(list(plt.yticks()[0]) + extraticks)
#plt.ylim(-0.05,0.9)
plt.xlabel('$\mathcal{R}_1(0)$')
plt.ylabel('relative error')
plt.grid()

ax=plt.twiny()
plt.xlim(1/2,15.5)
ax.set_xticks(F1sS1_axe)
ax.set_xticklabels(Omega(F1sS1_axe))
ax.set_xlabel(r'mean $\omega\times 100$')
#ax.text(16,0.93,r'$\times 10^{-2}$')

#plt.yscale("log")
	
plt.savefig('score_fct_1somesPolysomes_test.pdf',format='pdf')
	



plt.show()
