#####################################################################################
### Program plotting and fitting Ribo-Seq data of k-somes for the gene HIST1H2BM. ###
#####################################################################################

# Running with python3
# 2 choices have to be made on lines 17 and 18.

# Author: Carole Chevalier

# Libraries
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
import scipy.special as spp
from math import*
import numpy as np
from scipy.optimize import curve_fit

### CHOICES ####
window=19 #smoothing window (19 in the article)
Fitdata=True #Show the ballistic densities with fitted parameters

#--------------------- FUNCTIONS --------------------
#Function calculating densities and smoothing data
def Smooth_density(data,windowf,n_ribo):
	for i in range(len(data)):
		data[i]["Lissage"] = data[i]["Coverage"].rolling(windowf,min_periods=1,center=True).mean()
		reads=data[i]["Lissage"] #reads is the column 1 of dataRibo
		norm=reads.agg('sum') #sum all reads
		data[i]["Lissage"]=data[i]["Lissage"]/norm*n_ribo #normalisation to k
		v=data[i]["Lissage"];

#Function plotting data
def Multiple_profiles(data,graphs,n_ligne,n_row,n_ribo) :
    ballistic_colors=['blue','orange','green','red']
    ax.set_xlabel(r"$i$ (codon)")
    if norm1==False:
            ax.set_ylabel(r"${\rho}_k(i)$")
    else:
            ax.set_ylabel(r"$\bar{\rho}_k(i)$")
    length=len(data[graphs[0]]['Lissage'])
    if window==1:
            ax.plot(np.linspace(1,length,length),data[graphs[0]]["Lissage"],'+-',label=str(n_ribo)+"-some")
    else:
            ax.scatter(np.linspace(int(window/2),length-int(window/2),length-2*int(window/2)+1),data[graphs[0]]["Lissage"][int(window/2):length-int(window/2)+1],marker='+',label=str(n_ribo)+"-some")
            if Fitdata==True:
                ax.plot(Ballistic_model_homo_dimensionless_parameters(0.06*mean_popt,mean_popt/3600,length,n_ribo)[0],color=ballistic_colors[n_ribo-1],lw=2)
    ax.set_xlim(0,length-1)
		
#k-some model function for curve_fit
def Ballistic_homo_for_fit(x,T_L):
	omega=1/3600
	L=length
	alpha=0.06
	p=L/T_L
	s=omega+alpha
	Pk=omega/s*(alpha/s)**k_*spp.gammainc(k_+1,s*T_L)+(alpha*T_L)**k_/factorial(k_)*exp(-s*T_L)
	ballistic_ksome_density=omega/Pk/p*(alpha/s)**k_*(spp.gammainc(k_,s*T_L)-spp.gammainc(k_,s*x/p))+alpha/Pk/p*(alpha*T_L)**(k_-1)/factorial(k_-1)*exp(-s*T_L)	
	return ballistic_ksome_density
    
#k-some model function for plotting
def Ballistic_model_homo_dimensionless_parameters(tilde_alpha,tilde_omega,L,k):
	s=tilde_alpha+tilde_omega
	density=[0]*L
	polysome_density=[0]*L
	if norm1==True: norm=k
	else: norm=1
	Pk=tilde_omega/s*(tilde_alpha/s)**k*spp.gammainc(k+1,s)+(tilde_alpha)**k/factorial(k)*exp(-s)
	for j in range(L):
		density[j]=(tilde_omega/Pk/L*(tilde_alpha/s)**k*(spp.gammainc(k,s)-spp.gammainc(k,s*j/L))+tilde_alpha/Pk/L*(tilde_alpha)**(k-1)/factorial(k-1)*exp(-s))/norm
		polysome_density[j]=tilde_alpha*exp(-tilde_omega*j/L)/L
	return density,polysome_density
#----------------------------------------------------------------------------------------------------

#-----------Data reading----------
filename = glob('Monosome/NM_003521*')
dataRibo1 = [pd.read_csv(f,delimiter="\t",usecols=[1,2]) for f in filename]
Smooth_density(dataRibo1,window,1)
print(dataRibo1)

filename = glob('Disome/NM_003521*')
dataRibo2 = [pd.read_csv(f,delimiter="\t",usecols=[1,2]) for f in filename]
Smooth_density(dataRibo2,window,2)
print(dataRibo2)

filename = glob('Trisome/NM_003521*')
dataRibo3 = [pd.read_csv(f,delimiter="\t",usecols=[1,2]) for f in filename]
Smooth_density(dataRibo3,window,3)

filename = glob('Quadrisome/NM_003521*')
dataRibo4 = [pd.read_csv(f,delimiter="\t",usecols=[1,2]) for f in filename]
Smooth_density(dataRibo4,window,4)

Data=[dataRibo1,dataRibo2,dataRibo3,dataRibo4]

# ----------------- FITTING -----------------    
norm1=False
length=len(dataRibo1[0]['Lissage'])
mean_popt=0
for k_ in range(1,5):
    popt, pcov = curve_fit(Ballistic_homo_for_fit,np.linspace(int(window/2),length-int(window/2),length-2*int(window/2)+1),Data[k_-1][0]['Lissage'][int(window/2):length-int(window/2)+1])
    print(str(k_)+'-some: T(L) = '+str(popt))
    mean_popt+=popt
mean_popt=mean_popt/4
print('Mean T(L) = '+str(mean_popt))

# ----------------- PLOTTING ---------------	
xsubplot=1
ysubplot=1

fig, ax=plt.subplots(ysubplot,xsubplot,figsize=(5,3),constrained_layout=True)
        
l = [0]
for k_ in range(4):
    Multiple_profiles(Data[k_],l,ysubplot,xsubplot,k_+1)

if Fitdata==True:
    ax.plot(Ballistic_model_homo_dimensionless_parameters(0.06*mean_popt,mean_popt/3600,length,1)[1], c='k',lw=2, label='polysome')
    
plt.legend()

fig, ax=plt.subplots(ysubplot,xsubplot,figsize=(5,3),constrained_layout=True)
norm1=True
for k_ in range(4):
    reads=Data[k_][0]["Lissage"] #reads is the column 1 of dataRibo
    norm=reads.agg('sum') #sum all reads
    Data[k_][0]["Lissage"]=Data[k_][0]["Lissage"]/norm #normalisation to 1
    Multiple_profiles(Data[k_],l,ysubplot,xsubplot,k_+1)
ax.set_ylim(0,0.04)
plt.legend()

plt.show()
