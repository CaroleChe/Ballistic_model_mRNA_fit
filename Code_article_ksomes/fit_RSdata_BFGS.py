##################################################################################
### Program doing the fit method with 1some and polysome to find alpha and p_i ###
##################################################################################

# Author: Carole Chevalier

# Running with python3
# Fit only 1 mRNA at a time
# Put the name of the file that have to be read on "filename" l. 37. 
# Keep the name format of output files to plot the fit results with score_fct.py: 
#   - For same parameters but several p_i sequences, change the letter after 'p_iSeq1' (use 'p_iSeq1a', 'p_iSeq1b', 'p_iSeq1c', ...)
#   - For other parameters you can use 'p_iSeq2' for example

# Libraries:
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from scipy.optimize import minimize
import pandas as pd
from math import*
import scipy.special as spp
import time
import random as rand
import numba as nb
import statistics as st

start_time = time.time()

# MAKE THE CHOICES !
F1sS1=3 #R1(0)
Folder="" # To put the outputs in another folder than this script
FigureName="R-S_ballistic_p_iSeq1a_F1sS1_"+str(F1sS1) # Name of outputs
alpha=0.08 # Alpha which have to be retrieve 
drawfigure=False # show p_i errors, densities, ...

# DATA READING
filename="Ribo-Seq_ballistic_p_iSeq1a_F1sS1_"+str(F1sS1)+".dat" # k-somes densties normalized to 1 in this file
RiboSeq_data=pd.read_csv(filename,delimiter="\t",usecols=[1,2,3,4,5,6],engine='python',skipfooter=5)

# FUNCTION TO MINIMIZE 
def Density_ksomes_stochastic_death_for_minimize(parameters_to_fit,*arg):
    omega=arg[1]; L=arg[2]; # density_data : argument 0
    alpha=parameters_to_fit[0];
    density=np.zeros((4,L)); gap=0.0	
    s=omega+alpha
    
    # minimization of monosome
    k=1
    P_k=omega/s*(alpha/s)**k*spp.gammainc(k+1,s*parameters_to_fit[L])+(alpha*parameters_to_fit[L])**k/factorial(k)*exp(-s*parameters_to_fit[L])
    P_km1=omega/s*(alpha/s)**(k-1)*spp.gammainc(k,s*parameters_to_fit[L])+(alpha*parameters_to_fit[L])**(k-1)/factorial(k-1)*exp(-s*parameters_to_fit[L])
    
    I=(parameters_to_fit[1]*spp.gammainc(k,s*parameters_to_fit[1]))*factorial(k-1)+1/s*factorial(k)*(-spp.gammainc(k+1,s*parameters_to_fit[1]))
    density[k-1][0]=omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*parameters_to_fit[L])*parameters_to_fit[1]-I/factorial(k-1))+alpha*parameters_to_fit[1]/P_k*(alpha*parameters_to_fit[L])**(k-1)/factorial(k-1)*exp(-s*parameters_to_fit[L])
    gap+=(abs(density[k-1][0]-arg[0][k-1][0])/arg[0][k-1][0])**2
    for j in range(1,L):
            p_j=1/(parameters_to_fit[j+1]-parameters_to_fit[j])
            I=(parameters_to_fit[j+1]*spp.gammainc(k,s*parameters_to_fit[j+1])-parameters_to_fit[j]*spp.gammainc(k,s*parameters_to_fit[j]))*factorial(k-1)+1/s*factorial(k)*(spp.gammainc(k+1,s*parameters_to_fit[j])-spp.gammainc(k+1,s*parameters_to_fit[j+1]))
            density[k-1][j]=omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*parameters_to_fit[L])/p_j-I/factorial(k-1))+alpha/p_j/P_k*(alpha*parameters_to_fit[L])**(k-1)/factorial(k-1)*exp(-s*parameters_to_fit[L])
            gap+=(abs(density[k-1][j]-arg[0][k-1][j])/arg[0][k-1][j])**2
                            
    #minimization of polysome
    dens_poly=[]
    dens_poly.append((1-exp(-omega*parameters_to_fit[1]))/(1-exp(-omega*parameters_to_fit[L])))
    gap+=(abs(dens_poly[0]-arg[0][4][0])/arg[0][4][0])**2
    for j in range(1,L):
            dens_poly.append((exp(-omega*parameters_to_fit[j])-exp(-omega*parameters_to_fit[j+1]))/(1-exp(-omega*parameters_to_fit[L])))
            gap+=(abs(dens_poly[j]-arg[0][4][j])/arg[0][4][j])**2
                    
    return gap

# FUNCTION TO PLOT K-SOME DENSITIES OF THE MODEL
def Density_ksomes_stochastic_death(fitted_parameters, *arg):
	omega=arg[1]; L=arg[2]; # density_data : argument 0
	alpha=fitted_parameters[0];

	density=np.zeros((5,L)); 
	
	s=omega+alpha
	
	for k in range(1,5):
		P_k=omega/s*(alpha/s)**k*spp.gammainc(k+1,s*fitted_parameters[L])+(alpha*fitted_parameters[L])**k/factorial(k)*exp(-s*fitted_parameters[L])
		
		I=(fitted_parameters[1]*spp.gammainc(k,s*fitted_parameters[1]))*factorial(k-1)+1/s*factorial(k)*(-spp.gammainc(k+1,s*fitted_parameters[1]))
		density[k-1][0]=omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*fitted_parameters[L])*fitted_parameters[1]-I/factorial(k-1))+alpha*fitted_parameters[1]/P_k*(alpha*fitted_parameters[L])**(k-1)/factorial(k-1)*exp(-s*fitted_parameters[L])
		for j in range(1,L):
			p_j=1/(fitted_parameters[j+1]-fitted_parameters[j])
			I=(fitted_parameters[j+1]*spp.gammainc(k,s*fitted_parameters[j+1])-fitted_parameters[j]*spp.gammainc(k,s*fitted_parameters[j]))*factorial(k-1)+1/s*factorial(k)*(spp.gammainc(k+1,s*fitted_parameters[j])-spp.gammainc(k+1,s*fitted_parameters[j+1]))
			density[k-1][j]=omega/P_k*(alpha/s)**k*(spp.gammainc(k,s*fitted_parameters[L])/p_j-I/factorial(k-1))+alpha/p_j/P_k*(alpha*fitted_parameters[L])**(k-1)/factorial(k-1)*exp(-s*fitted_parameters[L])
	
	density[4][0]=(1-exp(-omega*fitted_parameters[1]))/(1-exp(-omega*fitted_parameters[L]))
	for j in range(1,L):
		density[4][j]=(exp(-omega*fitted_parameters[j])-exp(-omega*fitted_parameters[j+1]))/(1-exp(-omega*fitted_parameters[L]))
			
	return density
    

#print(RiboSeq_data)

# DECLARATIONS
name="some"
L=np.size(RiboSeq_data[str(1)+name])
xdata = np.linspace(1,L,L)

# Calcul of T_i
Ti_reels=[1/RiboSeq_data["jump rates"][0]]
for i in range(L-1):
	Ti_reels.append(Ti_reels[i]+1/RiboSeq_data["jump rates"][i+1])
print("T(L) = ",Ti_reels[L-1])

#Function calculating R1(0)=F_1(0)/S_1
def R1_approx(alphaTL,omegaTL):
	F1sS1=omegaTL/alphaTL*(exp(alphaTL)-1)
	return F1sS1

#------------------------ MINIMISATION-----------------------------
fichier_minimize = open(Folder+"minimize_"+FigureName+".dat", "w")
fichier_minimize.write("alpha_0\tT(L)_0\tchi2\talpha_fit\tT(L)_fit\n")
# ARGUMENTS LIST FOR MINIMISATION 1
norm_poly=sum(RiboSeq_data["poly"+name][:])
dens_RiboSeq=np.array([RiboSeq_data["1"+name][:],RiboSeq_data["2"+name][:],RiboSeq_data["3"+name][:],RiboSeq_data["4"+name][:],[i/norm_poly for i in RiboSeq_data["poly"+name][:]]])
args_list=(dens_RiboSeq, F1sS1*alpha/(exp(alpha*Ti_reels[L-1])-1), L)

for es_a in range(1):	
        for es_T in range(1):	
                # INITIAL GUESS and BOUNDS
                initial_parameters=np.array(0.03) # alpha
                bounds_array=[(0.01,0.2)] # alpha
                for i in range(L):
                        initial_parameters=np.append(initial_parameters,(i+1))# T(i)
                        bounds_array.append(((i+1)/10,(i+1)*10))
                #initial_parameters=np.append(initial_parameters,L/2)
                #bounds_array.append((50,200))	
                path_ = [initial_parameters]

                res=minimize(Density_ksomes_stochastic_death_for_minimize,initial_parameters,args_list,method='L-BFGS-B',bounds=bounds_array,options={'ftol': 1e-15, 'gtol': 1e-10,"maxiter":10**1,'maxfun':10**9,'disp': False})
                fichier_minimize.write("%g\t%g\t%g\t%g\t%g\n" % (initial_parameters[0],initial_parameters[L],res.fun,res.x[0],res.x[L]))
                
                print(res)
                print("Simu (%s,%s) temps execution = %s s" % (es_a,es_T,time.time() - start_time))
                fichier_minimize.write("\nSimu (%s,%s) temps execution = %s s\n\n" % (es_a,es_T,time.time() - start_time))
                fichier_minimize.write(str(res))

# ----------------------- DRAW FIGURES ----------------------------------
if drawfigure==True :	

        alpha_min, alpha_max, alpha_step = 0.03, 0.1, .002
        TL_min, TL_max, TL_step = bounds_array[L][0], bounds_array[L][1], 1

        x, y = np.meshgrid(np.arange(alpha_min, alpha_max+alpha_step, alpha_step), np.arange(TL_min, TL_max+TL_step, TL_step))
        #print([x,res.x[1:L], y])

        fixed_parameters=[]
        for i in range(1,L):
                fixed_parameters.append(res.x[i])
        #print(fixed_parameters)
                
        densities=[]; 
        for i in range(len(x)):
                densitiesX=[]
                for j in range(len(x[0])):
                        densitiesX.append(Density_ksomes_stochastic_death_for_minimize([x[i][j]]+fixed_parameters+[y[i][j]],*args_list))		
                densities.append(densitiesX)

        # PLOT THE DENSITIES
        plt.figure("ksomes_densities")
        gen_coord=np.linspace(1,L,L)
        for k in range(1):
                plt.plot(gen_coord,[i*(k+1) for i in RiboSeq_data[str(k+1)+name]],label=str(k+1)+'some')
                plt.plot(gen_coord,[i for i in Density_ksomes_stochastic_death(res.x, *args_list)[k]], "--",label=str(k+1)+'some')
        plt.legend()
        plt.savefig(Folder+"ksomes_densities_"+FigureName+".png")
        plt.figure("polysome_density")
        plt.plot(gen_coord,[i/sum(RiboSeq_data['poly'+name]) for i in RiboSeq_data['poly'+name]])
        plt.plot(gen_coord,Density_ksomes_stochastic_death(res.x, *args_list)[4],"--")
        plt.xlabel('codon')
        plt.ylabel('polysome normalized density')
        plt.savefig(Folder+"polysome_density_"+FigureName+".png")
                
        # PLOT JUMP RATES
        plt.figure("T(i)")
        plt.plot(gen_coord,Ti_reels)
        plt.plot(gen_coord,res.x[1:], "--")
        plt.savefig(Folder+"Ti_"+FigureName+".png")
        
        # PLOT ERRORS	
        rate_p=[]
        error=[]
        rate_p.append(1/res.x[1])
        error.append(abs(RiboSeq_data["jump rates"][0]-rate_p[0])/RiboSeq_data["jump rates"][0])
        for i in range(1,L):
                rate_p.append(1/(res.x[i+1]-res.x[i]))
                error.append(abs(RiboSeq_data["jump rates"][i]-rate_p[i])/RiboSeq_data["jump rates"][i])
                
        plt.figure("jump rates")
        gen_coord=np.linspace(1,L,L)
        plt.plot(gen_coord,RiboSeq_data["jump rates"],'o--',label='to find')
        plt.plot(gen_coord,rate_p,'o--',label='found')
        plt.xlabel('i')
        plt.ylabel('$p_i$')
        plt.legend()
        plt.savefig(Folder+"p_i_"+FigureName+".png")

        plt.figure("error")
        plt.plot(gen_coord,error,'o--')
        plt.xlabel('i')
        plt.ylabel('$p_i$ relative error')
        plt.savefig(Folder+"p_i_errors_"+FigureName+".png")
        
        fichier_minimize.write("\n RELATIVE ERRORS \n Mean p_i error: "+str(st.mean(error))+"\n max p_i error: "+str(max(error)))

fichier_minimize.close()
print("temps execution = %s s" % (time.time() - start_time))

#plt.show()

