# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 15:52:07 2016

@author: smates
"""
#from Tkinter import *
import numpy
import pylab
import KBar_data_inspector_support as ksb
import matplotlib.pyplot as plt

global average_emissivity
#global equilib


#expt='1824' 
#quickplot_volts_data(expt)
 

##############################
# this identifies wave timing from tests without samples
###############################
#bars_only_test='1824'
#window=400
#sg1_hi=60.0
#sg1_lo=0.0
#sg2_hi=-60.0
#sg2_lo=0.0
#calculate_wave_timing(bars_only_test,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window)


##########################################################################################################
#
# This loop calculates stres-strain data for the input list of experiments subject to the input data parameters from 
#
##########################################################################################################

#analysis_file='1018HT_inputfile.csv'
# This reads in the file containing experiment names and data analysis parameters for analyzing multiple experiments simultaneously
#compare_trans_forces(analysis_file)



#################################################
#
#   this calculates stress-strain curves for all experiments in the analysis file
#
#################################################

analysis_file='pure_iron_analysis_file.csv'
expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay=read_expts_for_analysis(analysis_file)
#SOMETHING GOES HERE, ASK BENNY

#####these arrays are used to store final wave analysis delay parameters for determining effect of wave delay time variation
num_expts=len(expt_list)
refl_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
trans_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
equilib=numpy.arange(0,num_expts,dtype=numpy.float)


increment=0
#plt.clf()    
for expt in expt_list:
    print('experiment:',expt)    
    strain,stress,strain_rate,refl_delay_final,trans_delay_final,criterion_final=calculate_stress_strain(str(int(expt_list[increment])),dia[increment],thk[increment],sg1_hi[increment],sg1_lo[increment],sg2_hi[increment],sg2_lo[increment],window[increment],refl_delay[increment],trans_delay[increment],gfoil[increment],bar_modulus[increment])          
    ##### these data are saved in an output file    
    refl_delay_f[increment]=refl_delay_final
    trans_delay_f[increment]=trans_delay_final
    equilib[increment]=criterion_final
    #plt.plot(strain,stress)
    increment+=1

data_outfile = open("wave_diag.csv","w")
data_outfile.write("expt, refl delay, trans delay, equilib criterion")
data_outfile.write('\n')       
data_to_write=numpy.column_stack((expt_list,refl_delay_f,trans_delay_f,equilib))
numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
data_outfile.close()  
##### show plots for all for stress-strain curves
#plt.show()  

##########################################################################################################
#
# For averaging stress-strain data at a single temperature
#
##########################################################################################################

##### specify strain range to average
#low_strain=0.05
#hi_strain=0.30
#average_stress_strain_data(expt_list,low_strain,hi_strain)


##########################################################################################################
#
# For thermal plotting analysis and plotting
#
##########################################################################################################
average_emissivity=0.25
calculate_and_plot_temperatures(expt_list,TC_type,N_sens)


