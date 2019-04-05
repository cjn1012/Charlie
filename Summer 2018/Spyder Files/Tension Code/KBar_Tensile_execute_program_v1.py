# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 15:52:07 2016

@author: smates
"""
from Tkinter import *
import numpy
import pylab
import KBar_Tensile_support_v2
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

analysis_file='Tension_Input_File.csv'
# This reads in the file containing experiment names and data analysis parameters for analyzing multiple experiments simultaneously

expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,sg3_hi,sg3_lo,gage_w,gage_thk,gage_len,spec_type,TC_type,N_sens,ibar_modulus,tbar_modulus,window,refl_delay,trans_delay,DIC,num_gages,tgage_type,spec_type,kistler,DIC_first,DIC_last,DIC_frame_rate,DIC_col=read_expts_for_analysis(analysis_file,1)

print(expt_list)

##############   Need to add 'T' to filename for tension file - only reading in test number into float

#compare_trans_forces(analysis_file,1)


#################################################
#
#   this calculates stress-strain curves for all experiments in the analysis file for tension
#
#################################################

#####these arrays are used to store final wave analysis delay parameters for determining effect of wave delay time variation
#num_expts=len(expt_list)
#print(num_expts)
#refl_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
#trans_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
#equilib=numpy.arange(0,num_expts,dtype=numpy.float)


incr=0
#plt.clf()    
for expt in expt_list:
#    strain,stress,strain_rate,refl_delay_final,trans_delay_final,criterion_final=calculate_stress_strain(str(int(expt_list[increment])),dia[increment],thk[increment],sg1_hi[increment],sg1_lo[increment],sg2_hi[increment],sg2_lo[increment],window[increment],refl_delay[increment],trans_delay[increment],gfoil[increment],bar_modulus[increment])          
    strain,stress,strain_rate=calculate_tensile_stress_strain(str(int(expt_list[incr])),gage_w[incr],gage_thk[incr],gage_len[incr],sg1_hi[incr],sg1_lo[incr],sg2_hi[incr],sg2_lo[incr],sg3_hi[incr],sg3_lo[incr],window[incr],refl_delay[incr],trans_delay[incr],DIC[incr],num_gages[incr],tgage_type[incr],spec_type[incr],kistler[incr],ibar_modulus[incr],tbar_modulus[incr],DIC_first[incr],DIC_last[incr],DIC_frame_rate[incr],DIC_col[incr])
    print(expt)
    
    # these two functions calculate and plot axial strain distributions averaged over cross-sections at specific avg strain levels
    DIC_time,DIC_strain,DIC_strain_std=calculate_dic_strain_time(str(int(expt_list[incr])),DIC_first[incr],DIC_last[incr],DIC_frame_rate[incr],gage_len[incr],DIC_col[incr])
    calculate_DIC_axial_strain_profiles(str(int(expt_list[incr])),DIC_first[incr],DIC_last[incr],DIC_frame_rate[incr],gage_len[incr],DIC_col[incr],DIC_strain)

    #plt.plot(strain,stress)
    incr+=1

#data_outfile = open("wave_diag.csv","w")
#data_outfile.write("expt, refl delay, trans delay, equilib criterion")
#data_outfile.write('\n')       
#data_to_write=numpy.column_stack((expt_list,refl_delay_f,trans_delay_f,equilib))
#numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
#data_outfile.close()  
##### show plots for all for stress-strain curves
#plt.show()  




#################################################
#
#   this calculates stress-strain curves for all experiments in the analysis file for compression
#
#################################################

#####these arrays are used to store final wave analysis delay parameters for determining effect of wave delay time variation
#num_expts=len(expt_list)
#refl_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
#trans_delay_f=numpy.arange(0,num_expts,dtype=numpy.float) 
#equilib=numpy.arange(0,num_expts,dtype=numpy.float)


#increment=0
#plt.clf()    
#for expt in expt_list:
#    strain,stress,strain_rate,refl_delay_final,trans_delay_final,criterion_final=calculate_stress_strain(str(int(expt_list[increment])),dia[increment],thk[increment],sg1_hi[increment],sg1_lo[increment],sg2_hi[increment],sg2_lo[increment],window[increment],refl_delay[increment],trans_delay[increment],gfoil[increment],bar_modulus[increment])          
    ##### these data are saved in an output file    
#    refl_delay_f[increment]=refl_delay_final
#    trans_delay_f[increment]=trans_delay_final
#    equilib[increment]=criterion_final
#    plt.plot(strain,stress)
#    increment+=1

#data_outfile = open("wave_diag.csv","w")
#data_outfile.write("expt, refl delay, trans delay, equilib criterion")
#data_outfile.write('\n')       
#data_to_write=numpy.column_stack((expt_list,refl_delay_f,trans_delay_f,equilib))
#numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
#data_outfile.close()  
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

#calculate_and_plot_temperatures(expt_list,TC_type,N_sens)


