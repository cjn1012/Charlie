# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 15:52:07 2016

@author: smates
"""
from Tkinter import *
import numpy
import pylab

global average_emissivity

#pylab.clf()  # this is needed to clear data from previous plots

#expt='3509' 
#preview_KBar_data_file('3424.txt')

# print('what the fuck is going on?')


#quickplot_volts_data('3325')
#quick_stress_strain('3327')

########################################################################################
#
#   Function to to quickly calculate a stress strain curve to see if the data are good
#
########################################################################################
def quick_stress_strain(expt):
    # test block to look at stress-strain calcualation
    #expt='3530'
 #   sg1_col,sg2_col,TC_col,N_col,N2_col=find_data_columns(expt)

    print(expt)
 
    dia=0.004
    thk=0.002
    #sg1_col=1
    #sg2_col=2
    sg1_hi=60.0
    sg1_lo=0.0
    sg2_hi=60.0
    sg2_lo=0.0
    window=400
   # inc_start,refl_start,trans_start=calculate_wave_timing(expt,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window)
   # print(inc_start,refl_start,trans_start)
    r_delay=0.000302
    t_delay=0.000302
    mod=190.0
    gfoil=1
    
    #print('i am doing nothing in this function')
    strain,stress,strain_rate,opt_refl_delay,opt_trans_delay,equilib=calculate_stress_strain(expt,dia,thk,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window,r_delay,t_delay,gfoil,mod)
     
    quickplot_xydata(strain,stress)
    return()

########################################################################################
#
#   Function to compare trans wave signals (in terms of force) to examine friction effects
#
########################################################################################
def compare_trans_forces(analysis_file):
    # test block to look at stress-strain calcualation
    #expt='3530'
 #   sg1_col,sg2_col,TC_col,N_col,N2_col=find_data_columns(expt)

    # open sg2 data for experiment
    #  convert to force and zero (cal volts, bar area and modulus)
    #  find trans pulse
#      write trans pulse to test#_force.csv file
#       average force.csv files and plot

#use same analysis input file as used for other items
    expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay=read_expts_for_analysis(analysis_file)

    increment=0
    plt.clf()    
    for expt in expt_list:
        test=str(int(expt)) 
        testfile=test + '.txt'
        output_file=test+'force_time.txt'
        sg1_col,sg2_col,dummy1,dummy2,dummy3=find_data_columns(test) #find sg2 column
        sg2_time,sg2_volts=load_timeseries_data(testfile,8,0,sg2_col)
        dt=sg2_time[1]-sg2_time[0]        
        sg2_volts=zero_timeseries_data(sg2_volts,100)# 100 is the number of data points to average  
        sg2_cal_volts=2*(sg2_hi[increment]-sg2_lo[increment])/1000.0
        sg2_cal=0.004854/sg2_cal_volts    
        sg2_strain=volts_to_strain_data(sg2_volts,sg2_cal)
        start_time,time_window,sg2_strain_window=find_data_pulse4(sg2_time,sg2_strain,window[increment],3)

        trans_max = numpy.amax(sg2_strain_window)    
        trans_min = numpy.amin(sg2_strain_window)
        if abs(trans_max) < abs(trans_min):
            sg2_strain_window=-1.0*sg2_strain_window



        sg2_force=sg2_strain_window*bar_modulus[increment]*1.0e9*3.14159/4*0.015   # assume 15 mm bar in computing force
        time_for_plot=time_window-start_time


        plt.plot(time_for_plot,sg2_force)
        write_text_file(output_file,time_for_plot,sg2_force,sg2_force)

       #print('i am here # ',increment)   
        increment+=1
    
   # data_outfile = open("wave_diag.csv","w")
   # data_outfile.write("expt, refl delay, trans delay, equilib criterion")
   # data_outfile.write('\n')       
   # data_to_write=numpy.column_stack((expt_list,refl_delay_f,trans_delay_f,equilib))
   # numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
   # data_outfile.close()  
    # show plots for all for stress-strain curves
    plt.show()  

    # call function to compute average force-time curve

    file_list=expt_list
    file_root='force_time.txt'
    low_x=0.0
    hi_x=0.0001
    delta_x=0.000001
    average_xy_data(file_list,file_root,low_x,hi_x,delta_x)


    return()



########################################################################################
#
#   Function to read a text file containing experiment files and analysis parameters from a single text file
#
########################################################################################

def read_expts_for_analysis(analysis_file,type=0):

    #type 0 is compression
    #type 1 is tension

    infile = open(analysis_file,"r")
    
    if type==0:
        expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay=numpy.loadtxt(infile,skiprows=1,delimiter=',',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)
        return(expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay)
    
    if type==1:
        expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,sg3_hi,sg3_lo,gage_w,gage_thk,gage_l,spec_type,TC_type,N_sens,ibar_modulus,tbar_modulus,window,refl_delay,trans_delay,DIC,num_gages,tgage_type,spec_type,kistler,DIC_first,DIC_last,DIC_frame_rate,DIC_col=numpy.loadtxt(infile,skiprows=1,delimiter=',',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),unpack=True)
        return(expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,sg3_hi,sg3_lo,gage_w,gage_thk,gage_l,spec_type,TC_type,N_sens,ibar_modulus,tbar_modulus,window,refl_delay,trans_delay,DIC,num_gages,tgage_type,spec_type,kistler,DIC_first,DIC_last,DIC_frame_rate,DIC_col)
    else:
        print('error on input file type')        

###########################################################################################################
#
#   this function is called to average sets of xy data between a low and high x value 
#
############################################################################################################


def average_xy_data(file_list,file_root,low_x,hi_x,delta_x=1):
    import numpy as np
    import math
    #import pickle
    from scipy import interpolate
    import pylab
    import matplotlib.pyplot as plt


    num = np.int_((hi_x-low_x)/delta_x)
    newx=np.linspace(low_x,hi_x,num=num)
    
    # this is the array used to compute the average of the stress-strain curves and statistics
    avg_y=np.array(np.linspace(low_x,hi_x,num=num))
    avg_y_stdev=np.array(np.linspace(low_x,hi_x,num=num))
    avg_y_stdom=np.array(np.linspace(low_x,hi_x,num=num))


    #this loop obtains stress, strain and temperature data for fitting - all the data goes into a single array
    index=0   
    for file_number in file_list:
        index = index + 1
        file_number=str(file_number)
        print(file_number[0:4])
        datafile=file_number[0:4] + file_root
        data_input = open(datafile,"r")
        x_data,y_data,dummy_data=np.loadtxt(data_input,skiprows=0,delimiter=',',usecols=(0,1,2),unpack=True)
        pylab.plot(x_data,y_data)
    # interpolate current data set at strains and accumulate into 
        interp=interpolate.interp1d(x_data,y_data)
        interp_y_data = interp(newx) 
        #plt.plot(strains,interp_stress_data,color=color,label=label)
        y_sum = np.array(interp_y_data)
        avg_y=avg_y.__add__(y_sum)
    # compute average stress at each interpolated strain value and plot average result (maybe with a thicker line?)   
    avg_y=avg_y.__truediv__(index)
    pylab.plot(newx,avg_y,'r--',linewidth=10.0,label="avg") 
    pylab.legend(loc='lower right')
    pylab.show()
    #pylab.close()

    # loop again through the data and assemble standard deviation and stdom's
    for file_number in file_list:
        file_number=str(file_number)
        print(file_number[0:4])
        datafile=file_number[0:4] + file_root
        data_input = open(datafile,"r")
        x_data,y_data,dummy_data=np.loadtxt(data_input,skiprows=0,delimiter=',',usecols=(0,1,2),unpack=True)
        interp=interpolate.interp1d(x_data,y_data)
        interp_y_data = interp(newx) 
        error=interp_y_data.__sub__(avg_y)
        error=error.__pow__(2)
        avg_y_stdev = avg_y_stdev.__add__(error)
        #print('I am calculating averages')        
        #rint(avg_stress[100])        
        #print(avg_stress_stdev[100])
    
    avg_y_stdev=avg_y_stdev.__truediv__(index-1)
    avg_y_stdev=avg_y_stdev.__pow__(0.5)
    root_index=math.sqrt(index)
    avg_y_stdom=avg_y_stdev.__truediv__(root_index)
    
   


    # write output into .csv format for importing into excel or other file
    
    output = open("avg_y.csv","w")
    output.write("root index: ")
    output.write(str(root_index))
    output.write('\n')      
    output.write("x, avg y, STDEV, STDOM")
    output.write('\n')       
    index=0
    for j in newx:
        val1 = str(newx[index])
        val2 = str(avg_y[index])
        val3 = str(avg_y_stdev[index])
        val4 = str(avg_y_stdom[index])

        output.write(val1)
        output.write(",")
        output.write(val2)
        output.write(",")
        output.write(val3)
        output.write(",")
        output.write(val4)
        output.write('\n')       
        index=index+1
        
    output.close()
    return()
###########################################################################################################
#
#   this function is called to average data sets at a single test temperature for uncertainty quantification
#
############################################################################################################


def average_stress_strain_data(exp_list,low_strain,hi_strain):
    import numpy as np
    import math
    #import pickle
    from scipy import interpolate
    import pylab
    import matplotlib.pyplot as plt
    
    pylab.close()

    delta_strain = 0.001 #sets density of interpolation
    num = np.int_((hi_strain-low_strain)/delta_strain)
    strains=np.linspace(low_strain,hi_strain,num=num)
    
    # this is the array used to compute the average of the stress-strain curves and statistics
    avg_stress=np.array(np.linspace(low_strain,hi_strain,num=num))
    avg_stress_stdev=np.array(np.linspace(low_strain,hi_strain,num=num))
    avg_stress_stdom=np.array(np.linspace(low_strain,hi_strain,num=num))


    #this loop obtains stress, strain and temperature data for fitting - all the data goes into a single array
    index=0   
    for exp_number in exp_list:
        index = index + 1
        exp_number=str(exp_number)
        print(exp_number[0:4])
        datafile=exp_number[0:4] + '_analyzed.csv'
        data_input = open(datafile,"r")
        strain_data,stress_data,strainrate_data=np.loadtxt(data_input,skiprows=0,delimiter=',',usecols=(0,1,2),unpack=True)
        pylab.plot(strain_data,stress_data)
    # interpolate current data set at strains and accumulate into 
        newx=strains
        interp=interpolate.interp1d(strain_data,stress_data)
        interp_stress_data = interp(newx) 
        #plt.plot(strains,interp_stress_data,color=color,label=label)
        stress_sum = np.array(interp_stress_data)
        avg_stress=avg_stress.__add__(stress_sum)
    # compute average stress at each interpolated strain value and plot average result (maybe with a thicker line?)   
    avg_stress=avg_stress.__truediv__(index)
    pylab.plot(strains,avg_stress,'r--',linewidth=10.0,label="avg") 
    pylab.legend(loc='lower right')
    pylab.show()
    #pylab.close()

    # loop again through the data and assemble standard deviation and stdom's
    for exp_number in exp_list:
        exp_number=str(exp_number)
        print(exp_number[0:4])
        datafile=exp_number[0:4] + '_analyzed.csv'
        data_input = open(datafile,"r")
        strain_data,stress_data,strainrate_data=np.loadtxt(data_input,skiprows=0,delimiter=',',usecols=(0,1,2),unpack=True)
        
        newx=strains
        interp=interpolate.interp1d(strain_data,stress_data)
        interp_stress_data = interp(newx)
        error=interp_stress_data.__sub__(avg_stress)
        error=error.__pow__(2)
        avg_stress_stdev = avg_stress_stdev.__add__(error)
        #print('I am calculating averages')        
        #rint(avg_stress[100])        
        #print(avg_stress_stdev[100])
    
    avg_stress_stdev=avg_stress_stdev.__truediv__(index-1)
    avg_stress_stdev=avg_stress_stdev.__pow__(0.5)
    root_index=math.sqrt(index)
    avg_stress_stdom=avg_stress_stdev.__truediv__(root_index)
    avg_stress_2stdom=2*avg_stress_stdev.__truediv__(root_index)
    
    #compute error bands for plotting and output - this one uses two standard deviation of the mean of the data
    avg_stress_hi=avg_stress.__add__(avg_stress_2stdom)
    avg_stress_lo=avg_stress.__sub__(avg_stress_2stdom)
    
    #plt.plot(strains,avg_stress,color='y',label="avg")
    plt.plot(strains,avg_stress,color='y',label='avg')
    plt.plot(strains,avg_stress_hi,color='b',label='hi')
    plt.plot(strains,avg_stress_lo,color='g',label='lo')
    plt.legend(loc='lower right')
    
    pylab.show()



    # write output into .csv format for importing into excel or other file
    
    output = open("avg_output.csv","w")
    output.write("root index: ")
    output.write(str(root_index))
    output.write('\n')      
    output.write("strain, avg stress, stress +2, stress -2, STDEV, STDOM")
    output.write('\n')       
    index=0
    for j in strains:
        val1 = str(strains[index])
        val2 = str(avg_stress[index])
        val3 = str(avg_stress_hi[index])
        val4 = str(avg_stress_lo[index])
        val5 = str(avg_stress_stdev[index])
        val6 = str(avg_stress_stdom[index])

        output.write(val1)
        output.write(",")
        output.write(val2)
        output.write(",")
        output.write(val3)
        output.write(",")
        output.write(val4)
        output.write(",")
        output.write(val5)
        output.write(",")
        output.write(val6)
        output.write('\n')   
        index=index+1
        
    output.close()
    return()

########################################################################################
#
#   Function to quickly plot raw volts data to check polarity
#
########################################################################################

def quickplot_volts_data(expt):
    
    sg1_col,sg2_col,TC_col,N_col,N2_col=find_data_columns(expt)
        
    filename= expt + '.txt'
    #filename='3411.txt'
#    filestream=open(filename,'r')

    sg1_time,sg1_volts=load_timeseries_data(filename,8,0,sg1_col)
    sg2_time,sg2_volts=load_timeseries_data(filename,8,0,sg2_col)
    
    dt_check=sg1_time[1]-sg1_time[0]
    print('dt check:',dt_check)

    pylab.clf()    # this is needed to clear data from previous plots
    pylab.plot(sg1_time, sg1_volts, color = 'b', label='SG1')
    pylab.plot(sg2_time, sg2_volts, color = 'r', label='SG2')
    pylab.legend()
    pylab.show()

    return()



    
########################################################################################
#
#   Function to identify the data columns in a Kolsky bar text file in tension experiment
#
########################################################################################
def find_data_columns_tension(expt):
    filename='T' + expt + '.txt'
    filestream=open(filename,'r')
    #read the line in the data file that contains column headings (line 6 usually)
    i=0        
    for lines in filestream:
        #print(lines)    
        if i == 6:
            heading=lines        
            #print('i am here')
        elif i > 10:
            break
        i+=1
    filestream.close()
    
    col_titles=heading.split("\t")
    #print('col_title:',col_titles)
    sg1_col=col_titles.index("SG2 IBar")   
    sg2_col=col_titles.index("SG3 TBar Axial")
#    sg3_col=col_titles.index("TBar Poisson")
 #   kistler_col=col_titles.index("Piezo Force Sense")
    
    flag=0
    try:
        TC_col=col_titles.index("TCA")
    except:
        print('No TCA column in KBar Data File!')
        flag=flag+1
    if flag==1:    
        try:
            TC_col=col_titles.index("TC A")
        except:
            flag=flag+1
            print('No TC column in KBar Data File!')
    if flag==2:    
        try:
            TC_col=col_titles.index("TC")
        except:
            flag=flag+1
            print('No TC column in KBar Data File!')
    if flag==3:
        TC_col=0
        print('No TC temperature for this test!')
    try:
        N_col=col_titles.index("NIMPY")
    except:
        N_col=0
        print('No NIMPY column in KBar Data File!')
    try:
        N2_col=col_titles.index("NIMPY 2")
    except:
        N2_col=0
        print('No NIMPY 2 column in KBar Data File!')

    try:
        sg3_col=col_titles.index("SG3")
    except:
        sg3_col=0
        print('No SG3 column in KBar Data File!')
    try:
        sg3_col=col_titles.index("TBar Poisson")
    except:
        sg3_col=0
        print('No TBar Poisson column in KBar Data File!')
    try:
        kistler_col=col_titles.index("Piezo Force Sense")
    except:
        kistler_col=0
        print('No Piezo Force Sense column in KBar Data File!')


    #print('SG1_col:',sg1_col)

    #return()
    return(sg1_col,sg2_col,sg3_col,kistler_col,TC_col,N_col,N2_col)


########################################################################################
#
#   Function to compute wave time from bars apart and bar to bar tests
#
########################################################################################
def calculate_wave_timing(test,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window=600):

    sg1_col,sg2_col,dummy1,dummy2,dummy3=find_data_columns(test)  # assume these calibrations will be the same for bars apart and together tests

    hi_rate_data=test + '.txt'

    #
    #sg1_col=int(sg1_col)
    #sg2_col=int(sg2_col)
    sg1_time,sg1_volts=load_timeseries_data(hi_rate_data,8,0,sg1_col)
    sg2_time,sg2_volts=load_timeseries_data(hi_rate_data,8,0,sg2_col)
    sg1_volts=zero_timeseries_data(sg1_volts,100)# 100 is the number of data points to average
    sg2_volts=zero_timeseries_data(sg2_volts,100)# 100 is the number of data points to average  
            
    sg1_cal_volts=2*(sg1_hi-sg1_lo)/1000.0
    sg2_cal_volts=2*(sg2_hi-sg2_lo)/1000.0
    
    sg1_cal=0.004854/sg1_cal_volts
    sg2_cal=0.004854/sg2_cal_volts 
    
    sg1_strain=volts_to_strain_data(sg1_volts,sg1_cal)
    sg2_strain=volts_to_strain_data(sg2_volts,sg2_cal)
    

    #Identify data pulses*********************************************   
    #pylab.clf()
    inc_start,inc_time,inc=find_data_pulse4(sg1_time,sg1_strain,window,1) #last number is 1=incident wave, 0=refl or trans
#    quickplot_xydata(inc_time,inc)
    
    refl_start,refl_time,refl=find_data_pulse4(sg1_time,sg1_strain,window,2)
#    quickplot_xydata(refl_time,refl)
    
    trans_start,trans_time,trans=find_data_pulse4(sg2_time,sg2_strain,window,3)
#    quickplot_xydata(trans_time,trans)

    refl_delay=refl_start-inc_start
    trans_delay=trans_start-inc_start
    
    #print(inc_start,refl_start,trans_start)    
    #refl_delay=(refl_start-inc_start)*1000000.0
    #trans_delay=(trans_start-inc_start)*1000000.0
    #print(refl_delay,trans_delay)

    
   #inc=inc*(-1.0)
   # trans=trans*(-1.0)

    #pylab.clf()    # this is needed to clear data from previous plots
    #pylab.plot(inc_time, inc, color = 'b', label='INC')
    #pylab.plot(inc_time, refl, color = 'r', label='REFL')
    #pylab.plot(inc_time, trans, color = 'y', label='TRANS')
    #pylab.legend()
    #pylab.show()
    
    #print(inc_time)
    
    return(inc_start,refl_start,trans_start)


########################################################################################
#
#   Function to perform KBar Mechanical Analysis and deliver stress-strain for a tension experiment using DIC
#
########################################################################################
def calculate_tensile_stress_strain(test,gage_w,gage_thk,gage_len,sg1_hi,sg1_lo,sg2_hi,sg2_lo,sg3_hi,sg3_lo,window=600,refl_delay=0.000340,trans_delay=0.000346,DIC=1,num_gages=3,tgage_type=0,spec_type=0,kistler=0,ibar_modulus=190,tbar_modulus=170,DIC_first=1,DIC_last=100,DIC_frame_rate=75000,DIC_col=7):

    from scipy import interpolate
    

    # conditionals   
    #      spec_type, or specimen type. 0=sheet, 1=round
    #      num_gages, or number of gages. 0=3 gages (ibar, tbar axial, tbar Poisson) ; 1=2 gages (ibar,tbar axial or combined)
    #      kistler. 0=no, 1=yes (need to include Kistler data calibration)
    #      DIC   0=yes, 1=no. Is there DIC data available
    #      tgage_type, or trans bar gage type. 0:axial+Poisson integrated; 1: axial only; 2:axial and Poisson separate
   
   
    gage_label=str(int(gage_len*1000.0))
    sg1_col,sg2_col,sg3_col,kistler_col,TC_col,N_col,N2_col=find_data_columns_tension(test)
    
   # hi_rate_data='T'+test + '_gage'+ gage_label+'.txt'
    hi_rate_data='T'+test + '.txt'
    output_file=test+'_gage'+gage_label+'_analyzed.csv'
    #
    sg1_time,sg1_volts=load_timeseries_data(hi_rate_data,8,0,sg1_col)
    sg2_time,sg2_volts=load_timeseries_data(hi_rate_data,8,0,sg2_col)
    
    sg1_volts=zero_timeseries_data(sg1_volts,100)# 100 is the number of data points to average
    sg2_volts=zero_timeseries_data(sg2_volts,100)# 100 is the number of data points to average  
            
    sg1_cal_volts=2*(sg1_hi-sg1_lo)/1000.0
    sg2_cal_volts=2*(sg2_hi-sg2_lo)/1000.0
    
    sg1_cal=0.004854/sg1_cal_volts
    sg2_cal=0.004854/sg2_cal_volts    
    sg1_strain=volts_to_strain_data(sg1_volts,sg1_cal)
    sg2_strain=volts_to_strain_data(sg2_volts,sg2_cal)
    
    #quickplot_xydata(sg2_time,sg1_volts)
    #Identify incident and trans data pulses*********************************************

    dt=sg1_time[1]-sg1_time[0]    
    inc_start,inc_time,inc=find_data_pulse4(sg1_time,sg1_strain,window,1) #last number is 1=incident wave, 0=refl or trans
    inc_start_index=int((inc_start-sg1_time[0])/dt)

#    quickplot_xydata(inc_time,inc)
    
    # pick reflected pulse based on previously determined delay time specified in the input file
    refl_start_index=int((inc_start+refl_delay-sg1_time[0])/dt)
    refl_time=sg1_time[refl_start_index:refl_start_index+window]
    refl=sg1_strain[refl_start_index:refl_start_index+window]

 #   plt.clf()    
 #   plt.plot(refl_time,refl,color='y',label='refl')
 #   plt.show()
   
 #   print(len(sg1_strain))
 #   print(window)
 #   print(len(refl))
 #   print(len(refl_time))

############################################
#    adjust incident pulse to make sure the inc pulse strain is 0 (within noise) at time zero, 
#       for cases of overlapping reflected strain with 1 m striker. should not affect other cases
########################################

 #   zero_tol=numpy.average(sg1_strain[0:100])
    zero_tol=0.00001   
   # print(zero_tol)
    # Uses a single-sided search routine that assumes the SG1 pulse is negative and initially the first strain value
    #   in the pulse is somewhere on the rising slope and not at zero. will not work if this is not true. need better robustness
    search_range=int(window/4)
    search_points=numpy.arange(0,search_range,1)
    shift_opt=0
    inc_opt=inc
    inc_time_opt=inc_time
    for i in search_points:
        inc_zero=abs(inc_opt[0])
        if inc_zero>zero_tol:
            shift_opt=shift_opt-1
        inc_opt=sg1_strain[inc_start_index+shift_opt:inc_start_index+shift_opt+window]
        inc_time_opt=sg1_time[inc_start_index+shift_opt:inc_start_index+shift_opt+window]
        
    inc=inc_opt   #set reflected wave to be used to analyze pulses equal to the optimal wave.

############################################
#    adjust reflected pulse to make sure the refl  strain is 0 (within noise) at time zero, 
#       for cases of overlapping reflected strain with 1 m striker. should not affect other cases
########################################

        
    #zero_tol=numpy.average(sg1_strain[0:100])
    zero_tol=0.00001   
   # print(zero_tol)
    # this search is two-sided and assumes the 1 m reflected strain pulse overlaps such that there is no flat zero
    #   level between incident and reflected pulses. Need to see how this works with non-overlapping pulses
    search_range=int(window/4)
    search_points=numpy.arange(0,search_range,1)
    shift_opt=0
    refl_opt=refl
    refl_time_opt=refl_time
    for i in search_points:
        refl_zero=refl_opt[0]
        if refl_zero>zero_tol:
            shift_opt=shift_opt-1
        if refl_zero<zero_tol:
            shift_opt=shift_opt+1
        #print('refl_zero:',refl_zero)
        #print('shift opt:',shift_opt)
        refl_opt=sg1_strain[refl_start_index+shift_opt:refl_start_index+shift_opt+window]
        refl_time_opt=sg1_time[refl_start_index+shift_opt:refl_start_index+shift_opt+window]

  #  print(len(refl_opt))

    refl=refl_opt   #set reflected wave to be used to analyze pulses equal to the optimal wave.
    
    #shift the incident wave the same amount as the reflected wave?
   # inc_start_index=int((inc_start-sg1_time[0])/dt)
   # inc_opt=sg1_strain[inc_start_index+shift_opt:inc_start_index+shift_opt+window]
   # inc=inc_opt
    
#    refl_start,refl_time,refl=find_data_pulse4(sg1_time,sg1_strain,window,2)
#    quickplot_xydata(refl_time,refl)
    
    trans_start,trans_time,trans=find_data_pulse4(sg2_time,sg2_strain,window,3)
    trans_start_index=int((trans_start-sg2_time[0])/dt)
    #print('trans_start_index:',trans_start_index)
    #print('sg2 strain at start index:',sg2_strain[trans_start_index])
   # quickplot_xydata(trans_time,trans)
   
    if num_gages==3:  #means there are three separate strain gages, meaning the transbar has a separate Poisson gage (transP)
        sg3_time,sg3_volts=load_timeseries_data(hi_rate_data,8,0,sg3_col)
        sg3_volts=zero_timeseries_data(sg3_volts,100)# 100 is the number of data points to average  
        sg3_cal_volts=2*(sg3_hi-sg3_lo)/1000.0
        sg3_cal=0.004854/sg3_cal_volts    
        sg3_strain=volts_to_strain_data(sg3_volts,sg3_cal)
        transP_start,transP_time,transP=find_data_pulse4(sg3_time,sg3_strain,window,3)
        #quickplot_xydata(transP_time,transP)
    else:
        transP=0.0
        
    if DIC == 0:   
        # Compute DIC strain by averaging over the nominal gage length, and interpolate result
#       to obtain a separate strain-time "wave" to identify
        DIC_time,DIC_strain,dummy=calculate_dic_strain_time(expt,DIC_first,DIC_last,DIC_frame_rate,gage_len,DIC_col) #area averaged DIC strain
        #DIC_time,DIC_strain=calculate_DIC_strain_virtualgage(expt,DIC_first,DIC_last,DIC_frame_rate,gage_len)
        
        DIC_length=len(DIC_time)        
        print('DIC Length:',DIC_length)    
        DIC_time_final=DIC_time[DIC_length-1]
        print(DIC_time)
        print('DIC time final:',DIC_time_final)
       # interpolate DIC data using time window of stress pulse
        time=numpy.arange(0,DIC_time_final,dt)    # creates time series based on dt for strain gages to interpolate DIC data which is at a lower data rate (larger native dt). this array is smaller than the wave analysis window
        newx=time
        interp=interpolate.interp1d(DIC_time,DIC_strain)
        interp_DIC_strain = interp(newx) 
        # find DIC start time
        DIC_start,DIC_time,DIC_newstrain=find_DIC_start(time,interp_DIC_strain,window) #last number is 1=incident wave, 0=refl or trans
    #    print('DIC start:',DIC_start)
    else:
        DIC_strain=0.0


    

    if kistler==1:  #means there are three separate strain gages, meaning the transbar has a separate Poisson gage (transP)
        kistler_time,kistler_volts=load_timeseries_data(hi_rate_data,8,0,kistler_col)
        kistler_volts=zero_timeseries_data(kistler_volts,100) # need to make sure this is necessary
        kistler_cal=1.0  # this is the volts-to-force calibration number that needs to be input here
        kistler_force=kistler_volts*kistler_cal
        kislter_start,kistler_time,kistler_force=find_data_pulse4(kistler_time,kistler_force,window,3) #not sure if the wave start will work here
#        quickplot_xydata(trans_time,trans)
    else:
        kistler_force=0.0


    common_time=DIC_time # this is for plotting all time-data together. if the start find rountines are consistent, the waves should be aligned
    
     # Perform Kolsky analysis to get stress-strain curves*************************************

    #pre-pulse is the initial zero level time for the pulse windows - in this code it 
    # no longer needs to be at a specific time because the compliance correciton is eliminated
    dt=common_time[1]-common_time[0]
    pre_pulse=int(0.0001/dt) #this needs to be tuned 
    #print('dt:',dt)
    
 
    #calculate sample strain, strain rate and stress histories using aligned waves    
    # First define analysis constants - all lengths in m

    ibar_diameter=0.026
    ibar_area=3.14159/4*ibar_diameter*ibar_diameter
 #   ibar_density=8048.0
    ibar_density=8070.0    
    iwave_speed=ibar_modulus*1.0e9/ibar_density
    iwave_speed=iwave_speed.__pow__(0.5)


    tbar_diameter=0.02
    tbar_area=3.14159/4*tbar_diameter*tbar_diameter
    tbar_density=2800.0
    twave_speed=tbar_modulus*1.0e9/tbar_density
    twave_speed=twave_speed.__pow__(0.5)
    
    if spec_type == 0:  #this designates a rectangular gage section
        sample_area=gage_w*gage_thk
    else:
        sample_area=3.14159/4*gage_w*gage_w # this is a round gage section
  
  
    #define strain rate and strain increment and stress vs time arrays using 1-wave analysis
    eng_strain=numpy.arange(0,window,1)*0.0

    eng_stress=tbar_area*tbar_modulus/sample_area*trans*1000 #[MPa]
    #true_strain=DIC_strain
    
    
    #this uses a 2 wave equation to estimate strain rate for our Tensile bar with different materials for ibar and tbar 
   # print('gage len:',gage_len)
    eng_strain_rate=(iwave_speed*(refl-inc)-twave_speed*trans)/gage_len
    eng_strain_increment=eng_strain_rate*dt
    # sum up engineering strain    
    sum=0.0
    #sum up total engineering strain
    window=int(window)
    for i in range(0,window,1):             
        sum+=eng_strain_increment[i] 
        eng_strain[i]=sum              

    #pylab.plot(eng_strain, eng_stress, color = 'r', label='Before Opt')
    

###############  
# Use the anticipated specimen elastic modulus to determine the optimal start time of the transmission wave
###############   
    theoretical_sample_modulus = 200.0  #in GPa
    theoretical_yield_stress = 1.0  #in GPa
    theoretical_stress_in_elastic_region = 0.75*theoretical_yield_stress
    theoretical_strain = theoretical_stress_in_elastic_region/theoretical_sample_modulus
    
    theo_strain=[0,.01]
    theo=theo_strain[1]*theoretical_sample_modulus*1000.0
    theo_stress=[0,theo]
    #pylab.plot(theo_strain, theo_stress, color = 'r', label=plot_label)

    
    stress_tol=0.020 #stress in GPa   
    strain_tol=0.0001 #strain for searching eng strain array 
    search_range=int(window/10) #select some fraction of the analysis window to calculate the optimal delay time for the trans wave
    search_points=numpy.arange(0,search_range,1)
    shift_opt=0   
    trans_opt=trans
    strain_index=0
    strain_index_final=0
    count=0
    for i in eng_strain:
        strain_check=abs(eng_strain[count]-theoretical_strain)
       # print('strain check:',strain_check)
        if strain_check < strain_tol:
            strain_index_final=count
           # print('final strain index found')
        strain_index+=1
        count+=1
   # print('strain_index after search:',strain_index_final) 
   # print('eng strain at strain index:',eng_strain[strain_index_final])
   # print('eng_stress_at_index:',eng_stress[strain_index_final])    
   # print('theor eng_stress_at_index:',theoretical_stress_in_elastic_region)    

    for j in search_points:    
        stress_check=eng_stress[strain_index_final]-theoretical_stress_in_elastic_region*1000.0 #in MPa
        if stress_check>stress_tol:
            shift_opt=shift_opt-1
        if stress_check<stress_tol:
            shift_opt=shift_opt+1
        #print('stress_check:',stress_check)
        #print('shift opt:',shift_opt)
        trans_opt=sg2_strain[trans_start_index+shift_opt:trans_start_index+shift_opt+window]
        eng_stress=tbar_area*tbar_modulus/sample_area*trans_opt*1000 #[MPa]

    trans=trans_opt   #set reflected wave to be used to analyze pulses equal to the optimal wave.
    
    # compute and plot wave equilibrium after adjusting for specimen modulus
    inc_plus_refl=inc+refl
 
####################
# once waves are aligned using the known Modulus method, adjust true strain-time using the DIC true strain averaged over the gage length
####################
    #first calculate true strain from engineering strain
    true_strain=numpy.arange(0,window,dtype=numpy.float)  #this just initializes a sequence of the right size   
    true_strain=numpy.log(1.0+eng_strain)

# If this is a fracture experiment, line-up wave and DIC strain-time based on wave start finder routine. Sometimes this doesn not work well however   
    true_strain_adj=true_strain
    DIC_strain_start,DIC_time,DIC_newstrain=find_DIC_start(common_time,DIC_newstrain,window) #using same routine to adjust true-strain time curve to DIC time curve
    print('DIC strain start:',DIC_strain_start) 
    true_strain_start,true_strain_adj_time,true_strain_adj=find_DIC_start(common_time,true_strain,window) #using same routine to adjust true-strain time curve to DIC time curve
    print('true strain start:',true_strain_start)    
    DIC_wave_time_shift=int((DIC_strain_start-true_strain_start)/dt)
    print('DIC Wave time shift:',DIC_wave_time_shift)
    
    #DIC_wave_time_shift=0

# shift DIC strain-time to correspond to wave strain-time
    DIC_strain_shifted=numpy.arange(0,window,1)*0.0   
    if DIC_wave_time_shift>=0:
        DIC_strain_shifted[0:window-DIC_wave_time_shift]=DIC_newstrain[DIC_wave_time_shift:window]
        DIC_strain_shifted[window-DIC_wave_time_shift+1:window]=0.0
    if DIC_wave_time_shift<0:
        DIC_strain_shifted[0:-DIC_wave_time_shift]=0.0
        DIC_strain_shifted[-DIC_wave_time_shift+1:window]=DIC_newstrain[0:window+DIC_wave_time_shift-1]

  #  plt.figure()
 #   plt.plot(common_time,DIC_newstrain,'o')
 #   plt.plot(common_time,DIC_strain_shifted,'+')   
 #   plt.plot(common_time,true_strain,'-')


  # further align DIC and wave strain-time data by looking at peaks - only for tension tests that do not fracture  
  #find peak values and indices of DIC strain and un-scaled wave strain 
  #  peak_index=numpy.argmax(DIC_strain_shifted)
  #  peak_strain=numpy.amax(DIC_strain_shifted)
  #  peak_index_wave=numpy.argmax(true_strain)
  #  peak_strain_wave=numpy.amax(true_strain)
   
  #  peak_shift=peak_index_wave-peak_index #this is the time index difference in the wave and DIC strain-time peaks
  #  print('peak shift:',peak_shift)
    
  #  peak_shift=0.0

# now shift the DIC strain-time to line up the peaks with the wave data   
 #   DIC_strain_shifted2=numpy.arange(0,window,1)*0.0   
 #   if peak_shift>=0:
 #       DIC_strain_shifted2[0:peak_shift]=0.0
 #       DIC_strain_shifted2[peak_shift+1:window]=DIC_strain_shifted[0:window-peak_shift-1]
 #   if peak_shift<0:
 #       DIC_strain_shifted2[0:window+peak_shift]=DIC_strain_shifted[-peak_shift:window]
 #       DIC_strain_shifted2[window+peak_shift+1:window]=0.0

    DIC_strain_shifted2=DIC_strain_shifted

 #   plt.plot(common_time,DIC_strain_shifted2,'^')       
 #   plt.show()
        
#    print('DIC peak strain:',peak_strain)    
#    print('DIC peak strain index:',peak_index)
        
    # compute average strain rate and use that to estimate the linear portion of the strain-time curve
    peak_index=numpy.argmax(DIC_strain_shifted2)
    peak_strain=numpy.amax(DIC_strain_shifted2)
 
    num_time_points = 200 # this is the number of time points used for the strain rate calculation. For 2 MHz this is 100 microseconds 
    estimated_strain_rate = (DIC_strain_shifted2[peak_index]-DIC_strain_shifted2[peak_index-num_time_points])/(float(num_time_points)*dt)
    estimated_range_for_comparison = int(peak_strain/estimated_strain_rate/dt)
    #print('length of DIC_newstrain:',len(DIC_newstrain))    
    #print('peak strain:',peak_strain)
    #print('dt:',dt)
   
    #print('peak strain index:',peak_index)
    #print('est strain rate:',estimated_strain_rate)
    #print('index for comparison:',estimated_range_for_comparison)
    
    # now adjust wave true-strain vs time to match DIC strain up to peak DIC strain (before fracture)

    # this tries to match the true-strain vs time slope between wave strain and DIC strain

    # this adjustment simply looks at the residual between DIC strain and wave strain for a linear portion of the strain-time curve   
    strain_range=numpy.arange(peak_index-estimated_range_for_comparison,peak_index,1)
    #opt_scale=numpy.arange(0.5,1.0,0.01)
    opt_scale=numpy.arange(0.5,1.2,0.01)
    DIC_plot_strain=DIC_strain_shifted2[peak_index-estimated_range_for_comparison:peak_index] #this is just for plotting
   #print(len(DIC_plot_strain))
   # print(len(strain_range))
    residual_sum=0.0
    residual_sum_old=1000.0
 #   plt.figure()
  #  plt.plot(strain_range,DIC_plot_strain,'o',label='DIC Strain')
    for j in opt_scale:
        residual_sum=0.0       
        for k in strain_range:
            residual=true_strain[k]*j-DIC_strain_shifted2[k]  # j is the scaling factor for wave strain
            residual=residual*residual
            residual=numpy.power(residual,0.5)
            residual_sum+=residual
        if residual_sum<residual_sum_old:
            residual_sum_old=residual_sum
            j_opt=j
        #print('residual_sum_old:',residual_sum_old)
        #print('residual_sum:',residual_sum)
        scaled_strain=true_strain[peak_index-estimated_range_for_comparison:peak_index]*float(j)
       # print('leng of scaled strain:',len(scaled_strain))
      #  print('leng of strain range:',len(strain_range))

#        plt.plot(strain_range,scaled_strain,'+',label='scaled wave')
#    plt.legend()    
#    plt.show()        
    print('j_opt:',j_opt)        
    true_strain_final=true_strain*j_opt            #compute final wave true strain and plot with original time
        
  #  print('DIC wave time shift:',DIC_wave_time_shift) 
    


# these are the final engineering and true stress-strain results based on DIC adjusted values
    eng_strain=numpy.exp(true_strain_final)-1.0
    true_stress=(1.0+eng_strain)*eng_stress
        
    #estimate true strain rate by computing instantaneous (current) sample length
    window=int(window)    
    true_strain_rate=numpy.arange(0,window,dtype=numpy.float)
    eng_strain_rate=numpy.arange(0,window,dtype=numpy.float)

    for i in range(0,window,1):             
        if i>0:        
            true_strain_rate[i]=(true_strain_final[i]-true_strain_final[i-1])/dt
            eng_strain_rate[i]=(eng_strain[i]-eng_strain[i-1])/dt
            
    #pick a stress value to see if it is positive or negative
   # stress_time=int(0.0001/dt)
   # stress_interrogation = true_stress[stress_time]
   #print('True stress at 100 µs :',stress_interrogation)
    
 
     
    #show plots showing various data 
    plt.clf()
    Title_1='Wave Equilibrium' 
    Title_2='DIC Strain Adjustment' 
    Title_3='Engineering Stress-Strain'
    Title_4='True Stress-Strain'
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].set_title(Title_1)
    axarr[0, 1].set_title(Title_2)
    axarr[1, 0].set_title(Title_3)
    axarr[1, 1].set_title(Title_4)
    elastic_label=str(theoretical_sample_modulus)+' GPa'

#    axarr[0, 0].plot(common_time,inc_plus_refl,color = 'b', label='inc+refl')
    axarr[0, 0].plot(common_time,trans,color = 'y', label='trans')
    axarr[0, 0].plot(common_time,inc,color = 'b', label='inc')
    axarr[0, 0].plot(common_time,refl,color = 'r', label='refl')
    axarr[0, 0].legend()

    axarr[0, 1].plot(common_time,true_strain,color = 'y', label='original true strain')
    axarr[0, 1].plot(common_time,DIC_strain_shifted2,color = 'b', label='DIC strain')       
    axarr[0, 1].plot(common_time,true_strain_final,color = 'r', label='final true strain')
    axarr[0, 1].legend()

    axarr[1, 0].plot(eng_strain, eng_stress, color = 'b', label='Eng Stress')
    axarr[1, 0].plot(eng_strain, eng_strain_rate, color = 'r', label='Eng Strain Rate')
    axarr[1, 0].plot(theo_strain, theo_stress, color = 'r', label=elastic_label)
    axarr[1, 0].legend()   
   
    axarr[1, 1].plot(true_strain, true_stress, color = 'b', label='True Stress')
    axarr[1, 1].plot(true_strain, true_strain_rate, color = 'r', label='True Strain Rate')
    axarr[1, 1].plot(theo_strain, theo_stress, color = 'r', label=elastic_label)
    axarr[1, 1].legend()   
    plt.show()
    
   

# write true and engineering stress-strain data to output file
    data_outfile = open(output_file,"w")
    data_outfile.write('DIC Strain Adj Factor:,')
    data_outfile.write(str(j_opt))
    data_outfile.write('\n')
    data_outfile.write('Time,True Strain,True Stress,True Strain Rate, Eng Strain, Eng Stress, Eng Strain Rate')
    data_outfile.write('\n')
    data_to_write=numpy.column_stack((common_time,true_strain_final,true_stress,true_strain_rate,eng_strain,eng_stress,eng_strain_rate))
    numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
    data_outfile.close()  

    
    #print('test number:',test)
    return(true_strain,true_stress,true_strain_rate)

    #return()
    
    ###################################################################################################

###################################################################################################
#
#            calculate average strain-time in the gage section of a tensile specimen from DIC data
#
###################################################################################################    

def calculate_dic_strain_time(expt,DIC_first,DIC_last,DIC_frame_rate,gage_l,DIC_col):

    #find initial gage length from first data file
    file_root='T'+str(int(expt))+'_0_'
    # read in first correlation data set to determine initial correlation (gage) length
    filenumber=str(int(DIC_first)).zfill(6)
    filename=file_root + filenumber + '.csv' 
    initial_DIC_datafile = open(filename,"r")
    dummy, initial_y=numpy.loadtxt(initial_DIC_datafile,delimiter=",",skiprows=1,usecols=(0,1),unpack=True)
    a=numpy.array(initial_y) #must cast arrays as numpy arrays before doing arithmetic
    y_max=a.max()
    y_min=a.min()
    initial_gage=y_max-y_min
    print(initial_gage)
    initial_DIC_datafile.close()
    
    #define interframe time
    frame_rate=float(DIC_frame_rate)
    dt = 1/frame_rate
    
    #define desired gage length in mm and other data file attributes needed
    desired_gage=gage_l*1000.0
    first_file_number=int(DIC_first)
    last_file_number=int(DIC_last)
    data_col_to_avg=int(DIC_col)
    
    #create DIC strain-time array for return to stress-strain calculator
    time_points=DIC_last-DIC_first
    time_final=float(time_points)*dt
    
    DIC_time=numpy.arange(0.0,time_final,dt)
    DIC_strain=numpy.arange(0.0,time_final,dt)
    DIC_strain_std=numpy.arange(0.0,time_final,dt)
    
    #define y limits to average within
    y_max_new = desired_gage
    y_max_new = y_max_new/2.0
    
     #open output file
    outfilename=file_root + 'avg.csv'
    output=open(outfilename,"w")
    output.write('initial gage: ')
    out1=str(initial_gage)
    output.write(out1)
    output.write('  desired gage: ')
    out1=str(desired_gage)
    output.write(out1)
    output.write('\n')
    output.write('filenumber, time [s], avg, std')
    output.write('\n')
    
    # loop through specific number of DIC data files and average the column selected
    # over the specified gage length
    time = 0.0
    incr=0
    for i in range(int(first_file_number),int(last_file_number)):
        filenumber=str(i).zfill(6)
        filename=file_root + filenumber + '.csv' 
        filestream = open(filename,"r")
        y_index, var_to_avg=numpy.loadtxt(filestream,delimiter=",",skiprows=1,usecols=(1,data_col_to_avg),unpack=True)
        
        # this structure is supposed to reduce the total DIC data into an array for averaging
        # by using a boolean array that has in it index values that satisfy the inequality
        j = numpy.array(len(y_index)*[False],bool)
        #print(y_index)
        count=0
        for k in y_index:
            #print(y_index[count])
            if(abs(y_index[count])<y_max_new):
                #print(count)
                j[count]=True
            count=count+1
        
        #here the sub-array called 'avg_array' is pulled from the original data array at 
        # index locations corresponding to true values in the j array, that were computed in the loop above hopefully
       # print(time,dt,frame_rate)        
        a=numpy.array(var_to_avg[j])
        avg_var=numpy.average(a)
        std_var=numpy.std(a)
        
        DIC_time[incr]=time
        DIC_strain[incr]=avg_var
        DIC_strain_std[incr]=std_var
        
        #convert data output items to strings for writing to output file
        filenum_str=str(filenumber)
        time_str=str(time)
        avg_var_str=str(avg_var)
        std_var_str=str(std_var)
        output.write(filenum_str)
        output.write(',')
        output.write(time_str)
        output.write(',')   
        output.write(avg_var_str)
        output.write(',')
        output.write(std_var_str)
        output.write('\n')
        filestream.close
        
        time = time + dt
        incr+=1
           
    output.close()

    
    return(DIC_time,DIC_strain,DIC_strain_std)
    
    
#--------------------END of DIC Analysis    
###########################################################################################################################
#
#   Function to identify DIC strain-time pulse start time for tension KBar experiments
#
############################################################################################################################
def find_DIC_start(time_data,data_data,window=400):

    data=numpy.array(data_data)
    time=numpy.array(time_data)
   # print(len(data))
   # print(len(time))
    dt=time_data[1]-time_data[0]
    #delay=0
    

    # define dummy arrays for DIC strain-time data with same dimensions as strain pulse window size
    max_time=window*dt
    time_dummy=numpy.arange(0.0,max_time,dt)
    data_dummy=numpy.arange(0.0,max_time,dt)*0.0    #initialize this array as all zeros        
    
    checksum_min= 150  # number of time points a pulse must be above a certain level to be considered a pulse
   # avg_level=numpy.average(data[0:5])
    #std_level=numpy.std(data[0:5])
    
    avg_level=0.00001
    #print('avg level:',avg_level)

    check_level=0.001  # checks level in strain to capture DIC strain-time pulse
  # delay=100
    delay=0
    increment=0
    checksum=0
    flag = 0
    
    #provisional start time is the first data point in time-data
    provisional_start_time=time_data[0]    
    start_time=provisional_start_time
    start=0
    for i in data:  
        if increment>delay:
            level_check = (data[increment]-avg_level)
            level_check=abs(level_check)
            slope = data[increment]-data[increment-1]
            slope = abs(slope/dt)
            #print ('level check',level_check,check_level)
            if slope < 5000: #this determines how many points are above the check level criteria
                if level_check > check_level:
                    if flag == 0:        
                        start = increment
                        flag = 1
                        print ('flag set to 1')
                    checksum+=1
        increment+=1
    start_time=start*dt
    incr=0    
    index=0
   # print(len(data_dummy))
   # print(len(data))
    #this part adds zero padding to make sure what??   
    for j in data:
        if incr >= start:
            if delay <= start:
                data_dummy[index+delay]=data[incr-delay]
            #data_dummy[incr-delay]=1
                index+=1
            else:
                data_dummy[index+delay]=data[incr]
              
            #print(index)
            #print(incr-delay)
        incr+=1   
    
   # print('provisional start time:',provisional_start_time)            
   # plt.figure()
   # plt.plot(time,data)
   # plt.show    
    
    return(start_time,time_dummy,data_dummy)

##############################################################################################
# This function will compute the average surface strain from Vic3D as a function
#    of length along the load axis to check uniformity of strain
###############################################################################################
def calculate_DIC_strain_virtualgage(expt,DIC_first,DIC_last,DIC_frame_rate,gage_l):
    #find initial gage length from first DIC data file
    file_root='T'+str(int(expt))+'_0_'
    # read in first correlation data set to determine initial correlation (gage) length
    filenumber=str(int(DIC_first)).zfill(6)
    filename=file_root + filenumber + '.csv' 
    initial_DIC_datafile = open(filename,"r")
    initial_x, initial_y=numpy.loadtxt(initial_DIC_datafile,delimiter=",",skiprows=1,usecols=(0,1),unpack=True)
    initial_x_array=numpy.array(initial_x) # used to calculate DIC subset spacing    
    initial_y_array=numpy.array(initial_y) #must cast arrays as numpy arrays before doing arithmetic
    y_max=initial_y_array.max()
    y_min=initial_y_array.min()
    initial_gage=y_max-y_min
    print(initial_gage)
    
    #estimate DIC subset spacing
    subset_spacing=initial_x_array[1]-initial_x_array[0]
    subset_spacing=abs(subset_spacing)
    print('subset spacing:',subset_spacing)    
    initial_DIC_datafile.close()
    
    #define interframe time
    frame_rate=float(DIC_frame_rate)
    dt = 1/frame_rate
    
    #define desired gage length in mm and other data file attributes needed
    desired_gage=gage_l*1000.0
    first_file_number=int(DIC_first)
    last_file_number=int(DIC_last)
  
    #define y limits to average within
    y_max_new = desired_gage
    y_max_new = y_max_new/2.0

    #open output file
    outfilename=file_root + '_DIC_vgage.csv'
    output=open(outfilename,"w")
    output.write('initial gage: ')
    out1=str(initial_gage)
    output.write(out1)
    output.write('  desired gage: ')
    out1=str(desired_gage)
    output.write(out1)
    output.write('\n')
    output.write('filenumber, time (microsec), y1_avg (mm), y2_avg (mm), Initial Gage (mm), U1-U2 (mm), Eng_strain, True_Strain ')
    output.write('\n')

    # loop through specific number of DIC data files and average the column selected
    # over the specified gage length
    time = 0.0
    num_frames = DIC_last - DIC_first
    
    #initialize arrays for output time-strain data
    DIC_time=numpy.arange(0,num_frames,1)*0.0
    DIC_strain=numpy.arange(0,num_frames,1)*0.0
    

    ####################################################
    #
    #   delta_y is the value used to identify a new row of DIC points along the length of the specimen. x is usually 
    #       transverse to the axis, and y is along the axis. Data averaging occurs over x for constant values of y,
    #       so we need to find when the end of a row of x's are by identifying a subset-level jump in the value of y
    #       This number is not exact due to slight tilting of the data. Should probably rotate the data to account for this but 
    #       the angle is usually small
    #
    ##################################################
    delta_y = 0.9*subset_spacing   
    frame_count=0

    for i in range(first_file_number,last_file_number):
        filenumber=str(i).zfill(6)
        filename=file_root + filenumber + '.csv' 
        filestream = open(filename,"r")
        y_index, V_displ=numpy.loadtxt(filestream,delimiter=",",skiprows=1,usecols=(1,4),unpack=True) #picking Y and V displacement in DIC data file
        
         # this structure is supposed to reduce the total DIC data into an array for averaging
        # by using a boolean array that has in it index values that satisfy the inequality   
    #    j = numpy.array(len(y_index)*[False],bool)  
       
        
        row_count=1
        point_count=0
        flag_first=1
        tol=0.5*delta_y
        low_count = 0
        low_pos_sum=0.0
        low_displ_sum=0.0
        hi_count = 0
        hi_pos_sum=0.0
        hi_displ_sum=0.0
        for k in y_index: # loop through all the DIC points in this frame
            #print("count:",count)
            k_search_low = k - y_max_new # for +y_max_new this is close to zero
            k_search_hi = k + y_max_new # for -y_max_new this is close to zero

           # print('k search low:',k_search_low)
            #print('k search hi:',k_search_hi)
            
            if(abs(k_search_low) <= tol): # we have the first of the positive y max new's
                low_count+=1
                low_pos_sum+=y_index[point_count]
                low_displ_sum+=V_displ[point_count]
                
            if(abs(k_search_hi) <= tol): # we have the first of the positive y max new's
                hi_count+=1
                hi_pos_sum+=y_index[point_count]
                hi_displ_sum+=V_displ[point_count]
                
            point_count+=1
            
        # compute position and displacement averages for this frame           
        low_pos=low_pos_sum/float(low_count)
        low_displ=low_displ_sum/float(low_count)
        hi_pos=hi_pos_sum/float(hi_count)
        hi_displ=hi_displ_sum/float(hi_count)
       # print('lo pos [mm]:',low_pos)            
       # print('hi pos [mm]:',hi_pos)            
       # print('lo displ [mm]:',low_displ)            
       # print('hi displ [mm]:',hi_displ)            
            
        initial_gage = low_pos - hi_pos
        gage_length_change = low_displ - hi_displ
        
        eng_strain = gage_length_change/initial_gage
        true_strain = numpy.log(1+eng_strain)
            
       # print('initial_gage:',initial_gage)
       # print('gage_length_change:',gage_length_change)
                
        # print data to output  
        filenum_str=str(filenumber)
        time_str=str(time)
        low_pos_str=str(low_pos)
        hi_pos_str=str(hi_pos)
        initial_gage_str=str(initial_gage)
        gage_length_change_str=str(gage_length_change)
        eng_strain_str=str(eng_strain)
        true_strain_str=str(true_strain)
        output.write(filenum_str)
        output.write(',')
        output.write(time_str)
        output.write(',')
        output.write(low_pos_str)
        output.write(',')
        output.write(hi_pos_str)
        output.write(',')   
        output.write(initial_gage_str)
        output.write(',')
        output.write(gage_length_change_str)
        output.write(',')
        output.write(eng_strain_str)
        output.write(',')
        output.write(true_strain_str)
        output.write('\n')
        
        print(frame_count,time,true_strain)
       # store time and DIC strain calculation in output arrays
        DIC_time[frame_count]=time
        DIC_strain[frame_count]=true_strain
        frame_count+=1        
        #increment time
        time = time + dt

    #close output file for axial strain data
    output.close()   
 
    return(DIC_time,DIC_strain)

##############################################################################################
# This function will compute the average surface strain from Vic3D as a function
#    of length along the load axis to check uniformity of strain
###############################################################################################
def calculate_DIC_axial_strain_profiles(expt,DIC_first,DIC_last,DIC_frame_rate,gage_l,DIC_col,DIC_strain):
    #find initial gage length from first DIC data file
    file_root='T'+str(int(expt))+'_0_'
    # read in first correlation data set to determine initial correlation (gage) length
    filenumber=str(int(DIC_first)).zfill(6)
    filename=file_root + filenumber + '.csv' 
    initial_DIC_datafile = open(filename,"r")
    initial_x, initial_y=numpy.loadtxt(initial_DIC_datafile,delimiter=",",skiprows=1,usecols=(0,1),unpack=True)
    initial_x_array=numpy.array(initial_x) # used to calculate DIC subset spacing    
    initial_y_array=numpy.array(initial_y) #must cast arrays as numpy arrays before doing arithmetic
    y_max=initial_y_array.max()
    y_min=initial_y_array.min()
    initial_gage=y_max-y_min
    print(initial_gage)
    
    #estimate DIC subset spacing
    subset_spacing=initial_x_array[1]-initial_x_array[0]
    subset_spacing=abs(subset_spacing)
    print('subset spacing:',subset_spacing)    
    initial_DIC_datafile.close()
    
    #define interframe time
    frame_rate=float(DIC_frame_rate)
    dt = 1/frame_rate
    
    #define desired gage length in mm and other data file attributes needed
    desired_gage=gage_l*1000.0
    first_file_number=int(DIC_first)
    last_file_number=int(DIC_last)
    data_col_to_avg=int(DIC_col)

    #define y limits to average within
    y_max_new = desired_gage
    y_max_new = y_max_new/2.0

    #open output file
    outfilename=file_root + '_axial_displ.csv'
    output=open(outfilename,"w")
    output.write('initial gage: ')
    out1=str(initial_gage)
    output.write(out1)
    output.write('  desired gage: ')
    out1=str(desired_gage)
    output.write(out1)
    output.write('\n')
    output.write('filenumber, time (microsec), row_count, y_avg (mm), eyy avg, std')
    output.write('\n')

    # loop through specific number of DIC data files and average the column selected
    # over the specified gage length
    time = 0.0
    ####################################################
    #
    #   delta_y is the value used to identify a new row of DIC points along the length of the specimen. x is usually 
    #       transverse to the axis, and y is along the axis. Data averaging occurs over x for constant values of y,
    #       so we need to find when the end of a row of x's are by identifying a subset-level jump in the value of y
    #       This number is not exact due to slight tilting of the data. Should probably rotate the data to account for this but 
    #       the angle is usually small
    #
    ##################################################
    delta_y = 0.9*subset_spacing 
    for i in range(first_file_number,last_file_number):
        filenumber=str(i).zfill(6)
        filename=file_root + filenumber + '.csv' 
        filestream = open(filename,"r")
        y_index, var_to_avg=numpy.loadtxt(filestream,delimiter=",",skiprows=1,usecols=(1,data_col_to_avg),unpack=True)
        
        # row index keeps track of what row we are on
        row_index = 1
        # this structure is supposed to reduce the total DIC data into an array for averaging
        # by using a boolean array that has in it index values that satisfy the inequality   
        j = numpy.array(len(y_index)*[False],bool)  
        row_count=1
        count=0
        flag_first=1
        for k in y_index:
            #print("count:",count)
            if(abs(k)<y_max_new):                  
                if(flag_first==0):
                    #print("row_count",row_count)                                                
                    if(row_count>1):
                        delta_y_actual = abs(k - y_index[count-1])
                        #print("count:",count,"k:",k,"delta_y_actual:",delta_y_actual)
                        if(abs(delta_y_actual)<delta_y):
                            #print("In a row",row_count)                    
                            j[count]=True
                            row_count=row_count+1                        
                        if(abs(delta_y_actual)>delta_y):
                           # print("end of row is reached", row_count,"time:",time)                    
                                #end of row is reached, average current row data and reset subselector and row counter                                        
                            # print averaged eyy to data file
                            a=numpy.array(var_to_avg[j])
                            avg_var=numpy.average(a)
                            std_var=numpy.std(a)
                            b=numpy.array(y_index[j])
                            avg_y=numpy.average(b)
                            std_y=numpy.std(b)
                            filenum_str=str(filenumber)
                            time_str=str(time)
                            rowcount_str=str(row_count)
                            avg_var_str=str(avg_var)
                            std_var_str=str(std_var)
                            avg_y_str=str(avg_y)
                            output.write(filenum_str)
                            output.write(',')
                            output.write(time_str)
                            output.write(',')
                            output.write(rowcount_str)
                            output.write(',')
                            output.write(avg_y_str)
                            output.write(',')   
                            output.write(avg_var_str)
                            output.write(',')
                            output.write(std_var_str)
                            output.write('\n')
                            # reset the subarray selector to True for the current element, which is the first element of the next row to average
                            j = numpy.array(len(y_index)*[False],bool)
                            j[count]=True
                            #reset row count to 2 since we are starting the next row and the first element has already been flagged above
                            row_count=2
                if(flag_first==1):
                    #print("Flag First Triggered")                
                    j[count]=True         
                    flag_first=0
                    row_count=2
                  # print("In flag triggered IF, count = ",count,"flag_first:",flag_first)
            count=count+1
        #write last data point in the file/time point   
        a=numpy.array(var_to_avg[j])
        avg_var=numpy.average(a)
        std_var=numpy.std(a)
        b=numpy.array(y_index[j])
        avg_y=numpy.average(b)
        std_y=numpy.std(b)
        filenum_str=str(filenumber)
        time_str=str(time)
        rowcount_str=str(row_count)
        avg_var_str=str(avg_var)
        std_var_str=str(std_var)
        avg_y_str=str(avg_y)
        output.write(filenum_str)
        output.write(',')
        output.write(time_str)
        output.write(',')
        output.write(rowcount_str)
        output.write(',')
        output.write(avg_y_str)
        output.write(',')   
        output.write(avg_var_str)
        output.write(',')
        output.write(std_var_str)
        output.write('\n')   
        #increment time
        time = time + dt
    #close output file for axial strain data
    output.close()
    
    # estiamte time points for automated plot of axial strain distributions at different strain levels
    
    # first pick number of strain levels to plot
    num_of_strains_to_plot=5
    DIC_strain=numpy.array(DIC_strain)
    peak_strain_index=numpy.argmax(DIC_strain) # this is the frame number for the peak strain
    peak_strain_frame=int(peak_strain_index+DIC_first) # this is the frame number for the peak strain
    peak_strain=DIC_strain[peak_strain_index]
    print('peak DIC strain:',peak_strain)
    print('peak strain index',peak_strain_index)
    print('peak strain frame',peak_strain_frame)
    delta_strain=peak_strain/float(num_of_strains_to_plot)  
    strain_per_frame=1000.0/frame_rate
    print('delta strain:',delta_strain)
    print('strain per frame:',strain_per_frame)
    indexes_to_plot=numpy.arange(0,num_of_strains_to_plot,1)
    
    filestream = open(outfilename,"r")
    frame,max_r,avg_y,avg_eyy=numpy.loadtxt(filestream,delimiter=",",skiprows=2,usecols=(0,2,3,4),unpack=True) # read file that was just written to get averaged axial data

    maxrows=int(initial_gage/subset_spacing) # this is the max number of y points expected
    print('maxrows:',maxrows)
    maxrows=float(maxrows)
    
    
    plt.figure()  #hopefully opens a new figure for plotting axial strain data
    index=0  # this counts the number of plots                 
    for x in indexes_to_plot:
        x_plot_data=numpy.arange(0.0,maxrows+1.0,1.0)
        y_plot_data=numpy.arange(0.0,maxrows+1.0,1.0)
        x_plot_data=numpy.array(x_plot_data)
        y_plot_data=numpy.array(y_plot_data)
        data_count=0 # this is a counter for the data for plotting
        # first calculte the framenum for plotting
        frame_num=peak_strain_frame-index*int(delta_strain/strain_per_frame)
        print('frame num:',frame_num)
        i=0 # this is the index in the data file        
        for j in frame:
            if j == frame_num:
                x_plot_data[data_count]=avg_y[i]
                y_plot_data[data_count]=avg_eyy[i]
                data_count+=1
            i+=1
            # slice data arrays to not plot zeros
            #print('data count:',data_count)
            x_plot_data2=x_plot_data[:data_count]
            y_plot_data2=y_plot_data[:data_count]
            plt.plot(x_plot_data2,y_plot_data2,'o')
        index+=1
        #print('x plot data 2')
        #print(x_plot_data2)
    plt.show()
    
    return()


###########################################################################################################################
#
#   Function to identify wave start times for Kbar analysis
#
############################################################################################################################
def find_data_pulse4(time_data,data_data,window_size=400,pulse_to_find=1):
    #this method looks for a sustained change in the level over some window
    # pulse to find is 1=incident, 2=reflected and 3=transmitted, 4=DIC pulse
    # incident pulse has a larger requirement for sustainled level (check level) due to strain gage noise problem
    # reflected pulse should be tensile (positive pulse level)
    # incident and transmitted pulses should both have a negative pulse level
    data=numpy.array(data_data)
    time=numpy.array(time_data)
    #print('length of numpy array')    
    #print(len(data))
    dt=time_data[1]-time_data[0]
   # print('dt = ',dt)
    #delay=int(0.000025/dt) # amount of time ahead of a pulse to include in the pulse
    delay=int(0.0001/dt) # amount of time ahead of a pulse to include in the pulse
    delay_inc=0.000000000
        
    
    checksum_min= 150  # number of time points a pulse must be above a certain level to be considered a pulse
    avg_level=numpy.average(data[0:delay])
    std_level=numpy.std(data[0:delay])

# different check levels are used based on the typical size of the various pulses
    if pulse_to_find == 1:
        check_level=0.0005
    if pulse_to_find == 2:
        check_level=0.0005
    if pulse_to_find == 3:    
        check_level=0.0001  # this is a check level in absolute volts to capture a pulse
    if pulse_to_find == 4:    
        check_level=0.001  # checks level in strain to capture DIC strain-time pulse
        delay=100
    #slope_data=numpy.array(data_data)
   # print('sign:',sign)
   # print('std level =',std_level)
    #calculate the expected value of start time

    increment=0
    checksum=0
    flag = 0
    
    #provisional start time is the first data point in time-data
    provisional_start_time=time_data[0]    
    start_time=provisional_start_time
    
    for i in data:  
        if increment>delay:
            level_check = (data[increment]-avg_level)
            level_check=abs(level_check)
            slope = data[increment]-data[increment-1]
            slope = abs(slope/dt)
            if pulse_to_find == 4:            
                print ('level check',level_check,check_level)
            if slope < 50: #this determines how many points are above the check level criteria
                if level_check > check_level:
                    if flag == 0:        
                        start = increment
                        flag = 1
                        if pulse_to_find == 4:            
                            print ('flag set to 1')
                       
                    checksum+=1
        increment+=1

    if checksum > checksum_min:     
        if pulse_to_find == 1: #incident pulse
            start_time = start*dt
            time_window = time_data[start-delay_inc:start-delay_inc+window_size]
            data_window = data_data[start-delay_inc:start-delay_inc+window_size]
        if pulse_to_find == 3: #trans pulse
            start_time = start*dt
            time_window = time_data[start-delay:start-delay+window_size]
            data_window = data_data[start-delay:start-delay+window_size]
        if pulse_to_find == 4: #DIC pulse
            start_time = start*dt
            time_window = time_data[start-delay:start-delay+window_size]
            data_window = data_data[start-delay:start-delay+window_size]



  #      try: 
   #         start_time = start*dt
   #     except:
   #         start = delay
  #          start_time = 0.0
  #          print('***wave start not found***')
        
   #     try:    
   #         time_window = time_data[start-delay:start-delay+window_size]
   #         data_window = data_data[start-delay:start-delay+window_size]
   #     except:
   #         time_window = 0
   #         data_window = 0
   #         print('***data window not found***')
     
    if start_time <= provisional_start_time:
        print('**************START TIME NOT FOUND************')
        print('pulse to find:',pulse_to_find)
        start_time=provisional_start_time
        start=0
        time_window = time_data[start:start+window_size]
        data_window = data_data[start:start+window_size]
            
     
    if pulse_to_find == 2: # this is a reflected pulse - the code assumes the first (incident) pulse has been found accurately
        #refl_delay_shift=0.0003/dt
        #refl_delay_shift=int(refl_delay_shift) # for use when trying to pick out reflected wave
        #refl_delay_shift=window_size+int(0.0003/(2.0*dt))
        #print('refl_delay_shift:',refl_delay_shift)
        #provisional start time is the first data point in time-data
       # print('window size:',window_size)
        #print('start:',start)

        increment = start+window_size# this jumps ahead in the data to find the reflected pulse
        #print('increment:',increment)
        provisional_start_time=time_data[increment]    
        start_time=provisional_start_time
        data_length=len(data)
        
        #this is the data subset that I am further subsampling based on the incident wave start found earlier
        #data_subset2=data[start:data_length]
        #time_subset2=time[start:data_length]

        data_wave1=data[start:start+window_size]
        time_wave1=time[start:start+window_size]


        data_subset=data[increment:data_length]
        time_subset=time[increment:data_length]
        new_avg_level=numpy.average(data_subset[0:delay])

        # plot part that I am searching for the reflected wave
        #pylab.clf()    # this is needed to clear data from previous plots
        #pylab.plot(time_subset,data_subset, color = 'r', label='subset sg1')
       # pylab.plot(time_wave1,data_wave1, color = 'b', label='wave1 sg1')
        #pylab.legend()
       # pylab.show()
        
        
        checksum=0
        flag=0
        increment2=0 #redefine increment for searching the data subset corresponding to the reflected wave
        time_shift=increment-increment2 # this allows me to pick out the start time in the orignal sg1 data set
        for i in data_subset:
            if increment2 > delay:              
                level_check = data_subset[increment2]-new_avg_level 
                level_check=abs(level_check)
                slope = data_subset[increment2]-data_subset[increment2-1]
                slope = abs(slope/dt)
                #print('level_check:',level_check)
                #print ('reflected slope:',slope)
                if slope < 50: #this eliminates counting spikes as pulse points
                    if level_check > check_level:
                        if flag == 0:        
                            start = increment2
                            flag = 1
                        checksum+=1
                        #print('checksum:',checksum)
            increment2+=1

        if checksum > checksum_min:     
            try: 
                start_time_refl = start*dt
                start_time=(start+time_shift)*dt
                #print('ive calculated a start time for refl wave:',start,start_time_refl,start_time)
            except:
                start = delay
                start_time = 0.0
                print('***wave start not found***')
            
            try:    
                time_window = time_subset[start-delay:start-delay+window_size]
                data_window = data_subset[start-delay:start-delay+window_size]
            except:
                time_window = 0
                data_window = 0
                print('***data window not found***')


     #   if start_time <= provisional_start_time:
      #      print('**************START TIME NOT FOUND************')
       #     print('pulse to find:',pulse_to_find)
        #    start_time=provisional_start_time
         #   start=0
          #  time_window = time_data[start:start+window_size]
           # data_window = data_data[start:start+window_size]
         
    
        
    ################ skip the first pulse and check for a second if the reflected pulse is being sought        
        
        
        
   # print('start time [s]:', start_time)
    return(start_time,time_window,data_window)


#######################################################################################
#
#   Preview KBar raw data text file in a window
#
#######################################################################################

def preview_KBar_data_file(filename):
    #function opens Kbar raw text data file for quick previewing to double check column orientation
    preview_characters = 1000  # defines number of characters to preview in input file
    f=open(filename,"r")    
    preview_text=f.read(preview_characters)
    # this stuff opens up a window for showing the text
    root=Tk()
    root.title('preview of Kbar Data file')
    geom = "1100x700"
    root.geometry(geom)
    text=Text(root)
    text.insert(INSERT,preview_text)
#    text.pack(cnf=expand)
    # need to figure a way to change the display width
    text.pack()
    root.mainloop()


def load_timeseries_data(filename,skiprows=8,col1=0,col2=1):
    #this function loads a time series data file, according to columns passed to it, and plots it if desired  
    time_data,volts_data=numpy.loadtxt(filename,skiprows=skiprows,usecols=(col1,col2),unpack=True)
    #pylab.clf()    
    #pylab.plot(time_data, volts_data, color = 'b', label='volts')
    #pylab.title("Raw KBar Data")
    #pylab.xlabel("Time [s]")
    #pylab.ylabel("Volts")
    #pylab.show()
    return(time_data,volts_data)
    
def zero_timeseries_data(volts_data,avg_window=100):
    baseline_start=0
    baseline_end=avg_window
    data=numpy.array(volts_data)
    baseline_avg = numpy.average(data[baseline_start:baseline_end])
    #print('avg volts: ',baseline_avg)
    #data=data.__sub__(baseline_avg)
    data=data-baseline_avg    
    return(data)

def quickplot_xydata(x_data,y_data):
    pylab.clf()  # this is needed to clear data from previous plots
    pylab.plot(x_data, y_data, color = 'b')
    pylab.title("Quickplot")
    pylab.xlabel("X Data")
    pylab.ylabel("Y Data")
    pylab.show()
    
def volts_to_strain_data(volts_data,cal_factor=1):
    data=numpy.array(volts_data)
    strain_data=data*cal_factor
    return(strain_data)
    
def volts_to_TC_temp(volts_data,type):
    data=numpy.array(volts_data)
    
    #this stuff was added to accomodate reading in analysis params from a file, which defaults to float input
    type=int(type)
    
    if type == 1:
        type = 'k'
    if type == 0:
        type = 'none'
    if type == 2:
        type='R'       

    
    if type == 'k':    
        # k-type thermocouple constants
        T0=725.0
        V0=30.17384
        p1=23.985623
        p2=-0.33624743
        p3=0.010218015
        p4=-0.000819765
        q1=-0.016475689
        q2=0.00040158272
        q3=-0.0000317366
    elif type == 'R':     
        # R-type thermocouple constants 760 C - 1275 C
        T0=1038.213
        V0=11.01476
        p1=74.66934
        p2=3.40907
        p3=-0.145112
        p4=0.0063077
        q1=0.05688
        q2=-0.0020513
        q3=0.000000
    elif type == 'none':    # just put values to avoid undefined variable useage error 
        T0=1.0
        V0=1.0
        p1=1.0
        p2=1.0
        p3=1.0
        p4=1.0
        q1=1.0
        q2=1.0
        q3=1.0
    elif type == 'bad':     # just put values to avoid undefined variable useage error
        T0=1.0
        V0=1.0
        p1=1.0
        p2=1.0
        p3=1.0
        p4=1.0
        q1=1.0
        q2=1.0
        q3=1.0
    else:
        print('Thermocouple Type Error')
    
    if type == 'none':    
        TC_data=23.0  #room temperature when test is done without heating
    elif type == 'bad':
        TC_data=0.0 #return a dummy temperature when TC signal is known to be bad
    else:
        numerator=((data-V0)*(p1+(data-V0)*(p2+(data-V0)*(p3+(data-V0)*p4))))
        denomenator=(1+(data-V0)*(q1+(data-V0)*(q2+(data-V0)*q3)))
        TC_data=T0+numerator/denomenator

    return(TC_data)

def volts_to_pyro_temp(volts_data,sensitivity=1):

    volts_data=numpy.array(volts_data)
    #print(sensitivity)
    if sensitivity == 1:
        #pyro constants 1uA/v
        pyro_a=779.69175
        pyro_b=-58.93766
        pyro_c=0.01115
    elif sensitivity == 10:
        #pyro constants 10uA/v
        pyro_a=954.13502
        pyro_b=-75.29034
        pyro_c=0.00402
    elif sensitivity == 50:
        #pyro constants 50uA/v
        pyro_a=1118.85429
        pyro_b=-101.56280
        pyro_c=0.00840
    elif sensitivity == 0: #case where the test is done below NIMPY sensitivity
        pyro_a=1.0
        pyro_b=1.0
        pyro_c=1.0
    else:
        print('NIMPY Sensitivity Error!!!')
        pyro_a=1.0
        pyro_b=1.0
        pyro_c=1.0

    nimpy_data = pyro_a-pyro_b*numpy.log(volts_data+pyro_c)-273.17
    return(nimpy_data)


def pick_stress_given_strain(filename,strain_for_plot):
    # idea is to pick stress from stress-strain curve at a given strain point for plotting thermal softening
    #open file with stress-strain data and find value of stress corresponding to input strain value
    data_input = open(filename,"r")
    #datafile is the file written by python stress-strain curve calculator
    strain_data,stress_data=numpy.loadtxt(data_input,skiprows=0,delimiter=',',usecols=(0,1),unpack=True)
    strain_data=numpy.array(strain_data)   
    stress_data=numpy.array(stress_data)    

    tol=0.005 # defines how close we need to get to the selected strain
    avg_window=5 # defines the number of points used for the stress average    
    increment=0
    flag = 0
    for i in strain_data:       
        if strain_data[increment] > 0.02:            
            try:
                finder=abs(strain_data[increment]-strain_for_plot)
            except:
                finder=0.0            
            if finder < tol:
                flag = increment
        increment+=1

    if flag == 0:
        print('Strain Selected is Out of Bounds!')
        stress_for_plot=0.0
    else:
        stress_for_plot = numpy.average(stress_data[flag-avg_window:flag+avg_window])  
        #stress_for_plot = stress_data[flag]
    data_input.close()
#    print(stress_for_plot)        
    return(strain_for_plot,stress_for_plot)

#def test_temperature_calculator(test,time_col=0,sg1_col=1,TC_col=2,N_col=3,N2_col=4,TC_type='k',N_sens=1):
def test_temperature_calculator(test,time_col=0,TC_type='k',N_sens=1):

    #DEFINE AVG EMISSIVITY FOR GROUP OF TESTS BASED ON PREVIOUS DATA FOR CASE WHERE TC IS BAD
    #avg_emissivity=0.3    
    global average_emissivity    
    avg_emissivity=average_emissivity
    
    test=int(test)
    test=str(test)

    sg1_col,sg2_col,TC_col,N_col,N2_col=find_data_columns(test)
    sg1_col=int(sg1_col)
    sg2_col=int(sg2_col)
    TC_col=int(TC_col)
    N_col=int(N_col)
    N2_col=int(N2_col)
    
    input_file=test+'.txt'
    data_input = open(input_file,"r")
    #time_data,N2_data,sg1_data,TC_data,N_data=numpy.loadtxt(data_input,skiprows=8,usecols=(0,1,2,3,4),unpack=True)
    time_data,sg1_data,TC_data,N_data,N2_data=numpy.loadtxt(data_input,skiprows=8,usecols=(time_col,sg1_col,TC_col,N_col,N2_col),unpack=True)
    dt=time_data[1]-time_data[0]
    #print(dt)

    time_array=numpy.array(time_data)
    TC_array=numpy.array(TC_data)
    N_array=numpy.array(N_data)
    N2_array=numpy.array(N2_data)



    #calculate the expected value of start time
    increment=0
    flag = 0
    for i in sg1_data:  
        if increment>3:
            slope = sg1_data[increment]-sg1_data[increment-3]
            slope = -slope/(3*dt)
            #print(slope)
            if slope > 750:
                if flag == 0:        
                    start = increment
                    flag = 1
        increment+=1
    start_time = start*dt
    # print('start time [s]:', start_time)
    
    # input the delay and average window in terms of number of time increments
    delay = 0
    window = 250

    #average the N, N2 and TC voltages in the window defined by start+delay to start+delay+window
    avg_N = numpy.average(N_array[start+delay:start+delay+window])
    avg_N2 = numpy.average(N2_array[start+delay:start+delay+window])
    avg_TC = numpy.average(TC_array[start+delay:start+delay+window])*1000.0
    
    N_subarray=N_array[start+delay:start+delay+window]
    time_subarray=time_array[start+delay:start+delay+window]
    
    #Diagnostic prints
    #print('avg N [V]:', avg_N)
    #print('avg N2 [V]:', avg_N2)
    #print('avg TC [mV]:',avg_TC)
    #pylab.plot( time_subarray, N_subarray, label='Data for Emissivity', linestyle='dotted', linewidth = 5 )
    #pylab.show()

    # define pyro constants    
    c2=0.014388
    lamda=0.0000012
   
 # calculate temperatures  from voltages depending on the case structure    
      
    if TC_type == 'bad': # assume NIMPY is good (N_sense not 0) if TC is bad
        emissivity = avg_emissivity
        avg_N_temp=volts_to_pyro_temp(avg_N,N_sens)
        avg_N2_temp=volts_to_pyro_temp(avg_N2,N_sens)
        avg_TC_temp=0.0
        N_trueT=1/((1/avg_N_temp)+lamda/c2*numpy.log(emissivity))
        N2_trueT=1/((1/avg_N2_temp)+lamda/c2*numpy.log(emissivity))
        avg_test_temp=(N_trueT+N2_trueT)/2.0
        temp_unc=abs(N_trueT-N2_trueT)
    elif TC_type == 'none': #do the same thing as if the TC is bad
        emissivity = avg_emissivity
        avg_N_temp=volts_to_pyro_temp(avg_N,N_sens)
        avg_N2_temp=volts_to_pyro_temp(avg_N2,N_sens)
        avg_TC_temp=0.0
        N_trueT=1/((1/avg_N_temp)+lamda/c2*numpy.log(emissivity))
        N2_trueT=1/((1/avg_N2_temp)+lamda/c2*numpy.log(emissivity))
        avg_test_temp=(N_trueT+N2_trueT)/2.0
        temp_unc=abs(N_trueT-N2_trueT)
    else: # assume TC, NIMPY and NIMPY 2 are present and correct
        avg_TC_temp=volts_to_TC_temp(avg_TC,TC_type)
        avg_N_temp=volts_to_pyro_temp(avg_N,N_sens)
        avg_N2_temp=volts_to_pyro_temp(avg_N2,N_sens) 
        if avg_TC_temp >= avg_N_temp:
            emissivity=numpy.exp(c2/lamda*(1/avg_TC_temp-1/avg_N_temp))
        if avg_TC_temp < avg_N_temp:
            print('TC Lower than NIMPY!')
            emissivity=avg_emissivity
        N_trueT=1/((1/avg_N_temp)+lamda/c2*numpy.log(emissivity))
        N2_trueT=1/((1/avg_N2_temp)+lamda/c2*numpy.log(emissivity))
        avg_test_temp=(N_trueT+N2_trueT)/2.0
        temp_unc=abs(N_trueT-N2_trueT)

# Catch for case where there are no NIMPY's. Not efficient, needs to be better        
    if N_sens == 0: #assume TC is good
        if TC_type == 'none':
            avg_test_temp=23.0 #room temperature
            temp_unc=3.0 #uncertainty in room temperature
            print('no good temperature available for this test!')
        elif TC_type == 'bad':
            avg_test_temp=23.0 #room temperature
            temp_unc=3.0 #uncertainty in room temperature
            print('no good temperature available for this test!')
        else:
            avg_TC_temp=volts_to_TC_temp(avg_TC,TC_type)
            avg_test_temp=avg_TC_temp
            temp_unc=3.0

       
       
    data_input.close()
    #Diagnostic prints
    #print('avg N impact [C]:', avg_N_temp)
    #print('avg N2 impact [C]:', avg_N2_temp)
    #print('avg TC impact [C]:',avg_TC_temp)
    #print('emissivity:',emissivity)
    #print('avg_test_temp :',avg_test_temp)
    #print('temp uncertainty:',temp_unc)


    return(avg_test_temp,temp_unc)

#------------------END OF IMPACT TEMPERATURE ESTIMATOR-------------------------------

    
      
    

#----------------------------------------------------------------------------       
def write_text_file(filename,array1,array2,array3):

    #parts = filename.split('.')
    #outfile = parts[0]+'_analyzed.txt'
    data_outfile = open(filename,"w")

    data_to_write=numpy.column_stack((array1,array2,array3))
    numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
    data_outfile.close()  
#-----------------------------------------------------------------



########################################################################################
#
#   Function to calculate and plot Kolsky bar test temperatures along with flow stresses
#
########################################################################################

def calculate_and_plot_temperatures(expt_list,TC_type,N_sens):

    #define a default emissivity
    average_emissivity=0.36  # put this in for cases where TC is no good
    
    total_data_points=len(expt_list)
    avg_T=numpy.arange(0,total_data_points,dtype=numpy.float)
    T_unc=numpy.arange(0,total_data_points,dtype=numpy.float)
    stress_for_plot=numpy.arange(0,total_data_points,dtype=numpy.float)
    stress_error=numpy.arange(0,total_data_points,dtype=numpy.float)
    
    
    
    #write data to output file
    filename='flowstress_vs_T.txt'
    data_outfile = open(filename,"w")
    
    #show a set of plots with stress-vs-temperature at a given values of strain
    strain_for_plot = [0.05,0.1,0.2,0.3] # value of total strain for plotting thermal softening
    Title_1='Strain = '+str(strain_for_plot[0]) 
    Title_2='Strain = '+str(strain_for_plot[1]) 
    Title_3='Strain = '+str(strain_for_plot[2]) 
    Title_4='Strain = '+str(strain_for_plot[3])
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].set_title(Title_1)
    axarr[0, 1].set_title(Title_2)
    axarr[1, 0].set_title(Title_3)
    axarr[1, 1].set_title(Title_4)
    
    count=0
    
    for strain in strain_for_plot:
        increment=0
        for expt in expt_list:
            expt=int(expt)        
            expt=str(expt)        
            ss_filename=expt + '_analyzed.txt'
            #print('am i here?',increment)
            #avg_T[increment],T_unc[increment]=test_temperature_calculator(expt,time_col[increment],sg1_col[increment],TC_col[increment],N_col[increment],N2_col[increment],TC_type[increment],N_sens[increment])
            avg_T[increment],T_unc[increment]=test_temperature_calculator(expt,0,TC_type[increment],N_sens[increment])
            print(expt,avg_T[increment],T_unc[increment]) 
            #print('i am here # ',increment)
            if T_unc[increment] <= 5.0 :
                T_unc[increment] = 5.0
            strain,stress_for_plot[increment]=pick_stress_given_strain(ss_filename,strain)
            stress_error[increment]=50.0 #need to calculate this on an experiment-by-experiment basis
            increment+=1
    #now plot the data at specified strain locations
    # write set of data to output file 
    
        data_outfile.write('Fixed Strain Value:')
        data_outfile.write(str(strain))
        data_outfile.write('\n')
        
        data_outfile.write('Expt#,AvgT,T Unc,Stress,Stress Unc')
        data_outfile.write('\n')
        
        array1=[int(x) for x in expt_list]
        array2=avg_T
        array3=T_unc
        array4=stress_for_plot
        array5=stress_error
        data_to_write=numpy.column_stack((array1,array2,array3,array4,array5))
        numpy.savetxt(data_outfile,data_to_write,delimiter=',',fmt='%-10.5e')
        data_outfile.write('\n')
    
    
        # now plot data on one of four plots, one for each strain value
        if count==0:
            axarr[0, 0].scatter(avg_T, stress_for_plot)
        if count==1:
            axarr[0, 1].scatter(avg_T, stress_for_plot)
        if count==2:
            axarr[1, 0].scatter(avg_T,stress_for_plot)
        if count==3:
            axarr[1, 1].scatter(avg_T, stress_for_plot) 
        count+=1
    # show flow stress vs temperature plot for all for selected strains 
    plt.show()





#-----------------------------------------------------------------    
