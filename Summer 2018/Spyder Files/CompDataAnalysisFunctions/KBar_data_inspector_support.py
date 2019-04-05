# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 15:52:07 2016

@author: smates
"""
#from Tkinter import *
import numpy
import pylab

global average_emissivity

#pylab.clf()  # this is needed to clear data from previous plots

#expt='3563' 
#preview_KBar_data_file('3389.txt')

#inc_start,refl_start,trans_start=calculate_wave_timing('3411',60.0,0.0,60.0,0.0,400)
#quickplot_volts_data(expt)
#quick_stress_strain('3793')

#exp_list=['3262','3263','3264']
#low_strain=0.03
#hi_strain=0.14
#average_stress_strain_data(exp_list,low_strain,hi_strain)

#quickplot_volts_data('1001')

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
    r_delay=0.000340
    t_delay=0.000346
    mod=170.0
    gfoil=1
    
    #print('i am doing nothing in this function')
    strain,stress,strain_rate,opt_refl_delay,opt_trans_delay,equilib=calculate_stress_strain(expt,dia,thk,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window,r_delay,t_delay,gfoil,mod)
     
    quickplot_xydata(strain,stress)
    return()


########################################################################################
#
#   Function to read a text file containing experiment files and analysis parameters from a single text file
#
########################################################################################

def read_expts_for_analysis(analysis_file):

    infile = open(analysis_file,"r")
   
    expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay=numpy.loadtxt(infile,skiprows=1,delimiter=',',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)
    #return(number,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay)
    #number,sg1_hi=numpy.loadtxt(infile,skiprows=1,usecols=(0,1),unpack=True)
    return(expt_list,sg1_hi,sg1_lo,sg2_hi,sg2_lo,dia,thk,gfoil,TC_type,N_sens,bar_modulus,window,refl_delay,trans_delay)

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
    
    
    sg1_volts=zero_timeseries_data(sg1_volts,avg_window=100)
    sg2_volts=zero_timeseries_data(sg2_volts,avg_window=100)

    pylab.clf()    # this is needed to clear data from previous plots
    pylab.plot(sg1_time, sg1_volts, color = 'b', label='SG1')
    pylab.plot(sg2_time, sg2_volts, color = 'r', label='SG2')
    pylab.legend()
    pylab.show()

    return()



########################################################################################
#
#   Function to identify the data columns in a Kolsky bar text file
#
########################################################################################
def find_data_columns(expt):
    filename= expt + '.txt'
    #filename='3411.txt'
    filestream=open(filename,'r')
    #read the line in the data file that contains column headings (line 6 usually)
    i=0        
    for lines in filestream:
        #print(lines)    
        if i == 4:
            heading=lines        
            #print('i am here')
        elif i > 10:
            break
        i+=1
    filestream.close()
    
    col_titles=heading.split("\t")
    #print('col_title:',col_titles)
    sg1_col=col_titles.index("SG1")   
    sg2_col=col_titles.index("SG2")
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
        flag=0
    try:
        N_col=col_titles.index("NIMPY")
    except:
        flag=flag+1
    if flag==1:
        try:
            N_col=col_titles.index("Pyro")
        except:
            N_col=0
            print('No NIMPY column in KBar Data File!')
            flag=0
    try:
        N2_col=col_titles.index("NIMPY 2")
    except:
        N2_col=0
        print('No NIMPY 2 column in KBar Data File!')

    #print('SG1_col:',sg1_col)

    #return()
    return(sg1_col,sg2_col,TC_col,N_col,N2_col)

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

  #  pylab.clf()    # this is needed to clear data from previous plots
  #  pylab.plot(inc_time, inc, color = 'b', label='INC')
  #  pylab.plot(inc_time, refl, color = 'r', label='REFL')
  #  pylab.plot(inc_time, trans, color = 'y', label='TRANS')
  #  pylab.legend()
  #  pylab.show()
    
    #print(inc_time)
    
    return(inc_start,refl_start,trans_start)


########################################################################################
#
#   Function to perform KBar Mechanical Analysis and deliver stress-strain
#
########################################################################################
def calculate_stress_strain(test,dia,thk,sg1_hi,sg1_lo,sg2_hi,sg2_lo,window=600,refl_delay=0.000340,trans_delay=0.000346,gfoil=1,bar_modulus=170):

    # select the test number
    #test='3190'
    #enter sample size, also will be in database
    
    sg1_col,sg2_col,TC_col,N_col,N2_col=find_data_columns(test)
    #print('sg2_col',sg2_col)
    
    #test=int(test)
    #test=str(test)
    hi_rate_data=test + '.txt'
    output_file=test+'_analyzed.csv'
    #
    #sg1_col=int(sg1_col)
    #sg2_col=int(sg2_col)
    sg1_time,sg1_volts=load_timeseries_data(hi_rate_data,8,0,sg1_col)
    sg2_time,sg2_volts=load_timeseries_data(hi_rate_data,8,0,sg2_col)
    sg1_volts=zero_timeseries_data(sg1_volts,100)# 100 is the number of data points to average
    sg2_volts=zero_timeseries_data(sg2_volts,100)# 100 is the number of data points to average  
            
    sg1_cal_volts=2*(sg1_hi-sg1_lo)/1000.0
    sg2_cal_volts=2*(sg2_hi-sg2_lo)/1000.0
    
    sg1_cal=-0.004854/sg1_cal_volts
    sg2_cal=-0.004854/sg2_cal_volts    
    sg1_strain=volts_to_strain_data(sg1_volts,sg1_cal)
    sg2_strain=volts_to_strain_data(sg2_volts,sg2_cal)
    
    #normalize sg1, sg2 time to zero on sg1 time
    
    sg1_time=sg1_time-sg1_time[0]
    sg2_time=sg1_time-sg1_time[0]
    
    #quickplot_xydata(sg2_time,sg1_volts)
    #Identify data pulses*********************************************
    inc_start,inc_time,inc=find_data_pulse4(sg1_time,sg1_strain,window,1) #last number is 1=incident wave, 0=refl or trans
    #quickplot_xydata(inc_time,inc)
    
    refl_start,refl_time,refl=find_data_pulse4(sg1_time,sg1_strain,window,2)
    #quickplot_xydata(refl_time,refl)
    
    trans_start,trans_time,trans=find_data_pulse4(sg2_time,sg2_strain,window,3)
    #quickplot_xydata(trans_time,trans)
    
#    pylab.clf()    
#    pylab.plot(inc_time, inc, color = 'y', label='inc')
 #   pylab.plot(inc_time, refl, color = 'b', label='refl')
 #   pylab.plot(inc_time, trans, color = 'b', label='trans')
 #   plt.legend(loc='lower right')
 #   pylab.show()
    
   # print('before calculating')
   # print('refl_delay:',refl_delay)
   # print('trans_delay:',trans_delay)

    #####################################################################    
    # ****Enable these to have the program identify wave start times****
    ######################################################################    
    #refl_delay=refl_start-inc_start
    #trans_delay=trans_start-inc_start
    print('after calculating')
    print('inc start[s]:',inc_start)    
    #print('refl start[s]:',refl_start)    
    #print('trans start[s]:',trans_start)    
    print('refl_delay:',refl_delay)
    print('trans_delay:',trans_delay)
    
    ######### check to see if trans delay is more than a small amount bigger than incident delay. If it is, we are capturing the second pulse
#    delay_check = abs(trans_delay-refl_delay)
#    if delay_check > 0.5*refl_delay:
#        print('delay check triggered')
#        trans_start,trans_time,trans=find_data_pulse3(sg2_time,sg2_strain,window,1,0) #try to find trans wave again, flipping polarity sign
#        new_trans_delay=trans_start-inc_start        
#        print('new trans delay:',new_trans_delay)
    
#    trans_delay=trans_start-inc_start

    # Perform Kolsky analysis to get stress-strain curves*************************************
    strain,stress,strain_rate,opt_refl_delay,opt_trans_delay,equilib=analyze_Kolsky_pulses(test,sg1_time,sg1_strain,sg2_strain,inc_start,refl_delay,trans_delay,window,dia,thk,gfoil,bar_modulus)

    print('after analyzing pulses')    
    print('opt_refl_delay:',opt_refl_delay)
    print('opt_trans_delay:',opt_trans_delay)


    write_text_file(output_file,strain,stress,strain_rate)
    #plt.plot(strain,stress)

    
    #print('test number:',test)
    return(strain,stress,strain_rate,opt_refl_delay,opt_trans_delay,equilib)

###########################################################################################################################
#
#   Function to identify wave start times for Kbar analysis
#
############################################################################################################################
def find_data_pulse4(time_data,data_data,window_size=400,pulse_to_find=1):
    #this method looks for a sustained change in the level over some window
    # pulse to find is 1=incident, 2=reflected and 3=transmitted
    # incident pulse has a larger requirement for sustainled level (check level) due to strain gage noise problem
    # reflected pulse should be tensile (positive pulse level)
    # incident and transmitted pulses should both have a negative pulse level
    #print('window size:',window_size)
    data=numpy.array(data_data)
    time=numpy.array(time_data)
    #print('length of numpy array')    
    #print(len(data))
    dt=time_data[1]-time_data[0]
    
    delay=int(0.000075/dt) # amount of time ahead of a pulse to include in the pulse
  #  print('delay :',delay)
    #delay=int(0.000050/dt) # amount of time ahead of a pulse to include in the pulse
    checksum_min= 50  # number of time points a pulse must be above a certain level to be considered a pulse
    avg_level=numpy.average(data[0:delay])
    std_level=numpy.std(data[0:delay])
    if pulse_to_find == 1:    
        check_level=0.0002  # this is a check level in absolute volts to capture a pulse
    else:
        check_level=0.0002
    #slope_data=numpy.array(data_data)
   # print('sign:',sign)
   # print('std level =',std_level)
    #calculate the expected value of start time
        
  #  print('avg level:',avg_level)    

    increment=0
    checksum=0
    flag = 0
    
    #provisional start time is the first data point in time-data
    provisional_start_time=time_data[0]    
    start_time=provisional_start_time
  #  print('provisional start time:',provisional_start_time)    
    
    for i in data:  
        if increment>delay:
            level_check = (data[increment]-avg_level)
          #  print('level check:',level_check)
            level_check=abs(level_check)
            slope = data[increment]-data[increment-1]
            slope = abs(slope/dt)
           # print ('slope:',slope)
            if slope < 50: #this determines how many points are above the check level criteria
                if level_check > check_level:
                    if flag == 0:        
                        start = increment
                        flag = 1
                    checksum+=1
  
        increment+=1
 #   print('pulse to find: ',pulse_to_find)


 #   print('pulse to find:',pulse_to_find)
 #   print('start:',start)
 #   print('checksum:',checksum)

    if checksum > checksum_min:     
        try: 
            start_time = start*dt
        except:
            start = delay
            start_time = 0.0
            print('pulse to find: ',pulse_to_find)
            print('***wave start not found***')
        
        try:    
            time_window = time_data[start-delay:start-delay+window_size]
            data_window = data_data[start-delay:start-delay+window_size]
        except:
            time_window = 0
            data_window = 0
            print('***data window not found***')
    
  #  print('start time:',start_time)
 #   print('provisional start time:',provisional_start_time)
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
        #print('increment, for reflected wave:',increment)
        provisional_start_time=time_data[increment] 
        #print('provisional start:',provisional_start_time)

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
        #print('new zero level:',new_avg_level)

        # plot part that I am searching for the reflected wave
       # pylab.clf()    # this is needed to clear data from previous plots
       # pylab.plot(time_subset,data_subset, color = 'r', label='subset sg1')
        #pylab.plot(time_wave1,data_wave1, color = 'b', label='wave1 sg1')
       # pylab.legend()
      #  pylab.show()
        
        
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
                #print('level_check:',level_check,check_level)
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

###################################################################################################
#
#            main function to calculate stress strain from Kolsky bar pulses
#
###################################################################################################    
def analyze_Kolsky_pulses(test,time_data,sg1_strain,sg2_strain,inc_start,refl_delay,trans_delay,window=400,dia=0.004,thk=0.002,gfoil=1,bar_modulus=170):    
    # gfoil=1 means the foil correction is needed. Otherwise it is not
    
    #pre-pulse is the initial zero level time for the pulse windows - in this code it 
    # no longer needs to be at a specific time because the compliance correciton is eliminated
    dt=time_data[1]-time_data[0]
    pre_pulse=int(0.000025/dt)
    #print('dt:',dt)

    #print(refl_delay)
    #print(trans_delay)          
    
    
    #shift the reflected and transmitted wave windows to line up with the incident window based on known inputs
    inc_start=inc_start-pre_pulse*dt
    refl_start=inc_start+refl_delay
    trans_start=inc_start+trans_delay
    
    inc_start_index=int(inc_start/dt)
    refl_start_index=int(refl_start/dt)
    trans_start_index=int(trans_start/dt)
        
    #create pulse arrays
    inc=sg1_strain[inc_start_index:inc_start_index+window]    
    refl=sg1_strain[refl_start_index:refl_start_index+window]    
    trans=sg2_strain[trans_start_index:trans_start_index+window]

    #ensure pulses are all zeroed
    zero_time_points=int(pre_pulse/2)
    inc=zero_timeseries_data(inc,zero_time_points)
    refl=zero_timeseries_data(refl,zero_time_points)
    trans=zero_timeseries_data(trans,zero_time_points)
    
   
    time=numpy.arange(0,window*dt,dt)    # creates time series just for plotting pulses together
    
    
    ####################################################################################################
    # check and adjust polarity of kolsky bar waves using min/max functions. Incident and Trans should be
    #   negative and reflected should be positive in all cases
    
   # inc_zero_level=5.0*numpy.average(inc[0:25])
   # refl_zero_level=5.0*numpy.average(inc[0:25])
   # trans_zero_level=5.0*numpy.average(inc[0:25])

    inc_max = numpy.amax(inc)    
    inc_min = numpy.amin(inc)
    refl_max = numpy.amax(refl)    
    refl_min = numpy.amin(refl)    
    trans_max = numpy.amax(trans)    
    trans_min = numpy.amin(trans)

    
    
    #print('inc_max:',inc_max,inc_zero_level)
    #print('refl_max:',refl_max,refl_zero_level)
    #print('trans_max:',trans_max,trans_zero_level)
    if abs(inc_max) > abs(inc_min):   
        inc=-1.0*inc
        sg1_strain=-1.0*sg1_strain
    if abs(refl_max) < abs(refl_min):
        refl=-1.0*refl
    if abs(trans_max) > abs(trans_min):
        trans=-1.0*trans
        sg2_strain=-1.0*sg2_strain
    
    ###################################################    

    #plot properly signed pulses   
   # pylab.clf()    # this is needed to clear data from previous plots
   # pylab.plot(time, inc, color = 'b', label='INC')
   # pylab.plot(time, trans, color = 'r', label='TRANS')
   # pylab.plot(time, refl, color = 'y', label='REFL')
   # pylab.legend()    
   # pylab.show()

    
    #create data for equlibrium plot
    inc_plus_refl = inc+refl
    equilib = inc_plus_refl-trans
    criterion_initial=numpy.std(equilib)

    
    
    #scan for a better equilibrium value and adjust waves to match this value
    scan_range = 100 #number of time points to scan   
    j_opt=0
    i_opt=0
    criterion_old=100
    for i in range(0,scan_range): 
        n=-scan_range/2+i        
        refl=sg1_strain[refl_start_index-n:refl_start_index+window-n]
        for j in range(0,scan_range):
            m=-scan_range/2+j  
            trans=sg2_strain[trans_start_index-m:trans_start_index+window-m]
            equilib=inc+refl-trans
            criterion = numpy.std(equilib) #this selects for flatness of the wave difference only
            if criterion < criterion_old:
                j_opt=j
                i_opt=i
                criterion_old=criterion
            #print('criterion: ', criterion)        
    # these optimum values are in units of time and are used to pick better time windows for
    #    the reflected and transmitted waves. the incident wave window is fixed
    n_opt=-scan_range/2+i_opt    
    m_opt=-scan_range/2+j_opt
    
    #this disables the scan for better equilibrium
    #n_opt=0
    #m_opt=0


    #refl_delay=inc_start_index-refl_start_index-n_opt
    #trans_delay=inc_start_index-trans_start_index-m_opt

    #refl_delay_index=inc_start_index-refl_start_index-n_opt
    #trans_delay_index=inc_start_index-trans_start_index-m_opt
    
    # assign optimized wave windows for analyzing stress-strain response    
    refl=sg1_strain[refl_start_index-n_opt:refl_start_index+window-n_opt]   
    trans=sg2_strain[trans_start_index-m_opt:trans_start_index+window-m_opt]
    
    opt_refl_delay=time_data[refl_start_index-n_opt]-time_data[inc_start_index]   
    opt_trans_delay=time_data[trans_start_index-m_opt]-time_data[inc_start_index]   

 #   print('n_opt:',n_opt)
  #  print('m_opt:',m_opt)    
  #  print('refl delay:',opt_refl_delay)
  #  print('trans_delay:',opt_trans_delay)

    
    #ensure pulses are all zeroed
    inc=zero_timeseries_data(inc,25)
    refl=zero_timeseries_data(refl,25)
    trans=zero_timeseries_data(trans,25)
 
    inc_plus_refl = inc+refl
    equilib = inc_plus_refl-trans
    criterion_final=numpy.std(equilib)
    #print('final equilib stdev:',criterion_final)
    
     
    #
    #  Plot this first time through a dataset to see how it looks
    #
 
    #pylab.clf()    # this is needed to clear data from previous plots
   # pylab.title("Wave Equlibrium")
   # pylab.xlabel("Time [s]")
   # pylab.ylabel("Strain")
       # pylab.legend()
   # pylab.show()
    
    
    #calculate sample strain, strain rate and stress histories using aligned waves    
    # First define analysis constants - all lengths in m
    bar_diameter=0.015
    bar_area=3.14159/4*bar_diameter*bar_diameter
    sample_area=3.14159/4*dia*dia
    density=8048.0
    #wave_speed=4584 # in m/s
    wave_speed=bar_modulus*1.0e9/density
    wave_speed=wave_speed.__pow__(0.5)
    #print('wave_speed:',wave_speed)
    #modulus=170 # in GPa
    #modulus=190 # in GPa
  
    #define strain rate and strain increment and stress vs time arrays using 1-wave analysis
    window=int(window)
    eng_strn=range(0,window,1)    
        
    eng_strn_rate=2*wave_speed/thk*refl
    eng_stress=-bar_area*bar_modulus/sample_area*trans*1000 #[MPa]
    eng_strn_increment=eng_strn_rate*dt
    
    # sum up engineering strain    
    sum=0.0
    #sum up total engineering strain
    for i in range(0,window,1):             
        sum+=eng_strn_increment[i] 
        eng_strn[i]=sum              
        
    #correct engineering strain for foil and compliance if needed; note foil and compliance are in terms of strain
    
    #call Zhao-Gary Punch Correction for All tests
    #
    ZG_Punch=numpy.arange(0,window,dtype=numpy.float) #define in case compliance correction is not invoked
    #ZG_Punch=ZG_Punch*0.0
    ZG_Punch=punch_correction(window,dt,eng_stress,dia,bar_modulus) #this calls the punch correction - need to figure out what to send to it
            
    #Calculate grafoil correction if doing pulse-heated tests
            
    foil=numpy.arange(0,window,dtype=numpy.float) #define in case foil correction is not invoked
    foil=foil*0.0

    #print(gfoil)
    # calculate foil correction if gfoil =1 meaning foil is used.
    if gfoil == 1:     
        foil=grafoil_correction(window,dt,eng_strn,eng_stress,thk,ZG_Punch)  
        #print('I am making a foil Correction!')
    
    eng_strn_corrected=numpy.arange(0,window,dtype=numpy.float)  #this just initializes an array of the right size
      
    for i in range(0,window,1):
        eng_strn_corrected[i]=eng_strn[i]-(ZG_Punch[i]+foil[i])/thk
        #eng_strn_corrected[i]=eng_strn[i]
        
    true_stress=(1-eng_strn_corrected)*eng_stress

    true_strn=numpy.arange(0,window,dtype=numpy.float)  #this just initializes a sequence of the right size   
    true_strn=-numpy.log(1-eng_strn_corrected)
   # print(eng_strn_corrected)
   
   # if foil-corrected true strains are negative or infinte, set them to zero    
    for i in range(0,window,1):
        if true_strn[i]>0:
            true_strn[i]=true_strn[i]
        else:
            true_strn[i]=0.0
        
    #estimate true strain rate by computing instantaneous (current) sample length
    #current_sample_l=range(0,window,1)  #this just initializes a sequence of the right size
 
    #This way to declare seems to solve the problems that occured using the earlier declaration method 
    current_sample_l=numpy.arange(0,window,dtype=numpy.float)
    for i in range(0,window,1):             
        current_sample_l[i]=numpy.exp(-true_strn[i])*thk
    true_strn_rate=2*wave_speed/current_sample_l*refl
    
    #pick a stress value to see if it is positive or negative
    stress_time=int(0.0001/dt)
    stress_interrogation = true_stress[stress_time]
    print('True stress at 100 µs :',stress_interrogation)
    
   
    #pylab.clf()    
    #pylab.plot(time, eng_strn, color = 'y', label='total')
    #pylab.plot(time, eng_strn_corrected, color = 'b', label='sample')
    #pylab.plot(time, foil, color = 'b', label='foil')
    #pylab.plot(time, ZG_Punch, color = 'b', label='ZG Punch')
    
    #pylab.plot(eng_strn, eng_stress, color = 'b', label='stress')
    


    #pylab.plot(eng_strn_corrected, eng_stress, color = 'y', label='eng stress')      
   # pylab.plot(true_strn, true_stress, color = 'y', label='stress')
   # pylab.legend()
   # pylab.show()



    #pylab.plot(true_strn, true_strn_rate, color = 'b', label='strain rate')

   # pylab.plot(time, true_stress, color = 'y', label='true')
    #pylab.title("Wave Equlibrium")
    #pylab.xlabel("Time [s]")
    #pylab.ylabel("Strain")
    #pylab.legend()

########## For Plotting Individual Curves

    #compute a nominal elastic loading curve for steel
    specimen_modulus = 200.0*1000.0 #in MPa
    max_elastic_stress=numpy.amax(true_stress)
    max_elastic_strain=max_elastic_stress/specimen_modulus
    max_strain=max_elastic_strain=0.005
    elastic_strain=[0,max_strain]
    elastic_stress=[0,max_strain*specimen_modulus]
  #  plt.clf()
    Title_1=str(test)+' Equilibrium'
    Title_2=str(test)+' True Stress-Strain'
    Title_3=str(test)+' True Strain-Rate'
    Title_4=str(test)+' Eng Stress-Strain' 

    plt.figure()
    
    plt.subplot(221)
    plt.title(Title_1)
    plt.plot(time,inc_plus_refl)
    plt.plot(time,trans)
    plt.plot(time,equilib)
    
    plt.subplot(222)
    plt.title(Title_2)
    plt.plot(true_strn,true_stress)
    plt.plot(elastic_strain,elastic_stress)

    plt.subplot(223)
    plt.title(Title_3)
    plt.plot(true_strn,true_strn_rate)
    
    plt.subplot(224)
    plt.title(Title_4)
    plt.plot(eng_strn_corrected,eng_stress)
    plt.plot(elastic_strain,elastic_stress)

    plt.show()
    
    
    
   
    return(true_strn,true_stress,true_strn_rate,opt_refl_delay,opt_trans_delay,criterion_final)
    #return(eng_strn_corrected,eng_stress,eng_strn_rate,opt_refl_delay,opt_trans_delay,criterion_final)
      
    #return(eng_strn,eng_stress,true_strn_rate)
    #return(inc,refl,trans)


#------------------END OF ANALYZE KOLSKY PULSES-------------------------------

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
    start = 0.0
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

    
      
    

def punch_correction(window,dt,eng_stress,dia,bar_modulus):
    # From Zhao Gary paper - used for all KBar compression tests to account for elastic punching of the specimen into 
    #   the end of the bar. Takes as input the elastic properties of the bar and of the specimen, and the size ratio, 
    #   and the specimen stress. It should be the specimen true stress, but for now it will be the engineering stress.
    #   ALSO - for heated tests need to pass the high temperature modulus of both specimen and bar to accurately account
    #   for punching.
    
    #E_bar=190000000000 #Young's modulus of bar tip in Pa (need to pass high T value)
    #E_bar=190000000000 #Young's modulus of bar tip in Pa (need to pass high T value)
    E_bar=bar_modulus*1.0e9        
    Poisson = 0.29 #Poisson's ratio for the bar tip
    D_bar=15.0/1000 #bar diameter in m
    D_sample = dia  #sample diameter passed from code
    
    # These are basically geometric constants in the Zhao Gary formula
    Hp = 1.416-(1.301-1.416)/(0.3-0.25)*(D_sample/D_bar-0.25) #this is an interpolation function
    Kp = Hp*16/(3*3.14159*3.14159)*(1-Poisson*Poisson)/(D_sample*E_bar)   
    
    punch_displ=numpy.arange(0,window,dtype=numpy.float)  #this just initializes a sequence of the right size    
    
    for i in range(0,window,1):             
        punch_displ[i]=2*Kp*eng_stress[i]*3.14159/4*(D_sample*D_sample)*1e6   #last factor is to scale eng stress to Pa
        #print(i,eng_stress[i])

    return(punch_displ)
#----------------------------------------------------------------------------    
    
#----------------------------------------------------------------------------       
def grafoil_correction(window,dt,eng_strn,eng_stress,sample_thickness,ZG_Punch):

    # foil iteration 1 uses engineering stress 
    eng_strn=numpy.array(eng_strn)       
    foil=numpy.arange(0,window,dtype=numpy.float)  #this just initializes a sequence of the right size
    for i in range(0,window,1):             
        foil[i]=(2.10098e-5)*eng_stress[i]**0.333        
        if foil[i] > 0.00026:        
            foil[i] = 0.00026   
      
    for i in range(0,window,1):
        eng_strn[i]=eng_strn[i]-(foil[i]+ZG_Punch[i])/sample_thickness
     
    
    # foil iteration 1: true-stress vs time        
    true_stress=numpy.arange(0,window,dtype=numpy.float)  #this just initializes a sequence of the right size
    true_stress=(1-eng_strn)*eng_stress
   

    # foil iteration 2 uses true stress        
    foil=range(0,window,1)  #this just initializes a sequence of the right size
    for i in range(0,window,1):             
        foil[i]=(2.10098e-5)*true_stress[i]**0.333        
        if foil[i] > 0.00026:        
            foil[i] = 0.00026   
      
    for i in range(0,window,1):
        eng_strn[i]=eng_strn[i]-(foil[i]-ZG_Punch[i]/1000)/sample_thickness

    true_stress=(1-eng_strn)*eng_stress

    # foil iteration 3 - final correction    
    foil=range(0,window,1)  #this just initializes a sequence of the right size
    for i in range(0,window,1):             
        foil[i]=(2.10098e-5)*true_stress[i]**0.333        
        if foil[i] > 0.00026:        
            foil[i] = 0.00026  
                        
    return(foil)
 #----------------------------------------------------------------------------       
   
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
    strain_for_plot = [0.05,0.075,0.1,0.15] # value of total strain for plotting thermal softening
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
            ss_filename=expt + '_analyzed.csv'
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
