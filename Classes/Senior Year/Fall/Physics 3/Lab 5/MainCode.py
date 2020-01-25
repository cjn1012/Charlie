# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:03:55 2019

@author: User
"""
############
#Problem 1
############













############
#Problem 2
############

import matplotlib.pyplot as plt
import math
import numpy as np
text_file = open("millikan.txt", "r")

lines = text_file.readlines()

x_array = []
y_array = []
for i in lines:
    x = i.split()
    print(x)
    x_array.append(float(x[0]))
    y_array.append(float(x[1]))
    x_array = list(map(float,x_array))
    y_array = list(map(float,y_array))
plt.plot(x_array,y_array,'o')
plt.plot(x_array,y_array)


def lineFit(x,y):
    print(x)
    Ex = (1/len(x))*sum(x)
    Ey = (1/len(y))*sum(y)
    
    Exx = (1/len(x))*sum(x*x)
    Exy = (1/len(y))*sum(x*y)
    
    m = (Exy-Ex*Ey)/(Exx-(Ex*Ex))
    return(m)
    
def lineFit2(x,y):
    Ex = (1/len(x))*sum(x)
    Ey = (1/len(y))*sum(y)
    
    Exx = (1/len(x))*sum(x*x)
    Exy = (1/len(y))*sum(x*y)
    
    b = (Exx*Ey-Ex*Exy)/(Exx-(Ex*Ex))
    return(b)    


m = lineFit(np.array(x_array),np.array(y_array))
b = lineFit2(np.array(x_array),np.array(y_array))
print(m)
print(b)  
fit_y = []
for i in x_array:
    fit_y.append(m*i + b)

plt.plot(x_array,fit_y)
plt.show()

# Part d in problem 2
e = 1.602176*pow(10,-19)
counter = [0,1,2,3,4,5]
h = []

v = x_array[1]
Volt = y_array[1]
work = 0
h = (e*Volt+work)/v
    
print(h)

def relativeError(xmeasured,xaccepted):
    Relative_Error = abs(float((xmeasured-xaccepted)/xaccepted))
    return(Relative_Error)
    
hreal = 6.62607015*pow(10,-34)
error = relativeError(h,hreal)
print(error)


############
#Problem 3
############











############
#Problem 4
############

text_file = open("ResistorandCurrent.txt", "r")

lines = text_file.readlines()

Resistors = lines[0:5]
Currents = lines[6:11]






















    