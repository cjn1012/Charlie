# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:05:48 2018

@author: cjn7
"""

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')
    
    
temp   = [-18,21,93,204,316,427,538,649,760,871,982,1093] # Degree C
spheat = [402,410,427,456,481,511,536,565,590,620,645,670]# J/Kg C

temp2 = [70,200,400,600,800,1000] # Degree C
modofe= [30.1,29.6,28.7,27.6,26.9,25.9] # ksi

import numpy as np
import matplotlib.pyplot as plt


m,b = np.polyfit(temp,spheat,1)
m2,b2 = np.polyfit(temp2,modofe,1)


plt.subplot(2, 1, 1)
plt.scatter(temp,spheat)
plt.xlabel('Temperature (Celcius)' )
plt.ylabel('Specific Heat (J/KgC)')
plt.title('Properties of Inconel 625 with Varying Temperature')
plt.xticks(np.arange(-200,1400,step=200))
abline(m,b)
plt.text(-200,600,'y=' + str(round(m,5)) + 'x' + str(round(b,4)))

plt.subplot(2, 1, 2)
plt.scatter(temp2,modofe)
plt.xlabel('Temperature (Celcius)' )
plt.ylabel('Elastic Modulus (ksi)')
plt.xticks(np.arange(-200,1400,step=200))
abline(m2,b2)
plt.text(-200,26,'y=' + str(round(m2,5)) + 'x' + str(round(b2,4)))