"""

Created by Robert Wolle, 9/19
Program that reads in data and fits it to Lorentzian peaks

"""


import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.

def Lorentzian(x, center, width, height, background=0):
    # Half height is at x = center +/- 1/2 * width  (full width at half max) 

    y = (height-background)*width**2/(4*(x-center)**2 + width**2) + background/2
    return y

def Gaussian(x, center, width, height, background=0):

    width = width*np.sqrt(np.log(4))
    y = background + (height-background)*np.exp(-2*(x-center)**2/width**2)
    return y


def make_DataArray(delimiter, data_columns, data_start, file_path):
    with open(file_path, "r") as file:
        data = file.read()

    data_values = []
    last_break = data_start
    for char in range(data_start,len(data)):
        if data[char] == delimiter or data[char] == '\n':
            data_values.append(float(data[last_break:char]))
            last_break = char+1

    # We can use the number of data points divided by number of columns to get the number of rows of data
    data_rows = int(len(data_values)/data_columns)

    data_organized = []
    for row in range(data_rows):
        data_organized.append(data_values[row*data_columns:(row+1)*data_columns])
    return data_organized

def nm_to_eV(wavelength):
    return 1240/wavelength


# User inputs:
delimiter = '\t'
data_columns = 3
data_start = 521     #character position at which the data begins
file_path = "/workspaces/Projects/Physics/PeakFitting/WSe2_MoSe2_NbSe2_100uW_PL_NoFlake2Test6.vdat"


DataArray = np.array(make_DataArray(delimiter, data_columns, data_start, file_path))
x_axis = nm_to_eV(DataArray[:,0])
y_axis = DataArray[:,2]

x_region = [1.4,1.75]

plt.plot(x_axis,y_axis,zorder=2)
plt.gca().invert_xaxis()
plt.xlim(x_region)
plt.title('PL Intensity off NbSe2, 10s exposure')
plt.xlabel('Energy (eV)')
plt.ylabel('Counts')


# Input peak parameters
n_peaks = 1
Background = 600

# Linear Background
Background = np.zeros(len(x_axis))
constBackground = 1000
slopeBackground = 0.2


# Formatted [[center_1,width_1,height_1], [center_2,width_2,height_2],...,[center_n,width_n,height_n]]
peaks = [[1.573, 0.045,16000], [1.663,0.038,0]]

def peakFits(peaks, n_peaks, background, x):
    peakFit = 0
    for n in range(n_peaks):
        peakFit = peakFit + Lorentzian(x, peaks[n][0], peaks[n][1], peaks[n][2], background)
    return peakFit
    

loss = 0
fit = np.zeros(len(x_axis)) + Background

for i in range(len(x_axis)):
    x_i = x_axis[i]
    if x_i >= x_region[0] and x_i <= x_region[1]:
        Background[i] = constBackground + slopeBackground*i
        fit[i] = peakFits(peaks, n_peaks, Background[i], x_i)
        loss = loss + (fit[i] - y_axis[i])**2/y_axis**2

plt.plot(x_axis,fit,zorder=1)
plt.text(1.525, 100, 'Peak at 1.573 eV', fontsize=11)
plt.savefig('/workspaces/Projects/Physics/PeakFitting/output.png')





