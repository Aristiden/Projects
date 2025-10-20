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
file_path = "/workspaces/Projects/Physics/PeakFitting/WSe2_MoSe2_NbSe2_100uW_PL_Flake2Test2.vdat"


DataArray = np.array(make_DataArray(delimiter, data_columns, data_start, file_path))
x_axis = nm_to_eV(DataArray[:,0])
y_axis = DataArray[:,2]

x_region = [1.4,1.75]

plt.plot(x_axis,y_axis,zorder=2)
plt.gca().invert_xaxis()
plt.xlim(x_region)
plt.title('PL Intensity on NbSe2, 10s exposure')
plt.xlabel('Energy (eV)')
plt.ylabel('Counts')


# Input peak parameters
n_peaks = 2
Background = 600

# Linear Background
Background = np.zeros(len(x_axis))
constBackground = 500
slopeBackground = 0.1


# Formatted [[center_1,width_1,height_1], [center_2,width_2,height_2],...,[center_n,width_n,height_n]]
peaks = [[1.574, 0.043,16500], [1.663,0.03,1400]]

def peakFits(peaks, n_peaks, background, x):
    peakFit = 0
    for n in range(n_peaks):
        peakFit = peakFit + Lorentzian(x, peaks[n][0], peaks[n][1], peaks[n][2], background)
    return peakFit
    

loss_previous = 0

peaks_i = np.array(peaks)
peaks_previous = np.array(peaks)


loss_diff = 0
features_diff = []
fit = np.zeros(len(x_axis)) + Background


def calculateFits(x_axis, x_region, peaks, n_peaks, Background, loss):
    for i in range(len(x_axis)):
                x_i = x_axis[i]
                if x_i >= x_region[0] and x_i <= x_region[1]:
                    Background[i] = constBackground + slopeBackground*i
                    fit[i] = peakFits(peaks, n_peaks, Background[i], x_i)
                    loss = loss + (fit[i] - y_axis[i])**2/y_axis[i]**2
    return fit, loss

fit,loss_previous = calculateFits(x_axis, x_region, peaks, n_peaks, Background, loss_previous)


iterations = 10
change_rate = 0.01      # iterations change value of the feature by change_rate*100% 





# First change:
loss_i = 0
for n in range(n_peaks):
    for k in range(3):
        peaks_i[n][k] = peaks_previous[n][k] + change_rate*peaks_previous[n][k]
        fit,loss_i = calculateFits(x_axis, x_region, peaks_i, n_peaks, Background, loss_i)

loss_diff = loss_i - loss_previous
features_diff = peaks_i - peaks_previous


peaks_previous = peaks_i
loss_previous = loss_i

# Successive iterations
for i in range(iterations):
    loss_i = 0
    for n in range(n_peaks):
        for k in range(3):
            print(np.linalg.norm(features_diff))
            peaks_i[n][k] = peaks_previous[n][k] + loss_diff/np.linalg.norm(features_diff)*peaks_previous[n][k]
    fit,loss_i = calculateFits(x_axis, x_region, peaks_i, n_peaks, Background, loss_i)
    loss_diff = loss_i - loss_previous
    features_diff = peaks_i - peaks_previous
    peaks_previous = peaks_i
    loss_previous = loss_i


print(loss_i)


plt.plot(x_axis,fit,zorder=1)
plt.text(1.525, 100, 'Peak at 1.573 eV', fontsize=11)
plt.savefig('/workspaces/Projects/Physics/PeakFitting/output.png')





