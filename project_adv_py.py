import numpy as np
import struct
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm



#Define the functions to pick the smp to start and stop counting the rise time 

def return_beg(y):
    ''' This function receives a range of the the signal around the pre_trg_smp
    and returns the first sample before rising. 

    This is done by calculating the
    derivative on each point of this segment of the signal, then, an array of bools
    are created to check if the derivative is growing. In general, the derivative is 
    greater than 2 when the signal starts rising. If there are 6 consecutive "trues"
    this is enough for the function to understand that the signal is rising and return
    the first sample (smp)'''

    dxdy = np.gradient(y)
    bool_arr = np.array([dxdy[i] > 2 for i in range(len(dxdy)-1)], dtype = bool)
    consecutive_true = np.convolve(bool_arr, np.ones(6), mode='valid') == 6
    start_smp = (np.argmax(consecutive_true) + 6) if any(consecutive_true) else 0
    return (490 + start_smp)-7

def return_end(x,y,smp2extrpl_x,smp2extrpl_y):
    '''This function receives a piece of the signal after the pre_trg_smp
    (x, y) and a piece of signal when the signal is decaying and returns 
    the final smp after rising.

    To do this the function performs a linear fit on this piece of signal
    when is decaying. Then, an extrapolation to the piece of signal after 
    the pre_trg_smp is performed to determine the last point after the signal 
    rise avoiding the overshot. After this the fuction finds the 
    closest value to the extrapolation and returns the sample number of this 
    value'''
    slope, intercept = np.polyfit(x,y,1)
    extrapolation = np.array(intercept + slope*smp2extrpl_x)
    smp_inter = np.argmin(np.abs(smp2extrpl_y - extrapolation))
    return smp_inter+500

# Define the Gaussian function for the fit
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * stddev ** 2))


# Find all files with the pattern "Mass_100_Energy_*.bin" in the directory
# file_pattern = '/media/integra/acd1b1cd-2cf9-4837-9619-4230c7a64b1d/Documents/2021_ILL/Analysis/Det_99_Binary/230321/Mass_100_Energy_*.bin'
file_pattern = '230313//Mass_100_Energy_*.bin'

# Initialize an empty list to hold the data
data = []

# Initialize arrays to store the means of the found rise times dists, the energy and their uncertainties
mean_rise_time_per_e = np.empty(shape=(0,), dtype=float)
mean_rise_time_per_e_unc = np.empty(shape=(0,), dtype=float)

mean_smp_start_per_e = np.empty(shape=(0,), dtype=float)
mean_smp_start_per_e_unc = np.empty(shape=(0,), dtype=float)

mean_smp_stop_per_e = np.empty(shape=(0,), dtype=float)
mean_smp_stop_per_e_unc = np.empty(shape=(0,), dtype=float)

energy = np.empty(shape=(0,), dtype=float)
energy_unc = np.empty(shape=(0,), dtype=float)

#A counter so that it can be traked in which file the program is
counter = 0 

# Loop over all files that match the pattern and read the data
for file_name in glob.glob(file_pattern):
    with open(file_name, 'rb') as f:
        buffer = f.read()
        f.close()
        print("Length of buffer is %d" % len(buffer))

    # Find the size
    size = str(int(len(buffer) / 8)) + 'd'
    print("size: ", size)

    # Unpack the data
    data = struct.unpack(size, buffer)

    # Find the number of signals in this file
    number_of_signals = int((len(buffer) - 9 * 8) / (8 * 2000))
    print(number_of_signals)

    #fill the energy and energy unc vectors with the information of the signals in the heading of the bin file
    energy = np.append(energy, data[4])
    energy_unc = np.append(energy_unc, data[5])

    #Remove the hading of the bin file to have only the signals
    data = np.array(data[9:])
    num_elements = number_of_signals * 2000
    #Reshape everything to a 2D array with the total number of signals per file times the smp_per_sgnal
    data = data[:num_elements].reshape(number_of_signals, 2000)

    #changing the name of "data" for something more meaningful to plot and creating an array of 2000 smp to plot as well
    y = data
    x = np.arange(0, 2000)

    ############# These steps are to check if the bin files are correct #########

    # # Compute the maximum value for each signal in y (to check if the binary files are correct, namely, there is one peak per file)
    # max_values = np.max(y[0:number_of_signals], axis=1)

    # # Compute the maximum value among all elements in y
    # global_max_value = np.max(y[0:number_of_signals,:])

    # # Create a histogram of the max values
    # num_bins = 50
    # plt.hist(max_values, bins=num_bins)
    # plt.xlabel('Max value')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of max values in signals')
    # plt.show()

    # print("Max value among all elements in y:", global_max_value)

    # Plot the first signal
    # print("Max value of first signal:", np.max(y))
    # if counter == 3:
    #     plt.plot(x, y[3179], color="red")
    #     plt.show()

    #############################################################################

    #declaring the arrays to start, stop and rise_time
    smp_start = np.empty(number_of_signals, dtype=int)
    smp_stop = np.empty(number_of_signals, dtype=int)
    rise_time = np.empty(number_of_signals, dtype=int)

    #for loop to calculate the rise time signal per signal
    for i in range(number_of_signals):
        smp_start[i] = return_beg(y[i,490:520])
        smp_stop[i] = return_end(x[1800:2000], y[i,1800:2000], x[500:600], y[i,500:600])
        rise_time[i] = smp_stop[i] - smp_start[i]

    # #Calculations for the histogram construction
    # # total number of bins
    # num_bins = int(np.max(rise_time) - np.min(rise_time))
    
    # #bin edges from the min bin to the max bin in rise time
    # bin_edges = np.arange(np.min(rise_time)-0.5, np.max(rise_time)+0.5)
    # bin_center = (0.5*(bin_edges[1:]+bin_edges[:-1]))
    # #hist containf the frecuencies per bin
    # hist, _ = np.histogram(rise_time, bins=bin_edges)

    # #guess for the fit
    # amplitude_guess = np.max(hist)
    # mean_guess = np.mean(rise_time)
    # sigma_guess = np.std(rise_time)

    # #To plot a function of the guess and check if the guess is good
    # x_gauss = np.arange(np.min(rise_time), np.max(rise_time), 0.1)
    # y_gauss = gaussian(x_gauss, amplitude_guess, mean_guess, sigma_guess)

    # # #Performs the gaussian fit
    # # coeff, cov = curve_fit(gaussian, bin_center, hist, p0=(amplitude_guess, mean_guess, sigma_guess))

    # # #The fit parameters are collected in this np arrays
    # # mean_rise_time_per_e = np.append(mean_rise_time_per_e, coeff[1])
    # # mean_rise_time_per_e_unc = np.append(mean_rise_time_per_e_unc, np.sqrt(cov[1,1 ]))

    # plt.hist(rise_time, bins = bin_edges)
    # plt.plot(bin_center, hist, "r.")
    # plt.plot(x_gauss, y_gauss, 'r--', label='Gaussian function')
    # plt.show()

    #################################################################################################

    num_bins = int(np.max(smp_start) - np.min(smp_start))
    
    #bin edges from the min bin to the max bin in rise time
    bin_edges = np.arange(np.min(smp_start)-0.5, np.max(smp_start)+0.5)
    bin_center = (0.5*(bin_edges[1:]+bin_edges[:-1]))
    #hist containf the frecuencies per bin
    hist, _ = np.histogram(smp_start, bins=bin_edges)

    #guess for the fit
    amplitude_guess = np.max(hist)
    mean_guess = np.mean(smp_start)
    sigma_guess = np.std(smp_start)

    # #To plot a function of the guess and check if the guess is good
    x_gauss = np.arange(np.min(smp_start), np.max(smp_start), 0.1)
    y_gauss = gaussian(x_gauss, amplitude_guess, mean_guess, sigma_guess)

    #Performs the gaussian fit
    coeff, cov = curve_fit(gaussian, bin_center, hist, p0=(amplitude_guess, mean_guess, sigma_guess))

    #The fit parameters are collected in this np arrays
    mean_smp_start_per_e = np.append(mean_smp_start_per_e, coeff[1])
    mean_smp_start_per_e_unc = np.append(mean_smp_start_per_e_unc, np.sqrt(cov[1,1 ]))

    # plt.hist(smp_start, bins = bin_edges)
    # plt.plot(bin_center, hist, "r.")
    # plt.plot(x_gauss, y_gauss, 'r--', label='Gaussian function')
    # plt.show()

    # #################################################################################################

    num_bins = int(np.max(smp_stop) - np.min(smp_stop))
    
    #bin edges from the min bin to the max bin in rise time
    bin_edges = np.arange(np.min(smp_stop)-0.5, np.max(smp_stop)+0.5)
    bin_center = (0.5*(bin_edges[1:]+bin_edges[:-1]))
    #hist containf the frecuencies per bin
    hist, _ = np.histogram(smp_stop, bins=bin_edges)

    #guess for the fit
    amplitude_guess = np.max(hist)
    mean_guess = np.mean(smp_stop)
    sigma_guess = np.std(smp_stop)

    # #To plot a function of the guess and check if the guess is good
    x_gauss = np.arange(np.min(smp_stop), np.max(smp_stop), 0.1)
    y_gauss = gaussian(x_gauss, amplitude_guess, mean_guess, sigma_guess)

    # Performs the gaussian fit
    coeff, cov = curve_fit(gaussian, bin_center, hist, p0=(amplitude_guess, mean_guess, sigma_guess))

    #The fit parameters are collected in this np arrays
    mean_smp_stop_per_e = np.append(mean_smp_stop_per_e, coeff[1])
    mean_smp_stop_per_e_unc = np.append(mean_smp_stop_per_e_unc, np.sqrt(cov[1,1 ]))

    # plt.hist(smp_stop, bins = bin_edges)
    # plt.plot(bin_center, hist, "r.")
    # plt.plot(x_gauss, y_gauss, 'r--', label='Gaussian function')
    # plt.show()

    #increase the file number
    counter = counter + 1

# # Create the plot with error bars
# plt.errorbar(energy, mean_rise_time_per_e, xerr=energy_unc, yerr=mean_rise_time_per_e_unc, fmt='o')

# # Add axis labels and title
# plt.xlabel('energy')
# plt.ylabel('mean_rise_time_per_e')
# plt.title('Rise time of the signals as a function of the energy for mass 100')

# # Display the plot
# plt.show()

####################################################################################################

# Create the plot with error bars
plt.errorbar(energy, mean_smp_start_per_e, xerr=energy_unc, yerr=mean_smp_start_per_e_unc, fmt='o')

# Add axis labels and title
plt.xlabel('energy')
plt.ylabel('mean_smp_start_per_e')
plt.title('Rise time of the signals as a function of the energy for mass 100')

# Display the plot
plt.show()

####################################################################################################

# Create the plot with error bars
plt.errorbar(energy, mean_smp_stop_per_e, xerr=energy_unc, yerr=mean_smp_stop_per_e_unc, fmt='o')

# Add axis labels and title
plt.xlabel('energy')
plt.ylabel('mean_smp_stop_per_e')
plt.title('Rise time of the signals as a function of the energy for mass 100')

# Display the plot
plt.show()