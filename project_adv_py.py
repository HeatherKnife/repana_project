import numpy as np
import struct
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm

# Find all files with the pattern "Mass_100_Energy_*.bin" in the directory
file_pattern = '/media/integra/acd1b1cd-2cf9-4837-9619-4230c7a64b1d/Documents/2021_ILL/Analysis/Det_97_Binary/230313/Mass_100_Energy_*.bin'

# Initialize an empty list to hold the data
data = []
mean_rise_time_per_e = np.empty(shape=(0,), dtype=float)
mean_rise_time_per_e_unc = np.empty(shape=(0,), dtype=float)
energy = np.empty(shape=(0,), dtype=float)
energy_unc = np.empty(shape=(0,), dtype=float)

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

    print("data[4] ", data[4])
    print("data[5] ", data[5])

    energy = np.append(energy, data[4])
    energy_unc = np.append(energy_unc, data[5])

    data = np.array(data[9:])
    num_elements = number_of_signals * 2000
    data = data[:num_elements].reshape(number_of_signals, 2000)

    # Stack all the data together into a single 3D array
    y = data
    x = np.arange(0, 2000)

    # Compute the maximum value for each signal in y
    max_values = np.max(y[0:number_of_signals], axis=1)

    # Compute the maximum value among all elements in y
    global_max_value = np.max(y[0:number_of_signals,:])

    # Create a histogram of the max values
    # num_bins = 50
    # plt.hist(max_values, bins=num_bins)
    # plt.xlabel('Max value')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of max values in signals')
    # plt.show()

    print("Max value among all elements in y:", global_max_value)

    # Plot the first signal

    def return_beg(y):
        dxdy = np.gradient(y)
        bool_arr = np.array([dxdy[i] > 2 for i in range(len(dxdy)-1)], dtype = bool)
        consecutive_true = np.convolve(bool_arr, np.ones(6), mode='valid') == 6
        start_smp = (np.argmax(consecutive_true) + 6) if any(consecutive_true) else 0
        return (490 + start_smp)-7

    def return_end(x,y,smp2extrpl_x,smp2extrpl_y):
        slope, intercept = np.polyfit(x,y,1)
        extrapolation = np.array(intercept + slope*smp2extrpl_x)
        smp_inter = np.argmin(np.abs(smp2extrpl_y - extrapolation))
        return smp_inter+500

    smp_start = np.empty(number_of_signals, dtype=int)
    smp_stop = np.empty(number_of_signals, dtype=int)
    rise_time = np.empty(number_of_signals, dtype=int)

    # Plot the first signal
    print("Max value of first signal:", np.max(y))
    # if counter == 3:
    #     plt.plot(x, y[3179], color="red")
    #     plt.show()

    for i in range(number_of_signals):
        smp_start[i] = return_beg(y[i,490:520])
        smp_stop[i] = return_end(x[1800:2000], y[i,1800:2000], x[500:600], y[i,500:600])
        rise_time[i] = smp_stop[i] - smp_start[i]

    print(rise_time[0:10])

    num_bins = int(np.max(rise_time) - np.min(rise_time))
    bin_edges = np.arange(np.min(rise_time)-0.5, np.max(rise_time)+0.5)
    bin_center = (0.5*(bin_edges[1:]+bin_edges[:-1]))
    hist, _ = np.histogram(rise_time, bins=bin_edges)

    amplitude_guess = np.max(hist)
    mean_guess = np.mean(rise_time)
    sigma_guess = np.std(rise_time)

    # Define the Gaussian function
    def gaussian(x, amplitude, mean, stddev):
        return amplitude * np.exp(-(x - mean) ** 2 / (2 * stddev ** 2))

    x_gauss = np.arange(np.min(rise_time), np.max(rise_time), 0.1)
    y_gauss = gaussian(x_gauss, amplitude_guess, mean_guess, sigma_guess)

    coeff, cov = curve_fit(gaussian, bin_center, hist, p0=(amplitude_guess, mean_guess, sigma_guess))

    mean_rise_time_per_e = np.append(mean_rise_time_per_e, coeff[1])
    mean_rise_time_per_e_unc = np.append(mean_rise_time_per_e_unc, np.sqrt(cov[1,1 ]))

    print("coeff: ", coeff)
    print("cov: ", cov)

    # plt.hist(rise_time, bins = bin_edges)
    # plt.plot(bin_center, hist, "r.")
    # plt.plot(x_gauss, y_gauss, 'r--', label='Gaussian function')
    # plt.show()

    counter = counter + 1
    print(counter)

print("energy ", energy)
print("energy_unc ", energy_unc)
print("mean_rise_time_per_e ", mean_rise_time_per_e)
print("mean_rise_time_per_e_unc ", mean_rise_time_per_e_unc)

# Create the plot with error bars
plt.errorbar(energy, mean_rise_time_per_e, xerr=energy_unc, yerr=mean_rise_time_per_e_unc, fmt='o')

# Add axis labels and title
plt.xlabel('energy')
plt.ylabel('mean_rise_time_per_e')
plt.title('Example plot with error bars')

# Display the plot
plt.show()