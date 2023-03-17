import numpy as np
import struct
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm

# Find all files with the pattern "Mass_100_Energy_*.bin" in the directory
file_pattern = '/media/integra/acd1b1cd-2cf9-4837-9619-4230c7a64b1d/Documents/2021_ILL/Analysis/Det_97_Binary/230313/Mass_100_Energy_*.bin'

# Initialize an empty list to hold the data
all_data = []

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

    data = np.array(data[9:])
    num_elements = number_of_signals * 2000
    data = data[:num_elements].reshape(number_of_signals, 2000)

    # Add the data to the list
    all_data.append(data)

# Stack all the data together into a single 3D array
y = np.vstack(all_data)
x = np.arange(0, 2000)

# Compute the maximum value for each signal in y
max_values = np.max(y[0:1352], axis=1)

# Compute the maximum value among all elements in y
global_max_value = np.max(y[0:1352,:])

# Create a histogram of the max values
num_bins = 50
plt.hist(max_values, bins=num_bins)
plt.xlabel('Max value')
plt.ylabel('Frequency')
plt.title('Histogram of max values in signals')
plt.show()

print("Max value among all elements in y:", global_max_value)

def return_beg(y):
    dxdy = np.gradient(y)
    bool_arr = np.array([dxdy[i+1]-dxdy[i] > 0 for i in range(len(dxdy)-1)], dtype = bool)  
    consecutive_true = np.convolve(bool_arr, np.ones(6), mode='valid') == 6
    start_smp = (np.argmax(consecutive_true) + 6) if any(consecutive_true) else None
    return (490 + start_smp)-7

def return_end(x,y,smp2extrpl_x,smp2extrpl_y):
    slope, intercept = np.polyfit(x,y,1)
    extrapolation = np.array(intercept + slope*smp2extrpl_x)
    smp_inter = np.argmin(np.abs(smp2extrpl_y - extrapolation))
    return smp_inter+500

smp_start = np.empty(1352, dtype=int)
smp_stop = np.empty(1352, dtype=int)

for i in range(1352):
    smp_start[i] = return_beg(y[i,490:520])
    smp_stop[i] = return_end(x[1800:2000], y[i,1800:2000], x[500:600], y[i,500:600])

rise_time = np.array(smp_stop - smp_start)
print("smp_start: ", smp_start[5])
print("smp_stop: ", smp_stop[5])
print(rise_time[0:10])

# Define the Gaussian function
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * stddev ** 2))

num_bins = int(np.max(rise_time) - np.min(rise_time))
mu,sigma = norm.fit(rise_time)
p0 = [1, mu, sigma]
coeff, cov = curve_fit(gaussian,  np.arange(len(rise_time)), rise_time, p0=p0)

# Extract the errors of mu and sigma from the diagonal of the covariance matrix
mu_err = np.sqrt(cov[1, 1])
sigma_err = np.sqrt(cov[2, 2])

print("mu_err: ", mu_err)
print("sigma_err: ", sigma_err)


# # Plot the first signal
# print("Max value of first signal:", np.max(y))
# plt.plot(x, y[0], color="red")
# plt.show()
