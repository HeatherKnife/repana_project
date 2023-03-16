import struct
import numpy as np
import matplotlib.pyplot as plt

testsite_array = []

#Read the file
with open('/media/integra/acd1b1cd-2cf9-4837-9619-4230c7a64b1d/Documents/2021_ILL/Analysis/Det_97_Binary/230313/Mass_100_Energy_*.bin', 'br') as f:
    buffer = f.read()
    f.close()
    print("Length of buffer is %d" % len(buffer))

#Find the size
size = str(int(len(buffer) / 8)) + 'd'

#Unpack the data
data = struct.unpack(size, buffer)

#test the heading
testsite_array.append(data[0:9])
print(testsite_array)

#Find total nbr of signals on each file
number_of_signals = int((len(buffer) - 9 * 8) / (8 * 2000))
y = np.array([data[i * 2000 + 9:(i + 1) * 2000 + 9] for i in range(number_of_signals)])
x = np.arange(0, 2000)

# Compute the maximum value for each signal in y
max_values = np.max(y, axis=1)

# Compute the maximum value among all elements in y
global_max_value = np.max(y)

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

smp_start = np.empty(number_of_signals, dtype=int)
smp_stop = np.empty(number_of_signals, dtype=int)

for i in range(number_of_signals):
    smp_start[i] = return_beg(y[i,490:520])
    smp_stop[i] = return_end(x[1800:2000], y[i,1800:2000], x[500:600], y[i,500:600])

rise_time = np.array(smp_stop - smp_start)
print("smp_start: ", smp_start[5])
print("smp_stop: ", smp_stop[5])
print(rise_time[0:10])

# Create a histogram of the rise time
num_bins = int((np.max(rise_time) - np.min(rise_time)))
plt.hist(rise_time, bins=num_bins)
plt.xlabel('Rise time')
plt.ylabel('Frequency')
plt.title('Histogram of the signals rise time')
plt.show()


# # Plot the first signal
# print("Max value of first signal:", np.max(y))
# plt.plot(x, y[0], color="red")
# plt.show()


