import numpy as np
import pygetdata as gd
import matplotlib.pyplot as plt
roach_path = 'roach_data/'
otherdata_path = 'xy_stage/'

def avg_arrays(how_long, array):
    counter = 0
    current_sum = 0.
    out_array = []
    for i in range(len(array)):
        current_sum += array[i]
        counter += 1
        if counter == how_long:
            counter = 0
            out_array.append(float(current_sum/how_long))
            current_sum = 0.
    return out_array
    
def avg_arrays(how_long, array):
    counter = 0
    current_sum = 0.
    out_array = []
    for i in range(len(array)):
        current_sum += array[i]
        counter += 1
        print counter, current_sum
        if counter == how_long:
            counter = 0
            print current_sum/how_long
            out_array.append(float(current_sum/how_long))
            current_sum = 0.
    return out_array

def IQtoMag(i_vals, q_vals):
    mag = []
    for i in range(len(i_vals)):
        mag.append(np.sqrt(np.square(i_vals[i])+np.square(q_vals[i])))
    return mag

def data_generator(x_array, y_array, value_array):
    out_array = []
    for i in range(len(x_array)):
        out_array.append([x_array[i], y_array[i], value_array[i]])
    return out_array