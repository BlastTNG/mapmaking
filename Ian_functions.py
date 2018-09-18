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

def fancier_avg(array, length):
    nsamples = float(len(array))/float(length)
    counter = 0.
    counter_limit = np.floor(nsamples)
    run_sum = 0.
    offset_check = float(np.mod(nsamples,1.))
    i_offset = 0.
    f_offset = float(np.mod(nsamples,1.))
    avg_array = []
    for i in range(len(array)):
        if counter < counter_limit:
            run_sum += array[i]
            #print counter, i
            counter += 1.
        else:
            #print counter, i, 'split value'
            run_sum += f_offset*array[i]
            avg_array.append(run_sum/nsamples)
            i_offset = 1.-f_offset
            run_sum = i_offset*array[i]
            if i_offset > offset_check:
                f_offset = 1.+offset_check-i_offset
                counter_limit = np.floor(nsamples)-1
                #print i_offset, f_offset, counter_limit, 'was bigger'
            else:
                f_offset = offset_check-i_offset
                counter_limit = np.floor(nsamples)
                #print i_offset, f_offset, counter_limit, 'was smaller'
            counter = 0.
    return avg_array

def IQtoPhase(Iarray, Qarray):
    out_array = []
    for i in range(len(Iarray)):
        if Iarray[i] != 0:
            out_array.append(np.arctan(Qarray[i]/Iarray[i]))
        else:
            out_array.append(np.pi/2.)
    return out_array