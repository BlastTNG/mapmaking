import dirfile_functions as df
import numpy as np
import matplotlib.pyplot as plt
import pygetdata as gd
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d
import pandas as pd
from scipy import signal
roach_path = 'roach_data/'
ancillary_path = 'xy_stage/'
old_path = '2012_data/'
bolo_path = '2012_data/bolo_data/'

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

#stolen from stack overflow to test
def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

#### below here are newer useful functions, above is preserved for records
def loadArbData(dirfile, file, file_type):
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    #print vectors
    for i in range (len(vectors)):
        if vectors[i] == file:
            if file_type == 'u16':
                values = d.getdata(file, gd.UINT16, num_frames = d.nframes)
            if file_type == 'u32':
                values = d.getdata(file, gd.UINT32, num_frames = d.nframes)
            if file_type == 's32':
                values = d.getdata(file, gd.INT32, num_frames = d.nframes)
    return np.asarray(values)


def arb_binning(ra, dec, bolo_data, pixel_size):
    #pixel size is in degrees so pass ra and dec in degrees
    # returns bolo data binned for map processing
    n_ra = np.int(np.ceil((ra.max()-ra.min())/pixel_size))+1
    n_dec = np.int(np.ceil((dec.max()-dec.min())/pixel_size))+1

    ra_bins = np.amin(ra)-pixel_size/2.+np.arange(n_ra+1)*pixel_size
    dec_bins = np.amin(dec)-pixel_size/2.+np.arange(n_dec+1)*pixel_size

    i_ind = np.digitize(ra,ra_bins)
    j_ind = np.digitize(dec, dec_bins)

    bolo_array = np.zeros((len(bolo_data),3))

    bolo_array[:,0] = bolo_data
    bolo_array[:,1] = i_ind
    bolo_array[:,2] = j_ind
    size = [n_ra, n_dec]
    return bolo_array, size

def ra_to_deg(ra):
    ra = ra*15
    return ra


# for posterity, the original processing function, use arb binning now
def arb_map_gen(ra, dec, bolo_data, pixel_size):
    #note that the bolo data ra and dec data passed must all be the same length
    # process the RA and DEC data into degrees first and interpolate.
    #pixel size is in degrees so pass ra and dec in degrees
    n_ra = np.int(np.ceil((ra.max()-ra.min())/pixel_size))
    n_dec = np.int(np.ceil((dec.max()-dec.min())/pixel_size))
    ra_bins, dec_bins = [], []
    ijbol_list = []
    size = [n_ra, n_dec]
    print n_ra
    print n_dec
    for i in range(0, n_ra):
        ra_bins.append(ra.min()+i*pixel_size)
    for j in range(0, n_dec):
        dec_bins.append(dec.min()+j*pixel_size)
    for n in range(0, len(bolo_data)):
        ra_test = np.abs(ra_bins-ra[n])
        dec_test = np.abs(dec_bins-dec[n])
        i_ind = np.where(ra_test == np.min(ra_test))
        j_ind = np.where(dec_test == np.min(dec_test))
        ijbol_list.append([bolo_data[n], i_ind[0][0], j_ind[0][0]])
    return ijbol_list, size