import scipy.signal as sgn
import numpy as np
import scipy.interpolate as interp
import time
import copy
import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def windowstats(stridearr):

    return np.mean(stridearr, axis=1), np.std(stridearr, axis=1)

def detect_peaks(data, thresh, iterations, valleys=True):

    count = 0
    y_mean = np.mean(data)
    y_std = np.std(data)
    final_mask = np.zeros(len(data), dtype = 'bool')
    comp = False
    while count < iterations and comp == False:
        print 'Test', y_std*thresh
        peakzero_max, = sgn.argrelmax(data, order = 30, mode = 'wrap')

        mask_max = np.zeros(len(data), dtype = 'bool')
        mask_max[peakzero_max] = True
        mask_max &= (data >= y_mean+thresh*y_std)

        temp_mask = mask_max.copy()

        if valleys == True:

            peakzero_min, = sgn.argrelmin(data, order = 30, mode = 'wrap')

            mask_min = np.zeros(len(data), dtype = 'bool')
            mask_min[peakzero_min] = True
            mask_min &= (data <= y_mean-thresh*y_std)

            temp_mask ^= mask_min
        
        peak_temp = np.where(temp_mask == True)[0]
        plt.plot(data)
        plt.plot(peak_temp, data[peak_temp], 'x')
        plt.show()
        print peak_temp
        le, re = width(data, peak_temp)
        exclude_mask = np.zeros(len(data), dtype = 'bool')
        for i in range(len(peak_temp)):
            exclude_mask[le[i]:re[i]] = True

        masked_array = np.ma.array(data, mask = exclude_mask)
        y_mean = np.ma.mean(masked_array)
        y_std = np.ma.std(masked_array)
        final_mask_temp = temp_mask.copy()
        comp = np.array_equal(final_mask, temp_mask)
        dim_temp = np.size(np.where(temp_mask == True)[0])
        
        if count != 0:
            dim_fin = np.size(np.where(final_mask == True)[0])
        else:
            dim_fin = 0
            ln_start = copy.copy(dim_temp)
        if ln_start > dim_temp:
            break
        else:
            if dim_temp >= dim_fin:
                final_mask = np.logical_or(final_mask, temp_mask)
            else:
                final_mask = temp_mask.copy()

        count += 1
        thresh *= 1.5

    return np.where(final_mask == True)[0]

def detect_valleys(data, thresh):

    y_mean = np.mean(data)
    y_std = np.std(data)

    peakzero_min, = sgn.argrelmin(data, order = 30, mode = 'wrap')

    mask_min = np.zeros(len(data), dtype = 'bool')
    mask_min[peakzero_min] = True
    mask_min &= (data <= y_mean-thresh*y_std)

    return np.where(mask_min == True)[0]
    
def width(data, peaks, loc=None):

    if loc==None:
        loc = np.mean(data)

    ledge = np.zeros(np.size(peaks), dtype=int)
    redge = np.zeros(np.size(peaks), dtype=int)
    for i in range(np.size(peaks)):

        if i == 0 and np.size(peaks) == 1:
            index1r = peaks[i]
            index2r = np.size(data)-1
            index1l = 0
            index2l = peaks[i]
        else:
            if i == 0:
                index1r = peaks[i]
                index2r = peaks[i+1]
                index1l = 0
                index2l = peaks[i]
            elif i != np.size(peaks) - 1:
                index1r = peaks[i]
                index2r = peaks[i+1]
                index1l = peaks[i-1]
                index2l = peaks[i]
            else:
                index1r = peaks[i]
                index2r = np.size(data)-1
                index1l = peaks[i-1]
                index2l = peaks[i]
        if peaks[i] == np.size(data)-1:
            rzero = np.size(data)-1
        else:
            f = interp.splrep(np.arange(len(data[index1r:index2r])), \
                                        data[index1r:index2r]-loc)
            rzero = interp.sproot(f, mest=10)
        if peaks[i] != 0:
            lzero_temp = np.zeros(1000)
            index1 = index1l
            count = 0
            print index1, index2l
            while len(lzero_temp) == 1000:
                f = interp.splrep(np.arange(len(data[index1:index2l])),\
                                                data[index1:index2l]-loc)
                lzero_temp = interp.sproot(f,mest=1000)
                if np.size(lzero_temp) <2:
                    print lzero_temp
                lzero = lzero_temp[-1]+index1
                index1 = int(np.floor(lzero_temp[-1]))+index1
            if np.size(lzero_temp) == 0:
                lzero = 0
        else:
            lzero = 0
        if isinstance(lzero, np.ndarray):
            fedge = lzero[-1]+peaks[i-1]
        else:
            fedge = lzero
        if isinstance(rzero, np.ndarray):
            if np.size(rzero) > 1:
                sedge = rzero[0]+peaks[i]
            else:
                sedge = rzero+peaks[i]
        else:
            sedge = rzero
        #print 'value', fedge, sedge, peaks[i]
        if peaks[i]>=fedge and peaks[i]<=sedge:
            ledge[i] = np.floor(fedge)
            redge[i] = np.ceil(sedge)
        elif peaks[i]<=fedge and peaks[i]<=sedge:
            ledge[i] = 0
            redge[i] = np.ceil(fedge)
        if peaks[i]>=fedge and peaks[i]>=sedge:
            ledge[i] = np.floor(fedge)
            redge[i] = np.size(data)-1

    return ledge, redge

def replace(data, ledge, redge):

    mask = np.zeros(len(data), dtype = bool)

    for i in range(len(ledge)):

        mask[ledge[i]:redge[i]+1] = True

    masked_array = np.ma.array(data, mask = mask)

    mean = np.ma.mean(masked_array)
    std = np.ma.std(masked_array)

    data[mask] = np.random.normal(mean, std, np.ma.count_masked(masked_array))

    return data

        
