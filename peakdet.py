import scipy.signal as sgn
import numpy as np
import scipy.interpolate as interp

def detect_peaks(data, thresh):

    y_mean = np.mean(data)
    y_std = np.std(data)

    peakzero_max, = sgn.argrelmax(data, order = 30, mode = 'wrap')

    mask_max = np.zeros(len(data), dtype = 'bool')
    mask_max[peakzero_max] = True
    mask_max &= (data >= y_mean+thresh*y_std)

    return np.where(mask_max == True)[0]

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
        print 'new cycle'

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
            while len(lzero_temp) == 1000:
                f = interp.splrep(np.arange(len(data[index1:index2l])),\
                                                data[index1:index2l]-loc)
                lzero_temp = interp.sproot(f,mest=1000)
                lzero = lzero_temp[-1]+index1
                index1 = int(np.floor(lzero_temp[-1]))+index1l
            if np.size(lzero_temp) == 0:
                lzero = 0
        else:
            lzero = 0
        if isinstance(lzero, np.ndarray):
            fedge = lzero[-1]+peaks[i-1]
        else:
            fedge = lzero
        if isinstance(rzero, np.ndarray):
            sedge = rzero[0]+peaks[i]           
        else:
            sedge = rzero
        print 'value', fedge, sedge, peaks[i]
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



        
