import numpy as np
import os
import scipy.signal as sgn
from scipy.interpolate import interp1d
import pygetdata as gd
from astropy import wcs, coordinates
from astropy.convolution import Gaussian2DKernel, convolve
import ConfigParser
import copy

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

class dataload():

    def __init__(self, filepath_map):

        self.filepath_map = filepath_map

    def conversion_type(self, file_type):

        if file_type == 'u16':
            gdtype = gd.UINT16
        elif file_type == 'u32':
            gdtype = gd.UINT32
        elif file_type == 's32':
            gdtype = gd.INT32

        return gdtype

    def map_param(self):

        model = ConfigParser.ConfigParser()
        model.read(self.filepath_map)
        sections = model.sections()

        for section in sections:

            if section.lower() == 'file repository':
                dirfile = model.get(section,'DIRFILE').split('#')[0].strip()
                detfile = model.get(section,'DETFILE').split('#')[0].strip()
                coor1type = model.get(section,'COOR1TYPE').split('#')[0].strip().lower()
                coor2type = model.get(section,'COOR2TYPE').split('#')[0].strip().lower()
                pointingfile = model.get(section,'POINTINGOFF').split('#')[0].strip() 
                dettable = model.get(section, 'DETTABLE').split('#')[0].strip()
                detpath = model.get(section, 'DETPATH').split('#')[0].strip()
                coorpath = model.get(section, 'COORPATH').split('#')[0].strip()
                
                detlist = np.loadtxt(detpath+detfile, unpack = 'True', dtype = str)


            elif section.lower() == 'map parameters':
                ctype = model.get(section,'CTYPE').split('#')[0].strip()   
                crpix = model.get(section,'CRPIX').split('#')[0].strip()
                crpix = np.array(crpix.split(',')).astype(float)
                crdelt = model.get(section,'CRDELT').split('#')[0].strip()
                crdelt = np.array(crdelt.split(',')).astype(float)
                crval = model.get(section,'CRVAL').split('#')[0].strip()
                crval = np.array(crval.split(',')).astype(float)
                conv = model.get(section,'CONV').split('#')[0].strip()
                stdev = float(model.get(section,'STDEV').split('#')[0])

                #Conversion of stdev from arcsec to pixel
                if conv.lower() != 'na':
                    stdev = float(stdev)/(np.mean(crdelt)*3600.)

            elif section.lower() == 'experiment parameters':

                detfreq = float(model.get(section, 'DETFREQ').split('#')[0])
                acsfreq = float(model.get(section, 'ACSFREQ').split('#')[0])
                det_dir_conv = model.get(section,'DET_DIR_CONV').split('#')[0].strip()
                det_dir_conv = np.array(det_dir_conv.split(',')).astype(float)
                coor1_dir_conv = model.get(section,'COOR1_DIR_CONV').split('#')[0].strip()
                coor1_dir_conv = np.array(coor1_dir_conv.split(',')).astype(float)
                coor2_dir_conv = model.get(section,'COOR2_DIR_CONV').split('#')[0].strip()
                coor2_dir_conv = np.array(coor2_dir_conv.split(',')).astype(float)
                frames = model.get(section,'FRAMES').split('#')[0].strip()
                frames = np.array(frames.split(',')).astype(int)
                det_samp_frame = float(model.get(section, 'DET_SAMP_FRAME').split('#')[0])
                acs_samp_frame = float(model.get(section, 'ACS_SAMP_FRAME').split('#')[0])
                det_file_type = model.get(section,'DET_FILE_TYPE').split('#')[0].strip()
                coor1_file_type = model.get(section,'COOR1_FILE_TYPE').split('#')[0].strip()
                coor2_file_type = model.get(section,'COOR2_FILE_TYPE').split('#')[0].strip()

        return dirfile, detlist, coor1type, coor2type, pointingfile, dettable, \
               detpath, coorpath, ctype, crpix, crdelt, crval, conv, stdev, \
               detfreq, acsfreq, det_dir_conv, coor1_dir_conv, coor2_dir_conv, frames, \
               det_samp_frame, acs_samp_frame, det_file_type, coor1_file_type, coor2_file_type
        
    def loaddata(self, filepath, file, file_type):
        if np.size(file) == 1: 
            d = gd.dirfile(filepath, gd.RDWR|gd.UNENCODED)
            vectors = d.field_list()
            print file
            for i in range (len(vectors)):
                if vectors[i] == file:
                    gdtype = self.conversion_type(file_type)
                    values = d.getdata(file, gdtype, num_frames = d.nframes)
            return np.asarray(values, dtype = 'int')
        else:
            d = gd.dirfile(filepath[0], gd.RDWR|gd.UNENCODED)
            vectors = d.field_list()
            len_det = len(d.getdata(vectors[detlist[0]], gd.UINT16, num_frames = d.nframes))
            values = np.zeros((len(filepath), len_det))

            for i in range(len(filepath)):
                d = gd.dirfile(filepath[i], gd.RDWR|gd.UNENCODED)
                vectors = d.field_list()
                for j in range(len(vectors)):
                    if vectors[j] == file[i]:
                        values[i,:] = np.asarray(d.getdata(vectors[file[i]], \
                                                 gdtype_det,num_frames = d.nframes))
                
            return values
    
    def loadfulldata(self, filepath, detlist, det_type, coord, coord_type1, coord_type2):
        d = gd.dirfile(filepath, gd.RDWR|gd.UNENCODED)
        vectors = d.field_list()
        len_det = len(d.getdata(vectors[detlist[0]], gd.UINT16, num_frames = d.nframes))
        data_value = np.zeros((len(detlist), len_det))

        gdtype_det = self.conversion_type(det_type)
        for i in range(len(detlist)):
            data_value[i,:] = np.asarray(d.getdata(vectors[detlist[i]], gdtype_det,\
                                         num_frames = d.nframes))

        gdtype_coord1 = self.conversion_type(coord_type1)
        gdtype_coord2 = self.conversion_type(coord_type2)
        if coord == 'RADEC':
            ra = np.asarray(d.getdata('ra', gdtype_coord1, num_frames = d.nframes), dtype = 'int')
            dec = np.asarray(d.getdata('dec', gdtype_coord2, num_frames = d.nframes), dtype = 'int')

            return data_value, ra, dec

        elif coord == 'ELAZ':

            el = np.asarray(d.getdata('el', gdtype_coord1, num_frames = d.nframes), dtype = 'int')
            az = np.asarray(d.getdata('az', gdtype_coord2, num_frames = d.nframes), dtype = 'int')

            return data_value, el, az

    def convert_dirfile(self, data_value, value1, value2=0.):

        data_value = value1*data_value+value2   

        return data_value

    def frame_zoom(self, data, sample_frame, fs, fps):

        frames = fps.copy()
        
        if sample_frame % 2 == 0.:
            frames[0] = fps[0]*sample_frame
            frames[1] = fps[1]*sample_frame+1
        else:
            frames[0] = fps[0]*sample_frame
            frames[1] = fps[1]*sample_frame+1

        if len(np.shape(data)) == 1:
            time = np.arange(len(data))/fs
            time = time[frames[0]:frames[1]]
            return time, data[frames[0]:frames[1]]
        else:
            time = np.arange(len(data[0, :]))/fs
            time = time[frames[0]:frames[1]]
            return  time, data[:,frames[0]:frames[1]]

    def coord_int(self, coord1, coord2, time_acs, time_det):

        coord1_int = interp1d(time_acs, coord1, kind='linear')
        coord2_int = interp1d(time_acs, coord2, kind= 'linear')

        return coord1_int(time_det), coord2_int(time_det)

class despike():

    def __init__(self, data, peakind, height, width, ledge, redge):

        self.data = data
        self.peakind = peakind
        self.height = height
        self.width = width
        self.ledge = ledge
        self.redge = redge

    def findpeak(self, hthres=0, pthres=0):

        '''
        hthresh and pthres are measured in how many std the height (or the prominence) 
        of the peak is computed. The height of the peak is computed with respect to 
        the mean of the signal        
        '''
        index = np.ones(1)
        index_final = np.array([], dtype = 'int')
        ledge = np.array([], dtype = 'int')
        redge = np.array([], dtype = 'int')
        count = 0
        while len(index) > 0:

            print index_final
            print ledge, redge

            if np.size(index_final) != 0:
                for i in range(len(index_final)):
                    data_masked[ledge[i]:redge[i]] = (data_masked[ledge[i]]+data_masked[redge[i]])/2.
            else:
                data_masked = self.data.copy()
            y_std = np.std(data_masked)
            y_mean = np.mean(data_masked)

            print y_std, y_mean
            
            if count >=2:
                sys.exit()
            if hthres != 0 and pthres == 0:
                index, param = sgn.find_peaks(data_masked, height = y_mean + hthres*y_std, distance=100)
            elif pthres != 0 and hthres == 0:
                index, param = sgn.find_peaks(data_masked, prominence = pthres*y_std)
            elif hthres != 0 and pthres != 0:
                index, param = sgn.find_peaks(data_masked, height = y_mean + hthres*y_std, \
                                              prominence = pthres*y_std)

            plt.plot(data_masked)
            plt.plot(index, self.data[index], 'x')
            plt.show()
            
            ledget = sgn.peak_widths(data_masked,index)[2]
            redget = sgn.peak_widths(data_masked,index)[3]

            index_final = np.append(index_final, index.astype(int))
            ledge = np.append(ledge, np.floor(ledget).astype(int))
            redge = np.append(redge, np.ceil(redget).astype(int))
            count += 1

        self.peakind = index_final.copy()

        return self.peakind

    def peak_width(self):
        
        param = sgn.peak_widths(self.data,self.peakind)
        self.width = param[0].copy()
        self.ledge = param[2].copy()
        self.redge = param[3].copy()

        return self.width, self.ledge, self.redge

    def replace_peak(self, hthres=5, pthres = 0):

        x_inter = np.array([], dtype = 'int')

        if np.size(self.peakind) == 0:
            self.findpeak(hthres=hthres, pthres=pthres)
        if np.size(self.width) == 0:
            self.peak_width()

        for i in range(0, len(self.peakind)):

            width = int(np.ceil(self.width[i]))
            if width <= 13:
                interval = 25
            elif width > 13 and width < 40:
                interval = width*2
            else:
                interval = width*3

            left_edge = int(np.floor(self.ledge[i]))
            right_edge = int(np.ceil(self.redge[i]))

            x_inter = np.append(x_inter, np.arange(left_edge, right_edge))
            self.data[left_edge:right_edge] = (self.data[left_edge]+\
                                                 self.data[right_edge])/2.

        final_mean = np.mean(self.data)
        final_std = np.std(self.data)
        final_var = np.var(self.data)

        p_stat = np.abs(final_mean/final_var-1.)

        if p_stat <=1e-2:
            '''
            This means that the variance and the mean are approximately the 
            same, so the distribution is Poissonian.
            '''
            mu = (final_mean+final_var)/2.
            y_sub = np.random.poisson(mu, len(x_inter))
        else:
            y_sub = np.random.normal(final_mean, final_std, len(x_inter))
        if np.size(y_sub) > 0:
            self.data[x_inter] = y_sub

        return self.data

class filterdata():

    def __init__(self, data, cutoff, fs):
        
        '''
        fs: sample frequency
        cutoff: cutoff frequency
        '''

        self.data = data
        self.cutoff = cutoff
        self.fs = fs
    
    def highpass(self, order):
        
        nyq = 0.5*self.fs
        normal_cutoff = self.cutoff / nyq
        b, a = sgn.butter(order, normal_cutoff, btype='highpass', analog=False)
        return b, a

    def butter_highpass_filter(self, order=5):
        b, a = self.highpass(order)
        filterdata = sgn.lfilter(b, a, self.data)
        return filterdata

    def cosine_filter(self, f):
        if f < .5*self.cutoff:
            return 0
        elif 0.5*self.cutoff < f  and f < self.cutoff:
            return 0.5-0.5*np.cos(np.pi*(f-0.5*self.cutoff)*(self.cutoff-0.5*self.cutoff)**-1)
        elif f > self.cutoff:
            return 1
    
    def fft_filter(self, window):

        if window == True:
            window = np.hanning(len(self.data))

            fft_data = np.fft.rfft(self.data*window)
        else:
            fft_data = np.fft.rfft(self.data)
        fft_frequency = np.fft.rfftfreq(np.size(self.data), 1/self.fs)

        vect = np.vectorize(self.cosine_filter)

        filtereddata = vect(fft_frequency)*fft_data

        return filtereddata

    def ifft_filter(self, window):
        
        ifft_data = np.fft.irfft(self.fft_filter(window=window), len(self.data))

        return ifft_data

class rotate():

    '''
    Pitch is a rotation around the x axis
    Yaw is a rotation around the z axis
    Roll is a rotation around the y axis
    To go from inertial to gondola rotate zxy (first yaw, then pitch and then roll)
    '''

    def __init__(self, yaw, pitch, roll):

        self.yaw = yaw
        self.pitch = pitch
        self.roll = roll

    def rotmatrix(self, yaw_mat = None, roll_mat = None, pitch_mat=None):
        if yaw_mat == None:
            yaw_mat = self.yaw.copy()
            roll_mat = self.roll.copy()
            pitch_mat = self.pitch.copy()
        yawMatrix = np.matrix([[np.cos(yaw_mat), -np.sin(yaw_mat), 0], \
                               [np.sin(yaw_mat), np.cos(yaw_mat), 0], \
                               [0, 0, 1]])

        rollMatrix = np.matrix([[np.cos(roll_mat), 0, -np.sin(roll_mat)],\
                                [0, 1, 0],\
                                [np.sin(roll_mat), 0, np.cos(roll_mat)]])

        pitchMatrix = np.matrix([[1, 0, 0],\
                                 [0, np.cos(pitch_mat), -np.sin(pitch_mat)],\
                                 [0, np.sin(pitch_mat), np.cos(pitch_mat)]]) 

        return pitchMatrix, rollMatrix, yawMatrix

    def offset_mat(self, yaw_off, pitch_off, roll_off, rot_mat=np.diag(np.ones(3))):

        pitch_off_mat = rotmatrix(yaw_mat = yaw_off, roll_mat = roll_off, pitch = pitch_off)[0]
        roll_off_mat = rotmatrix(yaw_mat = yaw_off, roll_mat = roll_off, pitch = pitch_off)[1]
        yaw_off_mat = rotmatrix(yaw_mat = yaw_off, roll_mat = roll_off, pitch = pitch_off)[2]

        rot1 = np.matmul(yaw_off_mat, rot_mat)
        rot2 = np.matmul(pitch_off_mat, rot1)
        rot3 = np.matmul(roll_off_mat, rot2)

        return rot3
    
    def offset_angle(self, yaw_off=0., pitch_off=0., roll_off=0., rot_mat=0.):

        if np.size(yaw_off) == 1:
            if np.greater(yaw_off,0.)==True or np.greater(roll_off,0.)==True or \
               np.greater(pitch_off,0.)==True:

                rot_matrix = offset_mat(yaw_off, pitch_off, roll_off)
        else:
            if np.any(np.greater(yaw_off,0.))==True or np.any(np.greater(pitch_off,0.))==True or \
                np.any(np.greater(roll_off,0.))==True:

                matrix = np.diag(np.ones(3))
                for i in len(yaw_off):
                    rot_matrix = offset_mat(yaw_off[i], pitch_off[i], roll_off[i], rot_mat=matrix)
                    matrix = rot_matrix.copy()

        if np.size(rot_mat) >= 3:

            rot_matrix = rot_mat

        pitch_off_final = np.arctan2(rot_matrix[1,2],np.sqrt(rot_matrix[1,0]**2+rot_matrix[1,1]**2))
        roll_off_final = np.arctan2(rot_matrix[0,2],rot_matrix[2,2])
        yaw_off_final = np.arctan2(rot_matrix[1,0],rot_matrix[1,1])

        return pitch_off_final, roll_off_final, yaw_off_final

    def finalcoord(self, yaw_off=0., pitch_off=0., roll_off=0.):

        cr = np.cos(self.roll)
        sr = np.sin(self.roll)
        cp = np.cos(self.pitch)

        self.yaw = self.yaw+2*np.arcsin(np.sin((yaw_off*cr+pitch_off*sr)/2.)/cp)
        self.roll = self.roll
        self.pitch = self.pitch+(-yaw_off*sr+pitch_off*cr)

        return self.pitch, self.roll, self.yaw
        
class detector():

    def __init__(self, data, responsivity, grid):

        self.data = data
        self.responsivity = responsivity
        self.grid = grid

    def calibrate(self):

        return self.data*self.responsivity

    def polangle(self, roll, hwp_angle):

        return self.grid-2*hwp_angle+roll

    def polyfit(self, n):

        x = np.arange(len(self.data))
        y = self.data.copy()

        poly = np.polyfit(x, y, n)
        p = np.poly1d(poly)

        return p(x)

class wcs_world():

    def __init__(self, ctype, crpix, crdelt, crval):

        self.ctype = ctype
        self.crdelt = crdelt
        self.crpix = crpix
        self.crval = crval

    def world(self, coord):
        
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = self.crpix
        w.wcs.cdelt = self.crdelt
        w.wcs.crval = self.crval
        if self.ctype.lower() == 'radec':
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        elif self.ctype.lower() == 'elaz':
            w.wcs.ctype = ["TLON-ARC", "TLAT-ARC"]

        world = w.wcs_world2pix(coord, 1)

        return world, w

class mapmaking(object):

    def __init__(self, data, weight, polangle, number, pixelmap):

        self.data = data
        self.weight = weight
        self.polangle = polangle
        self.number = number
        self.pixelmap = pixelmap

    def map_param(self, value = None, sigma=None, angle=None):

        if value == None:
            value = self.data.copy()
            if np.size(self.weight) > 1:
                sigma = self.weight.copy()
            else:
                sigma = copy.copy(self.weight)
            if np.size(self.polangle) > 1:
                angle = self.polangle.copy()
            else:
                angle = copy.copy(self.polangle)*np.ones(np.size(value))

        '''
        sigma is the inverse of the sqared white noise value, so it is 1/n**2
        '''
        x_map = self.pixelmap[:,0]   #RA 
        y_map = self.pixelmap[:,1]   #DEC

        x_map = np.round(x_map+np.abs(np.amin(x_map)))
        y_map = np.round(y_map+np.abs(np.amin(y_map)))

        x_len = np.amax(x_map)-np.amin(x_map)+1
        
        param = x_map+y_map*x_len
        param = param.astype(int)

        flux = value
        cos = np.cos(2.*angle)
        sin = np.sin(2.*angle)

        print type(param)
        print np.amax(x_map), np.amin(x_map)
        print np.amax(y_map), np.amin(y_map)
        print np.amax(param), np.amin(param)

        I_est_flat = np.bincount(param, weights=flux)*sigma
        Q_est_flat = np.bincount(param, weights=flux*cos)
        U_est_flat = np.bincount(param, weights=flux*sin)

        print len(I_est_flat)
        print 'count', len(np.bincount(param))

        N_hits_flat = 0.5*np.bincount(param)*sigma
        print len(N_hits_flat)
        c_flat = np.bincount(param, weights=0.5*cos)*sigma
        c2_flat = np.bincount(param, weights=0.5*cos**2)*sigma
        s_flat = np.bincount(param, weights=0.5*sin)*sigma
        s2_flat = N_hits_flat-c2_flat
        m_flat = np.bincount(param, weights=0.5*cos*sin)*sigma

        Delta = c_flat**2*(c2_flat-N_hits_flat)+2*s_flat*c_flat*m_flat-c2_flat*s_flat**2-\
                N_hits_flat*(c2_flat**2+m_flat**2-c2_flat*N_hits_flat)
        A = -(c2_flat**2+m_flat**2-c2_flat*N_hits_flat)
        B = c_flat*(c2_flat-N_hits_flat)+s_flat*m_flat
        C = c_flat*m_flat-s_flat*c2_flat
        D = -((c2_flat-N_hits_flat)*N_hits_flat+s_flat**2)
        E = c_flat*s_flat-m_flat*N_hits_flat
        F = c2_flat*N_hits_flat-c_flat**2

        return I_est_flat, Q_est_flat, U_est_flat, N_hits_flat, Delta, A, B, C, D, E, F, param

    def map_singledetector_Ionly(self, value=None, sigma=None, angle=None):

        value =self.map_param(value=value, sigma=sigma, angle=angle)

        I = np.zeros(len(value[0]))
        I[np.nonzero(value[0])] = value[0][np.nonzero(value[0])]/value[3][np.nonzero(value[3])]

        x_len = np.amax(self.pixelmap[:,0])-np.amin(self.pixelmap[:,0])
        y_len = np.amax(self.pixelmap[:,1])-np.amin(self.pixelmap[:,1])

        if len(I) < (x_len+1)*(y_len+1):
            valmax = (x_len+1)*(y_len+1)
            pmax = np.amax(value[-1])
            I_fin = 0.*np.arange(pmax+1, valmax)
            
            I = np.append(I, I_fin)

        # index, = np.where(I==0)
        # if np.size(index) >=1 :
        #     I[index] = 0.

        I_pixel = np.reshape(I, (y_len+1,x_len+1))

        return I_pixel, 1

    def map_multidetectors_Ionly(self):

        for i in range(self.number):
            mapvalues = map_singledetector_Ionly(value=self.data[i],sigma=self.weight[i],\
                                                 angle=self.polangle[i])

            I_pixel += mapvalues[0]

        return I_pixel, 2 

    def map_singledetector(self, value=None, sigma=None, angle=None):


        I_est_flat, Q_est_flat, U_est_flat, N_hits_flat, Delta, \
                    A, B, C, D, E, F, param = self.map_param(value=value, sigma=sigma,angle=angle)

        I_pixel_flat = (A*I_est_flat+B*Q_est_flat+C*U_est_flat)/Delta
        Q_pixel_flat = (B*I_est_flat+D*Q_est_flat+E*U_est_flat)/Delta
        U_pixel_flat = (C*I_est_flat+E*Q_est_flat+F*U_est_flat)/Delta

        I_pixel = np.reshape(I_pixel_flat, (len(self.pixelmap[:,0]),len(self.pixelmap[:,1])))
        Q_pixel = np.reshape(Q_pixel_flat, (len(self.pixelmap[:,0]),len(self.pixelmap[:,1])))
        U_pixel = np.reshape(U_pixel_flat, (len(self.pixelmap[:,0]),len(self.pixelmap[:,1])))

        return I_pixel, Q_pixel, U_pixel, 1

    def map_multidetectors(self):

        print 'This method gives realistic results only if the detector are calibrated'

        for i in range(self.number):

            mapvalues = map_singledetector(value=self.data[i],sigma=self.weight[i],\
                                           angle=self.polangle[i])
            
            I_map += mapvalues[0]
            Q_map += mapvalues[1]
            U_map += mapvalues[2]

        return I_map, Q_map, U_map, 2

    def convolution(self, std, map_value):

        kernel = Gaussian2DKernel(stddev=std)

        convolved_map = convolve(map_value, kernel)

        return convolved_map
        
class computeoffset():

    def __init__(self, data, angX_center, angY_center):

        self.data = data
        self.angX_center = angX_center 
        self.angY_center = angY_center

    def centroid(self, threshold=0.275):

        maxval = np.max(self.data)
        minval = np.min(self.data)
        y_max, x_max = np.where(self.data == maxval)

        lt_inds = np.where(self.data < threshold*maxval)
        gt_inds = np.where(self.data > threshold*maxval)

        weight = np.zeros((self.data.shape[1], self.data.shape[0]))
        weight[gt_inds] = 1.
        a = self.data[gt_inds]
        flux = np.sum(a)
        x_range = np.arange(0, self.data.shape[0])
        y_range = np.arange(0, self.data.shape[1])

        xx, yy = np.meshgrid(x_range, y_range)

        x_c = np.sum(xx*weight*self.data)/flux
        y_c = np.sum(yy*weight*self.data)/flux

        return np.rint(x_c), np.rint(y_c)

    def offset(self, wcs_trans, frame, threshold=0.275, return_pixel=True, altitude=0., lon=0., lat=0.):

        x_c, y_c = self.centroid(threshold=threshold)

        coord_centre = coordinates.Skycoord(self.angX_center, self.angY_center, frame, unit='deg')

        if return_pixel == True:
        
            x_map, y_map = wcs.utils.skycoord_to_pixel(coord, wcs_trans)

            x_off = x_map-x_c
            y_off = y_map-y_c

            return x_off, y_off
        
        else:
            coord = wcs.utils.pixel_to_skycoord(x_c, y_c, wcs_trans)
                        
            offset_angle = coord.transform_to(coord_centre.skyoffset_frame())

            if frame.lower() == 'elaz':

                return offset_angle.az, offset_angle.alt

            else:

                return offset_angle.ra, offset_angle.dec




        






        
