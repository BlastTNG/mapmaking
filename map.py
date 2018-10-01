import numpy as np 
import scipy.signal as sgn

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
        while len(index) > 0:

            mask = np.ones(len(self.data), dtype = 'Bool')
            mask[index_final] = False
            data_masked = self.data[mask]
            y_std = np.ma.std(data_masked)
            y_mean = np.ma.mean(data_masked)
            if hthres != 0 and pthres == 0:
                index, param = sgn.find_peaks(data_masked, height = y_mean + hthres*y_std)
            elif pthres != 0 and hthres == 0:
                index, param = sgn.find_peaks(data_masked, prominence = pthres*y_std)
            elif hthres != 0 and pthres != 0:
                index, param = sgn.find_peaks(data_masked, height = y_mean + hthres*y_std, \
                                              prominence = pthres*y_std)

            index_final = np.append(index_final, index)

        self.peakind = index_final.copy()

        return self.peakind

    def peak_width(self):
        
        param = sgn.peak_widths(self.data,self.peakind)
        self.width = param[0].copy()
        self.ledge = param[2].copy()
        self.redge = param[3].copy()

        return self.width, self.ledge, self.redge

    def replace_peak(self):

        x_inter = np.array([], dtype = 'int')

        if len(self.peakind) == 0:
            self.findpeak(hthres=5)
        if len(self.width) == 0:
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

        self.data[x_inter] = y_sub

        return self.data

class filterdata():

    def __init__(self, data, cutoff, fs):

        self.data = data
        self.cutoff = cutoff
        self.fs = fs
    
    def highpass(self, order):
        
        nyq = 0.5*self.fs
        normal_cutoff = self.cutoff / nyq
        b, a = sgn.butter(order, normal_cutoff, btype='highpass', analog=False)
        return b, a

    def butter_lowpass_filter(self, order=5):
        b, a = highpass(order)
        self.data = sgn.lfilter(b, a, self.data)
        return self.data

class mapmaking():

    def __init__(self, data, weight, polangle, number):

        self.data = data
        self.weight = weight
        self.polangle = polangle
        self.number = number

    def map_singledetector(self, value=self.data, sigma=self.weight, angle=self.polangle):

        param = i*len(i)+j

        flux = value.copy()/sigma
        cos = np.cos(2.*angle.copy())
        sin = np.sin(2.*angle.copy())

        I_est_flat = np.bincount(param, weight=flux)
        Q_est_flat = np.bincount(param, weight=flux*cos)
        U_est_flat = np.bincount(param, weight=flux*sin)

        N_hits_flat = np.bincount(param)/sigma
        c_flat = np.bincount(param, weight=0.5*cos)/sigma
        c2_flat = np.bincount(param, weight=0.5*cos**2)/sigma
        s_flat = np.bincount(param, weight=0.5*sin)/sigma
        s2_flat = N_hits_flat-c2_flat
        m_flat = np.bincount(param, weight=0.5*cos*sin)/sigma

        Delta = c_flat**2*(c2_flat-N_hits_flat)+2*s_flat*c_flat*m_flat-c2_flat*s_flat**2-\
                N_hits_flat*(c2_flat**2+m_flat**2-c2_flat*N_hits_flat)
        A = -(c2_flat**2+m_flat**2-c2_flat*N_hits_flat)
        B = c_flat*(c2_flat-N_hits_flat)+s_flat*m_flat
        C = c_flat*m_flat-s_flat*c2_flat
        D = -((c2_flat-N_hits_flat)*N_hits_flat+s_flat**2)
        E = c_flat*s_flat-m_flat*N_hits_flat
        F = c2_flat*N_hits_flat-c_flat**2

        I_pixel_flat = (A*I_est_flat+B*Q_est_flat+C*U_est_flat)/Delta
        Q_pixel_flat = (B*I_est_flat+D*Q_est_flat+E*U_est_flat)/Delta
        U_pixel_flat = (C*I_est_flat+E*Q_est_flat+F*U_est_flat)/Delta

        I_pixel = np.reshape(I_pixel_flat, (i,j))
        Q_pixel = np.reshape(Q_pixel_flat, (i,j))
        U_pixel = np.reshape(U_pixel_flat, (i,j))

        return I_pixel, Q_pixel, U_pixel

    def map_multidetectors(self):

        for i in range(self.number):

            mapvalues = map_singledetector(value=self.data[i],sigma=self.weight[i],\
                                           angle=self.polangle[i])
            
            I_map += mapvalues[0]
            Q_map += mapvalues[1]
            U_map += mapvalues[2]

        return I_map, Q_map, U_map





        
