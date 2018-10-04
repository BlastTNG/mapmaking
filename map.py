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
        b, a = highpass(order)
        self.data = sgn.lfilter(b, a, self.data)
        return self.data

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

    def rotmatrix(self, yaw_mat = self.yaw, roll_mat = self.roll, pitch_mat=self.pitch):
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
        
class mapmaking():

    def __init__(self, data, weight, polangle, number):

        self.data = data
        self.weight = weight
        self.polangle = polangle
        self.number = number

    def map_param(self, value=self.data, sigma=self.weight, angle=self.polangle):

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

        return I_est_flat, Q_est_flat, U_est_flat, Delta, A, B, C, D, E, F

    def map_singledetector(self, value=self.data, sigma=self.weight, angle=self.polangle):

        I_est_flat, Q_est_flat, U_est_flat, Delta, A, B, C, D, E, F = map_param(value=self.data, \
                                                                                sigma=self.weight,\
                                                                                angle=self.polangle)

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


        






        
