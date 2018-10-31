import map as mp
import numpy as np
import sys 
import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

map_file = sys.argv[1]

#Read Map Parameters 

d = mp.dataload(map_file)

param = d.map_param()

#Files
dirfile = param[0] 
detlist = param[1] 
coor1type = param[2] 
coor2type = param[3] 
pointingfile = param[4]
dettable = param[5]
detpath = param[6]
coorpath = param[7]
#Map Parameters
ctype = param[8] 
crpix = param[9]
crdelt = param[10]
crval = param[11] 
conv = param[12]
stdev = param[13]
#Experiment Parameters
detfreq = param[14]
acsfreq = param[15]
det_dir_conv = param[16] 
coor1_dir_conv = param[17] 
coor2_dir_conv = param[18] 
frames = param[19]
det_samp_frame = param[20] 
acs_samp_frame = param[21]
det_file_type = param[22]
coor1_file_type = param[23]
coor2_file_type = param[24]

if dirfile.lower() == 'na':

    if np.size(detlist) > 1:
        det_value_0 = d.loaddata(detpath, detlist[0], det_file_type)

        detTOD = np.zeros((len(detlist), len(det_value_0)))
        detTOD[0,:] = det_value_0

        for i in range(1,len(detlist)):
            detTOD[i,:] = d.loaddata(detpath, detlist[i], det_file_type)
    else:
        detTOD = d.loaddata(detpath, detlist, det_file_type)

    coord1 = d.loaddata(coorpath, coor1type, coor1_file_type)
    coord2 = d.loaddata(coorpath, coor2type, coor2_file_type)
else:

    detTOD, coord1, coord2 = d.loadfulldata(dirfile, detlist, det_file_type, coord_type1, \
                                            coor1_file_type, coor2_file_type)

detTOD = d.convert_dirfile(detTOD, det_dir_conv[0], det_dir_conv[1])
coord1 = d.convert_dirfile(coord1, coor1_dir_conv[0], coor1_dir_conv[1])
coord2 = d.convert_dirfile(coord2, coor2_dir_conv[0], coor2_dir_conv[1])

dettime, detTOD = d.frame_zoom(detTOD, det_samp_frame, detfreq, frames)

coord1time, coord1 = d.frame_zoom(coord1, acs_samp_frame, acsfreq, frames)
coord2time, coord2 = d.frame_zoom(coord2, acs_samp_frame, acsfreq, frames)

coord1, coord2 = d.coord_int(coord1, coord2, coord1time, dettime)

if coor1type.lower() == 'ra':
    coord1 = coord1*15.

wcsworld = mp.wcs_world(ctype, crpix, crdelt, crval)

w, proj = wcsworld.world(np.transpose(np.array([coord1,coord2])))

filterdat = mp.filterdata(detTOD, 0.2, detfreq) #0.1 is the frequency cutoff for the high pass filter 


cleanedata = filterdat.ifft_filter()


mapmaker = mp.mapmaking(detTOD, 1., 1., np.size(detlist), np.round(w).astype(int))

if np.size(detlist) > 1:
    finalI = mapmaker.map_multidetector_Ionly()

else:
    finalI = mapmaker.map_singledetector_Ionly()

if conv.lower() != 'na':

    finalI_conv = mapmaker.convolution(stdev, finalI)

fig = plt.figure()
fig.add_subplot(111, projection=proj)
plt.imshow(finalI[0], origin='lower', cmap=plt.cm.viridis)
plt.colorbar()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.show()