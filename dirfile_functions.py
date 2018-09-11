import numpy as np
import pygetdata as gd
import matplotlib.pyplot as plt
import pykst as kst

fs = 488.28125 # Hz

def loadIQAllChan(dirfile, nsec):
    nframes = np.int(nsec * fs)
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    ifiles = [v for v in vectors if v[0] == "i" and v != "INDEX"]
    qfiles = [v for v in vectors if v[0] == "q" and v != "INDEX"]
    ifiles = sorted(ifiles, key=lambda x: x.split('_')[1])
    qfiles = sorted(qfiles, key=lambda x: x.split('_')[1])
    ivals, qvals = [], []
    for i in range(len(ifiles)):
        ivals.append(d.getdata(ifiles[i], gd.FLOAT32, num_frames = nframes))
        qvals.append(d.getdata(qfiles[i], gd.FLOAT32, num_frames = nframes))
    d.close()
    ivals = np.asarray(ivals)
    qvals = np.asarray(qvals)
    return ivals, qvals

def loadIQsingleChan(dirfile, chan):
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    ifiles = [v for v in vectors if v[0] == "i" and v != "INDEX"]
    qfiles = [v for v in vectors if v[0] == "q" and v != "INDEX"]
    ifiles = sorted(ifiles, key=lambda x: x.split('_')[1])
    qfiles = sorted(qfiles, key=lambda x: x.split('_')[1])
    ivals = d.getdata(ifiles[chan], gd.FLOAT32, num_frames = d.nframes)
    qvals = d.getdata(qfiles[chan], gd.FLOAT32, num_frames = d.nframes)
    d.close()
    return np.asarray(ivals, dtype = 'float'), np.asarray(qvals, dtype = 'float')

def avgAllChan(dirfile):
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    ifiles = [v for v in vectors if v[0] == "i" and v != "INDEX"]
    qfiles = [v for v in vectors if v[0] == "q" and v != "INDEX"]
    ifiles = sorted(ifiles, key=lambda x: x.split('_')[1])
    qfiles = sorted(qfiles, key=lambda x: x.split('_')[1])
    ivals = np.zeros(d.nframes)
    qvals = np.zeros(d.nframes)
    count = 0
    for chan in range(len(ifiles)):
        print count
        ivals += d.getdata(ifiles[chan], gd.FLOAT32, num_frames = d.nframes)
        qvals += d.getdata(qfiles[chan], gd.FLOAT32, num_frames = d.nframes)
        count += 1
    d.close()
    return np.asarray(ivals, dtype = 'float'), np.asarray(qvals, dtype = 'float')

def loadChanPhase(dirfile, chan):
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    ifiles = [v for v in vectors if v[0] == "i" and v != "INDEX"]
    qfiles = [v for v in vectors if v[0] == "q" and v != "INDEX"]
    ifiles = sorted(ifiles, key=lambda x: x.split('_')[1])
    qfiles = sorted(qfiles, key=lambda x: x.split('_')[1])
    ivals = d.getdata(ifiles[chan], gd.FLOAT32, num_frames = d.nframes)
    qvals = d.getdata(qfiles[chan], gd.FLOAT32, num_frames = d.nframes)
    d.close()
    phi_rot = np.arctan2(np.mean(qvals), np.mean(ivals))
    s21 = ivals + 1j*qvals
    s21_rot = s21 * np.exp(-1j * phi_rot)
    phase = np.angle(s21_rot)
    return np.asarray(phase, dtype = 'float')

def loadChanPhaseRange(dirfile, chan, first_sample, num_samples):
    d = gd.dirfile(dirfile, gd.RDWR|gd.UNENCODED)
    vectors = d.field_list()
    ifiles = [v for v in vectors if v[0] == "i" and v != "INDEX"]
    qfiles = [v for v in vectors if v[0] == "q" and v != "INDEX"]
    ifiles = sorted(ifiles, key=lambda x: x.split('_')[1])
    qfiles = sorted(qfiles, key=lambda x: x.split('_')[1])
    ivals = d.getdata(ifiles[chan], gd.FLOAT32,
             first_sample = first_sample, num_samples = num_samples)
    qvals = d.getdata(qfiles[chan], gd.FLOAT32,
             first_sample = first_sample, num_samples = num_samples)
    d.close()
    phi_rot = np.arctan2(np.mean(qvals), np.mean(ivals))
    s21 = ivals + 1j*qvals
    s21_rot = s21 * np.exp(-1j * phi_rot)
    phase = np.angle(s21_rot)
    return np.asarray(phase, dtype = 'float')

def plotVectorKst(x, y, xname, yname):
    # start a kst session with the arbitrary name "NumpyVector"
    client=kst.Client("Vector")
    # copy the numpy arrays into kst
    V1 = client.new_editable_vector(x, name=xname) # the name is for the label
    V2 = client.new_editable_vector(y, name=yname) # the name is for the label
    # inside kst, create a curve, a plot, and add the curve to the plot
    c1 = client.new_curve(V1, V2)
    p1 = client.new_plot()
    p1.add(c1)
    return

