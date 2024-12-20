import os
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, filtfilt
from scipy import interpolate
import sys


# XSens data 

def load_XSens(filename):
    """Load the data from a file.

    Arguments:
        filename {str} -- File path

    Returns
    -------
    Pandas dataframe
        signal
    """

    with open(filename, 'r') as fileID:
        for i, line in enumerate(fileID):
            if line.startswith('PacketCounter'):
                break
        else:
            raise ValueError("\nPacketCounter has not been found in data_lb.")

    signal = pd.read_csv(filename, delimiter="[,\s]", skiprows=i, header=0)
    t = signal["PacketCounter"]
    t_0 = t[0]
    t_fin = t[len(t) - 1]

    time = [i for i in range(int(t_0), int(t_fin) + 1)]
    time_init_0 = [i for i in range(len(time))]
    d = {'PacketCounter': time_init_0}

    colonnes = signal.columns

    for colonne in colonnes[1:]:
        val = signal[colonne]
        f = interpolate.interp1d(t, val)
        y = f(time)
        d[colonne] = y.tolist()

    signal = pd.DataFrame(data=d)

    return signal


def import_XSens(filename, freq, start=0, end=200, order=8, fc=14):
    """Import and pre-process the data from a file.

    Arguments:
        filename {str} -- file path
        freq {int} -- acquisition frequency
        start {int} -- start of the calibration period
        end {int} -- end of the calibration period
        order {int} -- order of the Butterworth low-pass filter
        fc {int} -- cut-off frequency of the Butterworth low-pass filter

    Returns
    -------
    Pandas dataframe
        data
    """
    
    data = load_XSens(filename)
    
    data["FreeAcc_X"] = data["Acc_X"] - np.mean(data["Acc_X"][start:end])
    data["FreeAcc_Y"] = data["Acc_Y"] - np.mean(data["Acc_Y"][start:end])
    data["FreeAcc_Z"] = data["Acc_Z"] - np.mean(data["Acc_Z"][start:end])

    data = filter_sig(data, "Acc", order, fc, fe=freq)
    data = filter_sig(data, "FreeAcc", order, fc, fe=freq)
    data = filter_sig(data, "Gyr", order, fc, fe=freq)

    return data


def filter_sig(data, type_sig, order, fc, fe):
    """Application of Butterworth low-pass filter to a Dataframe

    Arguments:
        data {dataframe} -- pandas dataframe
        type_sig {str} -- "Acc", "Gyr" or "Mag"
        order {int} -- order of the Butterworth low-pass filter
        fc {int} -- cut-off frequency of the Butterworth low-pass filter

    Returns
    -------
    Pandas dataframe
        data
    """
    
    data[type_sig + "_X"] = low_pass_filter(data[type_sig + "_X"], order, fc, fe)
    data[type_sig + "_Y"] = low_pass_filter(data[type_sig + "_Y"], order, fc, fe)
    data[type_sig + "_Z"] = low_pass_filter(data[type_sig + "_Z"], order, fc, fe)

    return data
    

def low_pass_filter(sig, order=8, fc=14, fe=100):
    """Definition of a Butterworth low-pass filter

    Arguments:
        sig {dataframe} -- pandas dataframe
        order {int} -- order of the Butterworth low-pass filter
        fc {int} -- cut-off frequency of the Butterworth low-pass filter
        fe {int} -- acquisition frequency for the data
    Returns
    -------
    ndarray
        filter
    """
    
    f_nyq = fe / 2.  # Hz

    # definition of the Butterworth low-pass filter
    (b, a) = butter(N=order, Wn=(fc / f_nyq), btype='low', analog=False)

    # application
    return filtfilt(b, a, sig)


