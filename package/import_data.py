import pandas as pd
import numpy as np
import json
from scipy.signal import butter, filtfilt
from scipy import interpolate


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

    signal = pd.read_csv(filename, delimiter="[,\s]", skiprows=i, header=0, engine="python")
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


def import_XSens(path, freq, t_start=0.0, t_end=2.0, order=8, fc=14):
    """...
        t_start {float} -- start of the calibration period, in seconds
        t_end {float} -- end of the calibration period, in seconds
    """

    data = load_XSens(path)

    start = int(round(t_start * freq))
    end = int(round(t_end * freq))
    if end > len(data):
        raise ValueError(
            "The recording is shorter than the baseline window ({:.1f} s): a standing period is required at the beginning of the trial.".format(t_end))

    data["FreeAcc_X"] = data["Acc_X"] - np.mean(data["Acc_X"][start:end])
    data["FreeAcc_Y"] = data["Acc_Y"] - np.mean(data["Acc_Y"][start:end])
    data["FreeAcc_Z"] = data["Acc_Z"] - np.mean(data["Acc_Z"][start:end])

    data = filter_sig(data, "Acc", order, fc, freq)
    data = filter_sig(data, "FreeAcc", order, fc,  freq)
    data = filter_sig(data, "Gyr", order, fc, freq)

    return data


def filter_sig(data, type_sig, order, fc, freq):
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
    for axis in ("X", "Y", "Z"):
        data[type_sig + "_" + axis] = low_pass_filter(data[type_sig + "_" + axis], order, fc, freq)

    return data


def low_pass_filter(sig, order, fc, fe):
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
