import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matrixprofile as mp
from scipy.spatial.distance import cdist
from tslearn import metrics
import dtwalign
from scipy import stats

from package import deal_stride


def annotate_stride_estimation(data_1, data_2, pied, r=2, comp=["vide"], freq=100, output=0):
    """Plot the final figure for step detection and save the fig in the output folder as png file. 

    Parameters
    ----------
        steps_lim {dataframe} -- pandas dataframe with the detected gait events
        seg_lim {dataframe} -- pandas dataframe with phases events
        data_rf {dataframe} -- pandas dataframe with data from the right foot sensor
        data_lf {dataframe} -- pandas dataframe with data from the left foot sensor
        freq {int} -- acquisition frequency 
        output {str} -- folder path for output fig
    """
                                 
    gyr_estimation, acc_estimation, start_ref, end_ref = find_stride_estimation(data_1, data_2, pied, freq)
    
    len_estimation = len(gyr_estimation)
    gyr_ref, acc_ref, stride_ref_annotations = find_stride_ref(data, data_2, pied, len_estimation, freq=freq)

    s_y1 = np.array([1 * acc_estimation / (np.max(acc_estimation)), 1 * gyr_estimation / (np.max(abs(gyr_estimation)))])
    s_y1 = s_y1.transpose()

    s_y2 = np.array([1 * acc_ref / (np.max(acc_ref)), 1 * gyr_ref / (np.max(abs(gyr_ref)))])
    s_y2 = s_y2.transpose()

    path, sim = metrics.dtw_path(s_y1, s_y2, global_constraint="itakura", itakura_max_slope=r)

    patho_stride_annotations = deal_stride.annotate(path, stride_ref_annotations)

    comp = plot_annotate_stride_estimation(gyr_estimation, acc_estimation, patho_stride_annotations, s_y1, s_y2, path,
                                           pied, freq=freq, comp=comp, start=start_ref, output=output)

    return gyr_estimation, acc_estimation, patho_stride_annotations, comp
                             

def plot_annotate_stride_estimation(gyr_estimation, acc_estimation, patho_stride_annotations, s_y1, s_y2, path, pied,
                                   freq=100, comp=None, start=0, output=0):

    sz_1 = s_y1.shape[0]
    sz_2 = s_y2.shape[0]

    plt.close('all')

    fig = plt.figure(1, figsize=(8, 8))

    # axes definition
    left, bottom = 0.02, 0.1
    w_ts = h_ts = 0.2
    left_h = left + w_ts + 0.02
    width = 0.65 * sz_2 / max(sz_2, sz_1)
    height = 0.65 * sz_1 / max(sz_2, sz_1)
    bottom_h = bottom + height + 0.02

    rect_s_y = [left, bottom, w_ts, height]
    rect_gram = [left_h, bottom, width, height]
    rect_stride = [2 * left + left_h + 0.65, bottom, 0.65, height]
    rect_s_x = [left_h, bottom_h, width, h_ts]

    ax_gram = plt.axes(rect_gram)
    ax_stride = plt.axes(rect_stride)
    ax_s_x = plt.axes(rect_s_x, sharex=ax_gram)
    ax_s_y = plt.axes(rect_s_y)

    mat = cdist(s_y1, s_y2)

    ax_gram.imshow(mat, origin='lower')
    ax_gram.axis("off")
    ax_gram.autoscale(False)
    ax_gram.plot([j for (i, j) in path], [i for (i, j) in path], "w-",
                 linewidth=3.)

    ax_s_x.plot(np.arange(sz_2), 0.5 + s_y2[:, 0], "b-", linewidth=3.)
    ax_s_x.plot(np.arange(sz_2), - 0.5 + s_y2[:, 1], "y-", linewidth=3.)
    ax_s_x.axis("off")
    ax_s_x.set_xlim((0, sz_2 - 1))

    ax_s_y.plot(- 0.5 - s_y1[:, 0], np.arange(sz_1), "b-", linewidth=3.)
    ax_s_y.plot(0.5 - s_y1[:, 1], np.arange(sz_1), "y-", linewidth=3.)
    ax_s_y.axis("off")
    ax_s_y.set_ylim((0, sz_1 - 1))

    ax_stride.plot(acc_estimation / (np.max(acc_estimation)), "b-", linewidth=3.)
    ax_stride.plot(- 1 + gyr_estimation / (np.max(abs(gyr_estimation))), "y-", linewidth=3.)
    mi, ma = min(- 1 + gyr_estimation / (np.max(abs(gyr_estimation)))), max(acc_estimation / (np.max(acc_estimation)))

    ax_stride.vlines(patho_stride_annotations["HS"], mi, ma, 'black', label="Heel Strike")
    ax_stride.vlines(patho_stride_annotations["FF"], mi, ma, 'violet', label="Foot Flat")
    ax_stride.vlines(patho_stride_annotations["HO"], mi, ma, 'green', label="Heel Off")
    ax_stride.vlines(patho_stride_annotations["TO"], mi, ma, 'red', label="Toe Off")
    ax_stride.legend()
    ax_stride.grid()

    # figure save
    if pied == 1:
        titre = "steps_right.svg"
    if pied == 0:
        titre = "steps_left.svg"
    os.chdir(output)
    plt.savefig(titre, bbox_inches="tight")

    return comp


def find_stride_ref(data, data_autre, pied, len_estimation, freq=100):
    gyr_ref_decal, acc_ref_decal, stride_ref_decal_annotations = deal_stride.stride_sain_decal(int(len_estimation), freq)
    gyr_estimation, acc_estimation, p, q = find_stride_estimation(data, data_autre, pied, freq=freq)

    cout = []
    for j in range(0, len(gyr_estimation)):
        u = gyr_estimation / (np.max(abs(gyr_estimation)))
        v = gyr_ref_decal[str(j)] / (np.max(abs(gyr_ref_decal[str(j)])))
        w = acc_estimation / (np.max(abs(acc_estimation)))
        h = acc_ref_decal[str(j)] / (np.max(abs(acc_ref_decal[str(j)])))
        cout.append(stats.pearsonr(u, v[0:len(u)])[0] + stats.pearsonr(w, h[0:len(w)])[0])

    cout = np.array(cout)
    pic_cout = indexes(cout, thres=0.5 * np.max(cout), min_dist=1, thres_abs=True)
    if cout[0] > 0.5 * np.max(cout):
        pic_cout = np.append(pic_cout, 0)

    pic_cout = pic_cout[np.argsort(-cout[pic_cout])]
    pic_cout = pic_cout[:3]
    pic_cout.sort()

    if len(pic_cout) == 3:
        p1_3 = pic_cout[2] - pic_cout[0]
        p2_1 = pic_cout[0] + len(gyr_estimation) - pic_cout[1]
        p3_2 = pic_cout[1] + len(gyr_estimation) - pic_cout[2]
        if (p1_3 <= p2_1) & (p1_3 <= p3_2):
            decal_estim = pic_cout[1]
        else:
            if (p3_2 <= p2_1):
                decal_estim = pic_cout[0]
            else:
                decal_estim = pic_cout[2]
    else:
        if len(pic_cout) == 2:
            max_cout = max(cout[pic_cout[1]], cout[pic_cout[0]])
            cout_norm_0 = max(0, cout[pic_cout[0]] - 0.75 * max_cout)
            cout_norm_1 = max(0, cout[pic_cout[1]] - 0.75 * max_cout)
            if abs(pic_cout[0] - pic_cout[1]) < len(gyr_estimation) // 2:
                cout1 = cout_norm_1/(cout_norm_0 + cout_norm_1)
                decal_estim = int(min(pic_cout[0], pic_cout[1]) + abs(cout1 * (pic_cout[0] - pic_cout[1])))
            else:
                cout0 = cout_norm_0 / (cout_norm_1 + cout_norm_0)
                ecart = min(pic_cout[0], pic_cout[1]) + len(gyr_estimation) - max(pic_cout[0], pic_cout[1])
                decal_estim = int(max(pic_cout[0], pic_cout[1]) + ecart * cout0)
        else:
            if len(pic_cout) == 0:
                pic_cout = [0]
            decal_estim = int(np.median(pic_cout))

    decal_estim = decal_estim % len(gyr_estimation)

    return gyr_ref_decal[str(decal_estim)], acc_ref_decal[str(decal_estim)], stride_ref_decal_annotations[str(decal_estim)]


def find_stride_estimation(data_1, data_2, pied, freq=100):
                             
    x = data["Gyr_Y"]
    z = deal_stride.calculate_jerk_tot(data)
    t = data["PacketCounter"]

    window = int(len_stride_estimation(data_1, data_2, roll=1, freq=freq))

    # matrix profile
    mp_profile = mp.compute(x.to_numpy(), windows=window)

    av = []
    x_norm = x / np.max(abs(x))
    z_norm = z / np.max(z)
    for i in range(len(x) - window + 1):
        av.append(np.sum(abs(x_norm[i + window // 3: i + 2 * window // 3])) # ** 2
                  + np.sum(abs(z_norm[i + window // 3: i + 2 * window // 3])))

    av = av / np.max(av)

    mp_profile = mp.transform.apply_av(mp_profile, "custom", av)
    mp_profile = mp.discover.motifs(mp_profile, k=1, use_cmp=True)

    start_ref = mp_profile['motifs'][0]["motifs"][0]
    end_ref = start_ref + window

    return x[start_ref:end_ref].to_numpy(), z[start_ref:end_ref], start_ref, end_ref


def len_stride_estimation(data_1, data_2, roll=1, freq=100):
                            
    len_stride_data_1 = len_stride_one_side(data_1, roll=roll, freq=freq)
    len_stride_data_2 = len_stride_one_side(data_2, roll=roll, freq=freq)
    
    if len_stride_data_1 / len_stride_data_2 >= 1.5:
        return len_stride_data_2
    else:
        return len_stride_data_1


def len_stride_one_side(data, roll=1, freq=100):
    """Import and pre-process the data from a file.

    Arguments:
        data {pandas Dataframe} -- da
        start {int} -- start of the calibration period
        end {int} -- end of the calibration period
        order {int} -- order of the Butterworth low-pass filter
        fc {int} -- cut-off frequency of the Butterworth low-pass filter

    Returns
    -------
    Pandas dataframe
        data
    """
    
    x_11 = data["FreeAcc_X"]
    test_11 = x_11.to_numpy()
    x_12 = data["FreeAcc_Y"]
    test_12 = x_12.to_numpy()
    x_13 = data["FreeAcc_Z"]
    test_13 = x_13.to_numpy()
    x_2 = data["Gyr_Y"]
    test_2 = x_2.to_numpy()
    acf = (autocorr(test_11) / 3 + autocorr(test_12) / 3 + autocorr(test_13)) / 3 / 2 + autocorr(test_2) / 2

    y = pd.DataFrame(acf)
    y_mean = y.rolling(roll, center=True, win_type='cosine').mean()
    y_mean = y_mean.fillna(0)
    y_mean_np = y_mean.to_numpy().transpose()[0]

    index_pic = autocorr_indexes(y_mean_np[:len(y_mean_np) // 4], freq=freq)

    if len(index_pic) > 0:
      return index_pic[0]
    else:
        return 0


def autocorr_indexes(y, thres=0.7, min_dist=80, freq=100):
    """Find autocorrelation local maxima indexes.

    Parameters
    ----------
    y {ndarray} -- containing non-biased autocorrelation.
    thres {float} -- normalized threshold between [0., 1.]. Only the peaks with amplitude higher than the threshold will be detected.
    min_dist {int} -- minimum distance between two peaks. 
        
    Returns
    -------
    ndarray containing the numeric indexes of the non-biased autocorrelation that were detected.
    """
    
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    i = round(min_dist*freq/100)
    thres = thres * (np.max(y[i:]) - np.max(np.min(y[i:i + np.argmax(y[i:])]), 0)) + np.max(
        np.min(y[i:i + np.argmax(y[i:])]), 0)

    index_list = indexes(y, thres=thres, min_dist=min_dist, thres_abs=True)

    return index_list[(index_list > i)]
  

def autocorr(f):
    """Autocorrelation non-biased indicator.

    Parameters
    ----------
    f -- ndarray 
        1D data to compute autocorrelation.
        
    Returns
    -------
    acf -- ndarray
        Array containing non-biased autocorrelation.
    """
    N = len(f)
    fvi = np.fft.fft(f, n=2 * N)
    acf = np.real(np.fft.ifft(fvi * np.conjugate(fvi))[:N])
    d = N - np.arange(N)
    acf = acf / d  #  non biased indicator
    acf = acf / acf[0]
  
    return acf


def indexes(y, thres=0.3, min_dist=1, thres_abs=False):
    """Peak detection routine.

    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.

    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    thres_abs: boolean
        If True, the thres value will be interpreted as an absolute value, instead of
        a normalized threshold.

    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected.
        When using with Pandas DataFrames, iloc should be used to access the values at the returned positions.
    """
  
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    if not thres_abs:
        thres = thres * (np.max(y) - np.min(y)) + np.min(y)

    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    zeros, = np.where(dy == 0)

    # check if the signal is totally flat
    if len(zeros) == len(y) - 1:
        print("Signal seems to be flat !")
        return np.array([])

    if len(zeros):
        # compute first order difference of zero indexes
        zeros_diff = np.diff(zeros)
        # check when zeros are not chained together
        zeros_diff_not_one, = np.add(np.where(zeros_diff != 1), 1)
        # make an array of the chained zero indexes
        zero_plateaus = np.split(zeros, zeros_diff_not_one)

        # fix if leftmost value in dy is zero
        if zero_plateaus[0][0] == 0:
            dy[zero_plateaus[0]] = dy[zero_plateaus[0][-1] + 1]
            zero_plateaus.pop(0)

        # fix if rightmost value of dy is zero
        if len(zero_plateaus) and zero_plateaus[-1][-1] == len(dy) - 1:
            dy[zero_plateaus[-1]] = dy[zero_plateaus[-1][0] - 1]
            zero_plateaus.pop(-1)

        # for each chain of zero indexes
        for plateau in zero_plateaus:
            median = np.median(plateau)
            # set leftmost values to leftmost non zero values
            dy[plateau[plateau < median]] = dy[plateau[0] - 1]
            # set rightmost and middle values to rightmost non zero values
            dy[plateau[plateau >= median]] = dy[plateau[-1] + 1]

    # find the peaks by using the first order difference
    peaks = np.where(
        (np.hstack([dy, 0.0]) < 0.0)
        & (np.hstack([0.0, dy]) > 0.0)
        & (np.greater(y, thres))
    )[0]

    # handle multiple peaks, respecting the minimum distance
    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks
