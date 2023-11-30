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


def annotate_ref_stride(data_1, data_2, foot, r=2, freq=100, output=0):
    """Find, annote and plot a signal subset of the foot of interest to be considered as the reference stride. 
    Annotation of the model stride with the gait events (TO, HS, FF, HO). 
    
    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        foot {int} -- 0 for left, 1 for right
        r {float} -- size of the Itakura parallelogram to restrict the DTW
        freq {int} -- acquisition frequency (Hz)
        output {str} -- folder path for output fig

    Returns
    -------
    ndarray, ndarray, ndarray
       gyration time series of the reference stride {ndarray}
       jerk time series of the reference stride {ndarray}
       annotation of gait events of the reference stride in the trial {ndarray}
    """
                                 
    gyr_ref, acc_ref, start_ref, end_ref = find_ref_stride(data_1, data_2, foot, freq)
    
    len_ref = len(gyr_ref)
    gyr_model, acc_model, stride_model_annotations = find_model_stride(data_1, data_2, foot, len_ref, freq=freq)

    s_y1 = np.array([1 * acc_ref / (np.max(acc_ref)), 1 * gyr_ref / (np.max(abs(gyr_ref)))])
    s_y1 = s_y1.transpose()

    s_y2 = np.array([1 * acc_model / (np.max(acc_model)), 1 * gyr_ref / (np.max(abs(gyr_model)))])
    s_y2 = s_y2.transpose()

    path, sim = metrics.dtw_path(s_y1, s_y2, global_constraint="itakura", itakura_max_slope=r)

    ref_stride_annotations = deal_stride.annotate(path, stride_model_annotations)

    plot_annotate_ref_stride(gyr_ref, acc_ref, ref_stride_annotations, s_y1, s_y2, path,
                                           foot, freq=freq, start=start_ref, output=output)

    return gyr_ref, acc_ref, ref_stride_annotations
                             

def plot_annotate_ref_stride(gyr_ref, acc_ref, ref_stride_annotations, s_y1, s_y2, path, foot,
                                   freq=100, start=0, output=0):

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

    ax_stride.plot(acc_ref / (np.max(acc_ref)), "b-", linewidth=3.)
    ax_stride.plot(- 1 + gyr_ref / (np.max(abs(gyr_ref))), "y-", linewidth=3.)
    mi, ma = min(- 1 + gyr_ref / (np.max(abs(gyr_ref)))), max(acc_ref / (np.max(acc_ref)))

    ax_stride.vlines(ref_stride_annotations["HS"], mi, ma, 'black', label="Heel Strike")
    ax_stride.vlines(ref_stride_annotations["FF"], mi, ma, 'violet', label="Foot Flat")
    ax_stride.vlines(ref_stride_annotations["HO"], mi, ma, 'green', label="Heel Off")
    ax_stride.vlines(ref_stride_annotations["TO"], mi, ma, 'red', label="Toe Off")
    ax_stride.legend()
    ax_stride.grid()

    # figure save
    if foot == 1:
        titre = "stride_right.svg"
    if foot == 0:
        titre = "stride_left.svg"
    os.chdir(output)
    plt.savefig(titre, bbox_inches="tight")


def find_model_stride(data_1, data_2, foot, len_ref, freq=100):
    """Find a signal subset of the foot of interest to be considered as the model stride. 
    Annotation of the model stride with the gait events (TO, HS, FF, HO). 
    
    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        foot {int} -- 0 for left, 1 for right
        len_ref {int} -- size of the attended model stride
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    ndarray, ndarray, ndarray
       gyration time series of the found model stride {ndarray}
       jerk time series of the found model stride {ndarray}
       annotation of gait events of the found model stride in the trial {ndarray}
    """
    
    gyr_model_decal, acc_model_decal, stride_model_decal_annotations = deal_stride.stride_sain_decal(int(len_ref), freq)
    gyr_ref, acc_ref, p, q = find_ref_stride(data_1, data_2, foot, freq=freq)

    cout = []
    for j in range(0, len(gyr_ref)):
        u = gyr_ref / (np.max(abs(gyr_ref)))
        v = gyr_model_decal[str(j)] / (np.max(abs(gyr_model_decal[str(j)])))
        w = acc_ref / (np.max(abs(acc_ref)))
        h = acc_model_decal[str(j)] / (np.max(abs(acc_model_decal[str(j)])))
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
        p2_1 = pic_cout[0] + len(gyr_ref) - pic_cout[1]
        p3_2 = pic_cout[1] + len(gyr_ref) - pic_cout[2]
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
            if abs(pic_cout[0] - pic_cout[1]) < len(gyr_ref) // 2:
                cout1 = cout_norm_1/(cout_norm_0 + cout_norm_1)
                decal_estim = int(min(pic_cout[0], pic_cout[1]) + abs(cout1 * (pic_cout[0] - pic_cout[1])))
            else:
                cout0 = cout_norm_0 / (cout_norm_1 + cout_norm_0)
                ecart = min(pic_cout[0], pic_cout[1]) + len(gyr_ref) - max(pic_cout[0], pic_cout[1])
                decal_estim = int(max(pic_cout[0], pic_cout[1]) + ecart * cout0)
        else:
            if len(pic_cout) == 0:
                pic_cout = [0]
            decal_estim = int(np.median(pic_cout))

    decal_estim = decal_estim % len(gyr_ref)

    return gyr_model_decal[str(decal_estim)], acc_model_decal[str(decal_estim)], stride_model_decal_annotations[str(decal_estim)]


def find_ref_stride(data_1, data_2, foot, freq=100):
    """Find a signal subset of the foot of interest to be considered as the reference stride. 
    Selection using a matrix profile technique with an annotation vector. 
    
    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        foot {int} -- 0 for left, 1 for right
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    ndarray, ndarray, int, int
       gyration time series of the found reference stride {ndarray}
       jerk time series of the found reference stride {ndarray}
       beginning sample number of the found reference stride in the trial {int}
       ending sample number of the found reference stride in the trial {int}
    """

    # signals of interest: time, gyration in the axial plane, jerk norm
    t = data_1["PacketCounter"]
    x = data_1["Gyr_Y"]
    z = deal_stride.calculate_jerk_tot(data_1)

    # search window size : mean stride time estimation
    window = int(len_stride_estimation(data_1, data_2, roll=1, freq=freq))

    # matrix profile
    mp_profile = mp.compute(x.to_numpy(), windows=window)

    # annotation vector for matrix profile in order to promote swing phase in the center of the window
    av = []
    x_norm = x / np.max(abs(x))
    z_norm = z / np.max(z)
    for i in range(len(x) - window + 1):
        av.append(np.sum(abs(x_norm[i + window // 3: i + 2 * window // 3])) # ** 2
                  + np.sum(abs(z_norm[i + window // 3: i + 2 * window // 3])))
    av = av / np.max(av)
    mp_profile = mp.transform.apply_av(mp_profile, "custom", av)
    mp_profile = mp.discover.motifs(mp_profile, k=1, use_cmp=True)

    # extraction of the beginning and end of the sub-series which will be the reference stride
    start_ref = mp_profile['motifs'][0]["motifs"][0]
    end_ref = start_ref + window

    return x[start_ref:end_ref].to_numpy(), z[start_ref:end_ref], start_ref, end_ref


def len_stride_estimation(data_1, data_2, roll=1, freq=100):
    """Estimate the mean stride time from both feet data computed from autocorrelations. The first foot is the interest foot.  
    The second foot acts as a safety net, in case autocorrelation is limiting on the interest foot. 
    
    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        roll {int} -- size of the window center rolling. Default is 1, meaning no window rolling
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    int
        mean stride time estimation
    """
                            
    len_stride_data_1 = len_stride_one_side(data_1, roll=roll, freq=freq)
    len_stride_data_2 = len_stride_one_side(data_2, roll=roll, freq=freq)

    # if the two estimates are too far apart, it's likely that one of the peak detections is faulty (the two estimates should be equal).
    # if the foot of interest is defective, we take the lowest. 
    if len_stride_data_1 / len_stride_data_2 >= 1.5:
        return len_stride_data_2
    else:
        return len_stride_data_1


def len_stride_one_side(data, roll=1, freq=100):
    """Estimate the mean stride time from one foot data computed from autocorrelations.
    
    Arguments:
        data {pandas Dataframe} -- dataframe with data from one of the foot sensor
        roll {int} -- size of the window center rolling. Default is 1, meaning no window rolling
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    int
        index of the first peak of the time series computed from autocorrelations of FreeAcc_X, FreeAcc_Y, FreeAcc_Z, and Gyr_Y
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
