import os
import numpy as np
import matplotlib.pyplot as plt
import stumpy
from scipy.spatial.distance import cdist
from scipy.signal import find_peaks
from scipy import stats
from tslearn import metrics

from package import deal_stride


def annotate_ref_stride(data_1, data_2, foot, freq, r=2, output=0):
    """Find, annotate and plot a signal subset of the foot of interest to be considered as the reference stride.
    Annotation of the reference stride with the gait events (TO, HS, FF, HO).

    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        foot {int} -- 0 for left, 1 for right
        freq {int} -- acquisition frequency (Hz)
        r {float} -- maximum slope of the Itakura parallelogram constraining the DTW
        output {str} -- folder path for output fig

    Returns
    -------
    ndarray, ndarray, dict
       gyration time series of the reference stride
       jerk time series of the reference stride
       gait event indexes within the reference stride
    """

    # reference stride: extracted from the recording itself (step 2)
    gyr_ref, jerk_ref, start_ref, end_ref = find_ref_stride(data_1, data_2, foot, freq)

    # model stride: always the same healthy subject stride, stretched and shifted (step 3.1-3.2)
    gyr_model, jerk_model, stride_model_annotations = find_model_stride(gyr_ref, jerk_ref, freq)

    # data formatting: each channel is normalized by its own maximum so that both
    # contribute on a comparable scale to the local cost of the DTW
    s_y1 = np.array([jerk_ref / np.max(jerk_ref), gyr_ref / np.max(abs(gyr_ref))])
    s_y1 = s_y1.transpose()
    s_y2 = np.array([jerk_model / np.max(jerk_model), gyr_model / np.max(abs(gyr_model))])
    s_y2 = s_y2.transpose()

    # multidimensional DTW with dependence (DTW_D): the two channels form a single 2-D
    # sequence, the local cost being the Euclidean norm in R^2 (step 3.3)
    path, sim = metrics.dtw_path(s_y1, s_y2, global_constraint="itakura", itakura_max_slope=r)

    # annotate the reference stride with the path and the model stride
    ref_stride_annotations = deal_stride.annotate(path, stride_model_annotations)

    # plot the result
    plot_annotate_ref_stride(gyr_ref, jerk_ref, ref_stride_annotations, s_y1, s_y2, path,
                             foot, freq, start=start_ref, output=output)

    return gyr_ref, jerk_ref, ref_stride_annotations


def plot_annotate_ref_stride(gyr_ref, jerk_ref, ref_stride_annotations, s_y1, s_y2, path, foot, freq, start=0, output=0):
    """Plot the reference stride and its construction.
    Annotation of the reference stride with the gait events (TO, HS, FF, HO).

    Arguments:
        gyr_ref {ndarray} -- array with the gyration from the reference stride
        jerk_ref {ndarray} -- array with the jerk from the reference stride
        ref_stride_annotations {dict} -- index of the 4 gait events of the reference stride: HS, FF, HO, TO
        s_y1 {ndarray} -- transposed 2-D reference stride sequence
        s_y2 {ndarray} -- transposed 2-D model stride sequence
        path {list} -- list of pairs giving the DTW correspondence between both sequences
        foot {int} -- 0 for left, 1 for right
        freq {int} -- acquisition frequency (Hz)
        output {str} -- folder path for output fig

    Returns
    -------
    fig
    """

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

    # matrix computation
    mat = cdist(s_y1, s_y2)

    ax_gram.imshow(mat, origin='lower')
    ax_gram.axis("off")
    ax_gram.autoscale(False)
    ax_gram.plot([j for (i, j) in path], [i for (i, j) in path], "w-", linewidth=3.)

    ax_s_x.plot(np.arange(sz_2), 0.5 + s_y2[:, 0], "b-", linewidth=3.)
    ax_s_x.plot(np.arange(sz_2), - 0.5 + s_y2[:, 1], "y-", linewidth=3.)
    ax_s_x.axis("off")
    ax_s_x.set_xlim((0, sz_2 - 1))

    ax_s_y.plot(- 0.5 - s_y1[:, 0], np.arange(sz_1), "b-", linewidth=3.)
    ax_s_y.plot(0.5 - s_y1[:, 1], np.arange(sz_1), "y-", linewidth=3.)
    ax_s_y.axis("off")
    ax_s_y.set_ylim((0, sz_1 - 1))

    ax_stride.plot(jerk_ref / (np.max(jerk_ref)), "b-", linewidth=3.)
    ax_stride.plot(- 1 + gyr_ref / (np.max(abs(gyr_ref))), "y-", linewidth=3.)
    mi, ma = min(- 1 + gyr_ref / (np.max(abs(gyr_ref)))), max(jerk_ref / (np.max(jerk_ref)))

    # annotations computation
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


def find_model_stride(gyr_ref, jerk_ref, freq):
    """Stretch the model stride so that it matches the duration of the reference stride,
    and circularly shift it so that its swing phase is aligned with that of the reference.

    Arguments:
        gyr_ref {ndarray} -- gyration time series of the reference stride
        jerk_ref {ndarray} -- jerk time series of the reference stride
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    ndarray, ndarray, dict
       gyration time series of the aligned model stride
       jerk time series of the aligned model stride
       gait event indexes of the aligned model stride
    """

    # model stride, stretched to the reference stride duration
    gyr_model_stretch, jerk_model_stretch, stride_model_stretch_annotations = \
        deal_stride.model_stride_offset(int(len(gyr_ref)), freq)

    # goal: find the circular shift that maximizes the joint Pearson correlation
    # between the model stride and the reference stride
    correlation = []
    for j in range(0, len(gyr_ref)):
        u = gyr_ref / np.max(abs(gyr_ref))
        v = gyr_model_stretch[str(j)] / np.max(abs(gyr_model_stretch[str(j)]))
        w = jerk_ref / np.max(abs(jerk_ref))
        h = jerk_model_stretch[str(j)] / np.max(abs(jerk_model_stretch[str(j)]))
        correlation.append(stats.pearsonr(u, v[0:len(u)])[0] + stats.pearsonr(w, h[0:len(w)])[0])

    correlation = np.array(correlation)
    stretch_estim = int(np.argmax(correlation)) % len(gyr_ref)

    return (gyr_model_stretch[str(stretch_estim)],
            jerk_model_stretch[str(stretch_estim)],
            stride_model_stretch_annotations[str(stretch_estim)])


def find_ref_stride(data_1, data_2, foot, freq):
    """Find a signal subset of the foot of interest to be considered as the reference stride.
    Selection using a matrix profile technique corrected by an annotation vector.

    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        foot {int} -- 0 for left, 1 for right
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    ndarray, ndarray, int, int
       gyration time series of the found reference stride
       jerk time series of the found reference stride
       beginning sample number of the found reference stride in the trial
       ending sample number of the found reference stride in the trial
    """

    # signals of interest: gyration in the sagittal plane, jerk norm
    x = data_1["Gyr_Y"].to_numpy(dtype=float)
    z = np.asarray(deal_stride.calculate_jerk_tot(data_1, freq), dtype=float)

    # search window size: mean stride duration estimation (step 1)
    window, autocorr_value = len_stride_estimation(data_1, data_2, freq)
    window = int(window)

    if window < 3 or window > len(x):
        raise ValueError(
            "The estimated stride duration ({} samples) is not compatible with the length of the "
            "recording ({} samples).".format(window, len(x)))

    # matrix profile and annotation vector (step 2.1 and 2.2)
    mp_values = matrix_profile(x, window)
    av = annotation_vector(x, z, window)
    cmp_values = corrected_matrix_profile(mp_values, av)

    # extraction of the beginning and end of the sub-series which will be the reference stride
    start_ref = int(np.argmin(cmp_values))
    end_ref = start_ref + window

    return x[start_ref:end_ref], z[start_ref:end_ref], start_ref, end_ref


def matrix_profile(x, window):
    """Exact matrix profile of a time series, using the z-normalized Euclidean distance.

    Computed with stumpy.stump, which implements the exact STOMP algorithm with the default
    exclusion zone of ceil(window / 4) samples: a subsequence is never matched with its own
    immediate neighbours, which would otherwise be trivially similar.

    Arguments:
        x {ndarray} -- 1-D time series
        window {int} -- subsequence length, in samples

    Returns
    -------
    ndarray of size len(x) - window + 1 with the distance to the nearest neighbour
    """

    profile = stumpy.stump(np.asarray(x, dtype=float), m=int(window))

    return np.asarray(profile[:, 0], dtype=float)


def annotation_vector(gyr, jerk, window):
    """Annotation vector promoting the presence of the swing phase in the center of the window.

    For each candidate position, the normalized angular velocity and the normalized jerk are
    summed over the central third of the window. Both are maximal during the swing phase, so
    windows whose central third contains the swing phase are favoured. Standing periods, where
    both signals are almost zero, receive a value close to 0 and are therefore discarded, even
    though they are the most self-similar part of the recording.

    Arguments:
        gyr {ndarray} -- sagittal angular velocity
        jerk {ndarray} -- total jerk magnitude
        window {int} -- subsequence length, in samples

    Returns
    -------
    ndarray of size len(gyr) - window + 1 with values in [0, 1]
    """

    gyr_norm = np.abs(np.asarray(gyr, dtype=float))
    jerk_norm = np.abs(np.asarray(jerk, dtype=float))
    gyr_norm = gyr_norm / np.max(gyr_norm)
    jerk_norm = jerk_norm / np.max(jerk_norm)

    # cumulative sums allow the sliding sums to be computed in one pass
    cumsum = np.concatenate(([0.0], np.cumsum(gyr_norm + jerk_norm)))

    positions = np.arange(len(gyr_norm) - window + 1)
    start = positions + window // 3
    end = positions + (2 * window) // 3
    av = cumsum[end] - cumsum[start]

    max_av = np.max(av)
    if max_av <= 0:
        raise ValueError(
            "The annotation vector is identically zero: the recording does not seem to contain "
            "any movement.")

    return av / max_av


def corrected_matrix_profile(mp_values, av):
    """Correct a matrix profile with an annotation vector (Dau and Keogh, 2017).

    A subsequence with av = 0 gets a corrected value greater than or equal to the maximum of the
    matrix profile, and can therefore never be preferred to a subsequence with av = 1, whatever
    its repetition score. Between these two extremes the penalty is continuous and monotonic.

    Arguments:
        mp_values {ndarray} -- matrix profile
        av {ndarray} -- annotation vector, same size, values in [0, 1]

    Returns
    -------
    ndarray with the corrected matrix profile
    """

    mp_values = np.asarray(mp_values, dtype=float)

    if len(mp_values) != len(av):
        raise ValueError("The matrix profile and the annotation vector must have the same size.")

    finite = np.isfinite(mp_values)
    if not finite.any():
        raise ValueError("The matrix profile could not be computed on this recording.")

    max_mp = float(np.max(mp_values[finite]))
    mp_clean = np.where(finite, mp_values, max_mp)

    return mp_clean + (1 - av) * max_mp


def len_stride_estimation(data_1, data_2, freq):
    """Estimate the mean stride duration from both feet, computed from autocorrelations.
    The first foot is the foot of interest. The second foot acts as a safety net, in case the
    autocorrelation is limiting on the foot of interest.

    Arguments:
        data_1 {pandas Dataframe} -- dataframe with data from the foot sensor of interest
        data_2 {pandas Dataframe} -- dataframe with data from the foot sensor of the other side
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    int, float
        mean stride duration estimation, in samples
        value of the combined ACF at that lag
    """

    len_stride_data_1, autocorr_1 = len_stride_one_side(data_1, freq)
    len_stride_data_2, autocorr_2 = len_stride_one_side(data_2, freq)

    # no peak found on either side
    if len_stride_data_1 == 0 and len_stride_data_2 == 0:
        raise ValueError(
            "No stride duration could be estimated from the autocorrelation of either foot. "
            "The recording may be too short or may not contain a walking sequence.")

    # no peak found on one side: the other estimate is used
    if len_stride_data_1 == 0:
        return len_stride_data_2, autocorr_2
    if len_stride_data_2 == 0:
        return len_stride_data_1, autocorr_1

    # If the two estimates are too far apart, one of the peak detections is likely faulty: during
    # steady walking both should be nearly equal. The typical failure mode is that the first peak
    # is missed and a harmonic near 2*l is returned, hence the asymmetric test.
    if len_stride_data_1 / len_stride_data_2 >= 1.5:
        return len_stride_data_2, autocorr_2

    return len_stride_data_1, autocorr_1


def len_stride_one_side(data, freq):
    """Estimate the stride duration from one foot, as the first peak of a combined autocorrelation.

    Arguments:
        data {pandas Dataframe} -- dataframe with data from one of the foot sensors
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    int, float
        lag of the first peak of the combined ACF, in samples (0 if no peak is found)
        value of the combined ACF at that lag
    """

    acc_x = data["FreeAcc_X"].to_numpy(dtype=float)
    acc_y = data["FreeAcc_Y"].to_numpy(dtype=float)
    acc_z = data["FreeAcc_Z"].to_numpy(dtype=float)
    gyr_y = data["Gyr_Y"].to_numpy(dtype=float)

    # combined ACF from unbiased autocorrelations: the sagittal angular velocity carries half of
    # the total weight, the three acceleration axes sharing the other half
    acf = autocorr(gyr_y) / 2 + (autocorr(acc_x) + autocorr(acc_y) + autocorr(acc_z)) / 6

    # peaks are searched over the first quarter of the ACF only: the variance of the unbiased
    # estimator grows with the lag, since only N - i terms contribute to it
    acf_search = acf[:len(acf) // 4]
    index_pic = autocorr_indexes(acf_search, freq)

    if len(index_pic) > 0:
        return int(index_pic[0]), float(acf_search[index_pic[0]])

    return 0, 0


def autocorr_indexes(y, freq, thres=0.7, min_lag=0.8):
    """Find the local maxima of an autocorrelation function above an adaptive amplitude threshold.

    The threshold is built from the local dynamic range of the ACF beyond the minimum lag:
    T = lambda + thres * (y[i_max] - lambda), where y[i_max] is the highest ACF value for lags
    >= i_0 and lambda the lowest value between i_0 and i_max. It is therefore relative to the
    shape of the ACF rather than to an absolute amplitude, which makes it insensitive to the
    overall regularity of the gait.

    Arguments:
        y {ndarray} -- unbiased autocorrelation, already restricted to the lags of interest
        freq {int} -- acquisition frequency (Hz)
        thres {float} -- relative amplitude threshold in [0., 1.]
        min_lag {float} -- minimum lag, in seconds; shorter lags are discarded

    Returns
    -------
    ndarray with the indexes (lags, in samples) of the detected local maxima
    """

    y = np.asarray(y, dtype=float)

    i_0 = int(round(min_lag * freq))
    if len(y) <= i_0 + 1:
        return np.array([], dtype=int)

    # local dynamic range of the ACF beyond i_0
    i_max = i_0 + int(np.argmax(y[i_0:]))
    segment = y[i_0:i_max]
    lambda_min = float(np.min(segment)) if segment.size > 0 else float(y[i_0])
    threshold = lambda_min + thres * (float(y[i_max]) - lambda_min)

    peaks, _ = find_peaks(y, height=threshold)

    return peaks[peaks > i_0]


def autocorr(f):
    """Unbiased autocorrelation estimator.

    Computed through the Wiener-Khinchin theorem: the signal is zero-padded to 2N, its power
    spectrum is obtained by FFT, the inverse FFT is truncated to its first N values, divided
    term by term by (N - i), and finally rescaled so that acf[0] = 1.

    Arguments:
        f {ndarray} -- 1-D data to compute the autocorrelation from

    Returns
    -------
    acf {ndarray} -- array containing the unbiased autocorrelation
    """

    N = len(f)
    fvi = np.fft.fft(f, n=2 * N)
    acf = np.real(np.fft.ifft(fvi * np.conjugate(fvi))[:N])
    d = N - np.arange(N)
    acf = acf / d  # unbiased estimator
    acf = acf / acf[0]

    return acf
