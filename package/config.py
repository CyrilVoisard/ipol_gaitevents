"""Numerical parameters of the MP-DTW algorithm.
"""

# ---------------------------------------------------------------------------
# Step 1: stride duration estimation
# ---------------------------------------------------------------------------

#: Weight of the sagittal angular velocity in the combined ACF, Equation (acf).
#: The three accelerations share the remaining weight (1 - W_GYR) equally.
W_GYR = 0.5

#: Relative amplitude threshold for the first ACF peak, in [0, 1].
ALPHA_ACF = 0.7

#: Minimum plausible stride duration, in seconds: lower lags are discarded.
MIN_LAG_S = 0.8

#: Fraction of the ACF searched for the first peak. The upper bound of the search range is
#: ACF_SEARCH_FRAC * N, where N is the number of samples.
ACF_SEARCH_FRAC = 0.25

#: Ratio above which the contralateral stride duration is preferred, Equation (ratio).
#: Halfway between agreement (1) and period doubling (2).
TWO_FOOT_RATIO = 1.5

# ---------------------------------------------------------------------------
# Step 2: reference stride extraction
# ---------------------------------------------------------------------------

#: Bounds of the window fraction used by the annotation vector, Equation (av).
#: The default (1/3, 2/3) is the central third of the candidate window.
AV_START_FRAC = 1.0 / 3.0
AV_END_FRAC = 2.0 / 3.0

# ---------------------------------------------------------------------------
# Steps 3 and 4: Dynamic Time Warping
# ---------------------------------------------------------------------------

#: Maximum slope of the Itakura parallelogram constraining the warping path.
ITAKURA_MAX_SLOPE = 2

# ---------------------------------------------------------------------------
# Step 4: gait events identification
# ---------------------------------------------------------------------------

#: Minimum correlation for a segment to be a stride candidate, Equation (candidates).
C_MIN = 0.35

#: Minimum ratio between the standard deviation of a candidate segment and that of the
#: reference stride, Equation (cost).
MU_MIN = 0.1

#: Non-maximum suppression radius for the candidate selection, as a fraction of the stride
#: duration.
NMS_FRAC = 0.5

# ---------------------------------------------------------------------------
# Pre-processing
# ---------------------------------------------------------------------------

#: Baseline (gravity) estimation window, in seconds from the beginning of the recording.
BASELINE_START_S = 0.0
BASELINE_END_S = 2.0

#: Butterworth low-pass filter: order and cut-off frequency (Hz).
FILTER_ORDER = 8
FILTER_FC = 14

#: Duration of the boxcar window smoothing the jerk, in seconds.
JERK_WIN_S = 0.09


def as_dict():
    """Return every parameter as a dictionary, to be logged alongside a run."""
    return {k: v for k, v in globals().items()
            if k.isupper() and not k.startswith("_")}
