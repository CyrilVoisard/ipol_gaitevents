import os
import pandas as pd
import numpy as np
from functools import lru_cache
from scipy.signal import butter, filtfilt
from scipy import interpolate, stats
import re

pd.options.mode.chained_assignment = None  # default='warn'


# =========================================================
# ---------------------- UTILS -----------------------------
# =========================================================

def find_t0(data_list, mode="auto"):
    if mode == "mtw":
        return min(data["PacketCounter"].iloc[0] for data in data_list)
    elif mode == "dot":
        return max(data["Time"].iloc[0] for data in data_list)
    else:
        raise ValueError(f"Mode inconnu : {mode}")


def interpolate_signal(x, y, new_x):
    f = interpolate.interp1d(
        x,
        y,
        bounds_error=False,
        fill_value="extrapolate"
    )
    return f(new_x)


def quaternion_to_gyration(data):
    data['Gyr_X'], data['Gyr_Y'], data['Gyr_Z'] = quaternion_to_euler_angle_vectorized(data['dq_W'], data['dq_X'],
                                                                                       data['dq_Y'], data['dq_Z'])
    data.rename(columns={'dv[1]': 'Acc_X', 'dv[2]': 'Acc_Y', 'dv[3]': 'Acc_Z', }, inplace=True)
    return data


def quaternion_to_euler_angle_vectorized(w, x, y, z):
    ysqr = y * y
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + ysqr)
    x = np.degrees(np.arctan2(t0, t1))
    t2 = +2.0 * (w * y - z * x)
    t2 = np.clip(t2, a_min=-1.0, a_max=1.0)
    y = np.degrees(np.arcsin(t2))
    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (ysqr + z * z)
    z = np.degrees(np.arctan2(t3, t4))
    return x, y, z


# =========================================================
# ------------------ OUTLIERS ------------------------------
# =========================================================

def detect_outliers_zscore(series, threshold=10):
    z_scores = np.abs(stats.zscore(series))
    return z_scores > threshold


def clean_sensor_data(df, columns, threshold=10, disp=False):
    df_clean = df.copy()
    outliers_count = {}

    for col in columns:
        outliers = detect_outliers_zscore(df_clean[col], threshold)
        outliers_count[col] = np.sum(outliers)

        if disp and outliers_count[col] > 0:
            print(f"{col}: {outliers_count[col]} outliers")

        # interpolation au lieu de suppression
        df_clean.loc[outliers, col] = np.nan
        df_clean[col] = df_clean[col].interpolate(limit_direction="both")

    return df_clean, outliers_count



# =========================================================
# --------- LECTURE MTW : cache + moteur C ----------------
# =========================================================
# Un essai fait parser les mêmes fichiers 9 fois (4 capteurs bruts + la tête une
# 5e fois pour 'mode', puis 4 relectures via import_xsens). Le parsing brut est
# identique dans les deux chemins : on le fait une fois et on rend une copie.
# Le moteur C remplace le moteur python (regex ligne à ligne) sans changer le
# découpage : le séparateur retenu est celui qui donne le même nombre de champs
# que la regex [,\s] d'origine sur la ligne d'en-tête.

def _mtw_header(filename):
    """(index de la ligne d'en-tête, ligne d'en-tête, fréquence)."""
    freq = None
    with open(filename, 'r', errors='replace') as f:
        for i, line in enumerate(f):
            if freq is None:
                match = re.search(r'(\d+\.?\d*)Hz', line)
                if match:
                    freq = float(match.group(1))
            if line.startswith('PacketCounter'):
                return i, line, freq
    raise ValueError("PacketCounter non trouvé.")


def _mtw_sep(header):
    ref = len([t for t in re.split(r"[,\s]", header.rstrip("\n")) if t != ""])
    for sep in (",", "\t", ";"):
        if len([t for t in header.rstrip("\n").split(sep) if t != ""]) == ref:
            return sep
    return r"\s+"


@lru_cache(maxsize=8)
def _mtw_parse_cached(filename, mtime_ns, size):
    i, header, freq = _mtw_header(filename)
    return pd.read_csv(filename, skiprows=i, sep=_mtw_sep(header), engine='c'), freq


def read_mtw(filename):
    """Parsing brut d'un fichier MTw, mis en cache (clé : chemin + date + taille)."""
    st = os.stat(filename)
    df, freq = _mtw_parse_cached(str(filename), st.st_mtime_ns, st.st_size)
    return df.copy(), freq


# =========================================================
# ---------------------- LOADERS -------------------------
# =========================================================

def load_xsens_raw(filename, disp=False, mode="auto"):
    """
    Wrapper pour charger MTw ou DOT
    """

    if mode == "auto":
        # détection simple basée sur colonnes
        with open(filename, 'r') as f:
            header = f.read(1000)

        if "PacketCounter" in header:
            mode = "mtw"
        elif "Time" in header:
            mode = "dot"
        else:
            raise ValueError("Impossible de détecter le type de fichier", header)

    if mode == "mtw":
        signal, freq = load_xsens_mtw_raw(filename, disp)
        return signal, freq, mode

    elif mode == "dot":
        signal, freq = load_xsens_dot_raw(filename, disp)
        return signal, freq, mode

    else:
        raise ValueError(f"Mode inconnu : {mode}")


def load_xsens_mtw_raw(filename, disp=False):
    return read_mtw(filename)


def load_xsens_dot_raw(filename, disp=False):
    df = pd.read_csv(filename)

    # --- Standardisation noms colonnes ---
    df = df.rename(columns={
        "SampleTimeFine": "Time",
        "acc_x": "Acc_X",
        "acc_y": "Acc_Y",
        "acc_z": "Acc_Z",
        "gyr_x": "Gyr_X",
        "gyr_y": "Gyr_Y",
        "gyr_z": "Gyr_Z",
    })

    # --- Sécurité ---
    required = ["Time", "Acc_X", "Acc_Y", "Acc_Z", "Gyr_X", "Gyr_Y", "Gyr_Z"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Colonnes manquantes : {missing}")

    time = df["Time"].values - df["Time"].iloc[0]
    dt = np.median(np.diff(time))
    df["PacketCounter"] = np.round(time / dt, 0)

    print("raw freq : ", 1/dt)

    return df, int(round(1 / dt, 0))


def load_xsens(filename, t_delay=None, disp=False, mode="auto"):
    """
    Wrapper pour charger MTw ou DOT
    """

    if mode == "auto":
        # détection simple basée sur colonnes
        with open(filename, 'r') as f:
            header = f.read(500)

        if "PacketCounter" in header:
            mode = "mtw"
        elif "Time" in header:
            mode = "dot"
        else:
            raise ValueError("Impossible de détecter le type de fichier")

    if mode == "mtw":
        if t_delay is None:
            t_delay = np.inf
        return load_xsens_mtw(filename, t0_min=t_delay, disp=disp)

    elif mode == "dot":
        if t_delay is None:
            t_delay = 0
        return load_xsens_dot(filename, t_start=t_delay, disp=disp)

    else:
        raise ValueError(f"Mode inconnu : {mode}")


def load_xsens_mtw(filename, t0_min=np.inf, disp=False):
    # -------- Read header + load data (cache + moteur C) --------
    signal, freq = read_mtw(filename)

    # -------- Clean outliers --------
    sensor_cols = [c for c in signal.columns if c.startswith(('Acc_', 'Gyr_', 'Mag_'))]
    signal, outliers_count = clean_sensor_data(signal, sensor_cols, disp=disp)

    if disp and any(outliers_count.values()):
        print("Outliers corrigés :", outliers_count)

    # -------- PacketCounter --------
    # copie explicite : avec pandas récent, .values peut être en lecture seule
    packet_counter = np.array(signal["PacketCounter"].values, copy=True)

    # -------- Overflow correction --------
    if np.any(np.diff(packet_counter) < 0):
        if disp:
            print("Overflow détecté")
        overflow_idx = np.where(np.diff(packet_counter) < 0)[0][0] + 1
        packet_counter[overflow_idx:] += 65536
        signal["PacketCounter"] = packet_counter   # l'original modifiait la colonne en place

    # -------- Delay correction --------
    if t0_min != np.inf:
        delay = packet_counter[0] - t0_min
    else:
        delay = 0

    packet_counter_corr = packet_counter - delay

    # -------- Detect missing packets --------
    diff = np.diff(packet_counter_corr)
    missing = np.any((diff != 1) & (diff != -65535))

    # -------- Interpolation if needed --------
    if missing or delay != 0:
        if disp:
            print("Interpolation appliquée")

        new_time = np.arange(packet_counter_corr[0],
                             packet_counter_corr[-1] + 1)

        interpolated_data = {}

        for col in signal.columns:
            interpolated_data[col] = interpolate_signal(
                packet_counter_corr,
                signal[col].values,
                new_time
            )

        signal = pd.DataFrame(interpolated_data)
        signal["PacketCounter"] = np.arange(len(signal))

    else:
        signal["PacketCounter"] = np.arange(len(signal))

    # -------- Time vector --------
    if freq is not None:
        signal["Time"] = signal["PacketCounter"] / freq

    return signal, freq


def load_xsens_dot(filename, t_start=0, disp=False):
    freq = 0

    signal_raw, freq_raw = load_xsens_dot_raw(filename)
    signal = signal_raw[signal_raw["Time"] >= t_start].reset_index(drop=True)
    if disp:
        print(f"Suppression des données avant {t_start} ms : {len(signal)} échantillons restants sur {len(signal_raw)}")

    # -------- Clean outliers --------
    sensor_cols = [c for c in signal.columns if c.startswith(('Acc_', 'Gyr_', 'Mag_'))]
    signal, outliers_count = clean_sensor_data(signal, sensor_cols, disp=disp)

    if disp and any(outliers_count.values()):
        print("Outliers corrigés :", outliers_count)

    # -------- Delay correction --------
    time = signal["Time"].values
    time_corr = time - t_start

    # fréquence réelle
    dt = np.median(np.diff(time_corr))
    i = 1
    while freq < 40:
        i = i * 10
        freq = round(i / dt, 0)
    if disp:
        print(f"Fréquence estimée : {freq} Hz, facteur de correction : {i}, fréquence brute : {freq_raw} Hz")
    packet_counter = np.round(time_corr * freq / i, 0)

    # -------- Overflow correction --------
    if np.any(np.diff(packet_counter) < 0):
        print("Overflow détecté")
        overflow_idx = np.where(np.diff(packet_counter) < 0)[0][0] + 1
        print(f"Index d'overflow : {overflow_idx}, PacketCounter avant : {packet_counter[overflow_idx-1]}, après : {packet_counter[overflow_idx]}")

    # -------- Detect missing packets --------
    diff = np.diff(packet_counter)
    missing = np.any((diff != 1) & (diff != -65535))

    # -------- Interpolation if needed --------
    if missing :
        if disp:
            print("Interpolation appliquée")

        new_time = np.arange(packet_counter[0],
                             packet_counter[-1] + 1)

        interpolated_data = {}

        for col in signal.columns:
            interpolated_data[col] = interpolate_signal(
                packet_counter,
                signal[col].values,
                new_time
            )

        signal = pd.DataFrame(interpolated_data)
        signal["PacketCounter"] = np.arange(len(signal))

    else:
        signal["PacketCounter"] = packet_counter

    # -------- Time vector --------
    if freq is not None:
        signal["Time"] = signal["PacketCounter"] / freq

    return signal, freq


# =========================================================
# ---------------------- PROCESSING MTW -------------------
# =========================================================

def import_xsens(
        path,
        t_delay=None,
        correction_acc=True,
        start=0,
        end=2,
        filtre=True,
        ordre=8,
        fc=14,
        freq_default=100,
        disp=False,
        mode="auto"
):
    data, freq = load_xsens(path, t_delay, disp=disp, mode=mode)

    default = False
    if freq is None:
        freq = freq_default
        default = True

    end = int(end * freq)

    # -------- Remove bias --------
    if correction_acc:
        for axis in ["X", "Y", "Z"]:
            data[f"FreeAcc_{axis}"] = data[f"Acc_{axis}"] - np.mean(data[f"Acc_{axis}"][start:end])
            data[f"Gyr_{axis}"] = data[f"Gyr_{axis}"] - np.mean(data[f"Gyr_{axis}"][start:end])

    # -------- Filtering --------
    if filtre:
        for sig in ["Acc", "FreeAcc", "Gyr"]:
            data = filtre_sig(data, sig, ordre, fc, freq)

    return data, freq, default


# =========================================================
# ---------------------- FILTERS ---------------------------
# =========================================================

def filtre_sig(data, prefix, ordre, fc, freq):
    for axis in ["X", "Y", "Z"]:
        col = f"{prefix}_{axis}"
        data[col] = data[col].fillna(0)
        data[col] = filtre_passe_bas(data[col], ordre, fc, freq)
    return data


def filtre_passe_bas(sig, ordre=8, fc=14, fe=100):
    if fc >= fe / 2:
        raise ValueError("Cutoff frequency must be < Nyquist")

    b, a = butter(N=ordre, Wn=fc / (fe / 2), btype='low')
    return filtfilt(b, a, sig)
