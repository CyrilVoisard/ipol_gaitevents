import os
import json
import warnings
import numpy as np

from load_package import xsens, events, base_loaders
from semiogram_package import semio_plot
from rotagram_package import rotagram as rot
from gait_events_package import gait_events_detection as ge

from segmentation_package import uturn_detection_sal as u_sal
from segmentation_package import segmentation
from segmentation_package import visualisation as viz

from quality_package import general as quality


sensors_dict_path = "/Users/cyril/Library/Mobile Documents/com~apple~CloudDocs/Borelli/4 - Classes et fonctions/sensors_dictionnaires.json"
with open(sensors_dict_path) as json_dict:
    sensors_dict = json.load(json_dict)

gait_events_dict_path = "/Users/cyril/Library/Mobile Documents/com~apple~CloudDocs/Borelli/4 - Classes et fonctions/steps_dictionnaires.json"


class Experience_XS(object):
    # -----------------------------------------------------------------------------------------------
    # ---------------------------------- Initialisation ---------------------------------------------
    # -----------------------------------------------------------------------------------------------

    def __init__(self, root, id_exp, freq, segment=False, segment_method=None, disp=False):
        # ---------------------------------- Data from the name -------------------------------------
        self.root = root
        self.id_exp = id_exp
        self.freq = freq
        self.segment_method = segment_method or {"ge": "mpdtw", "uturn": "sal"}
        self.disp = disp

        # -------- Metadata parsing --------
        self.name = id_exp[0:6]
        self.patho = id_exp[7:10]
        self.date = id_exp[11:19] if len(id_exp) > 11 else "00000000"
        self.type = id_exp[20:26]

        # ---------------------------------- Data from the files ------------------------------------
        self.files_list = [f for f in os.listdir(root) if f.startswith(id_exp)]

        # raw data and processed data
        self.data_raw, self.t0, self.mode = self._load_raw_data()
        self.still_before = 200
        self.data = self._load_processed_data()

        # metadate
        self.metadata = self._load_metadata()

        # -------- Segmentation --------
        if segment:
            self._run_segmentation()

        # -------- Quality --------
        if segment:
            self._compute_quality()

        # -------- Plot --------
        if disp and segment:
            self._plot()

    # =========================================================
    # ---------------------- LOADERS ---------------------------
    # =========================================================

    def _find_file(self, sensor_name):
        matches = [f for f in self.files_list if f.endswith(sensor_name)]
        if not matches:
            raise FileNotFoundError(f"Sensor file not found: {sensor_name}")
        return os.path.join(self.root, matches[0])

    def _try_sensors(self, func):
        """Try all sensor configurations"""
        errors = []
        for sensors in sensors_dict:
            try:
                return func(sensors)
            except Exception as e:
                errors.append([sensors, str(e)])
        print("All sensor configurations failed:")
        print(self.files_list)
        for s, err in errors:
            print(f"{s} -> {err}")

        raise ValueError("No valid sensor configuration found")

    def _load_raw_data(self):
        def load(s):

            # print(self._find_file(s["tete"]))

            data = {
                "HE": xsens.load_xsens_raw(self._find_file(s["tete"]))[0],
                "LB": xsens.load_xsens_raw(self._find_file(s["tronc"]))[0],
                "RF": xsens.load_xsens_raw(self._find_file(s["pied_droit"]))[0],
                "LF": xsens.load_xsens_raw(self._find_file(s["pied_gauche"]))[0],
            }

            mode = xsens.load_xsens_raw(self._find_file(s["tete"]))[2]

            """
            data = {
                "HE": self.loader.load_raw(self._find_file(s["tete"]))[0],
                "LB": self.loader.load_raw(self._find_file(s["tronc"]))[0],
                "RF": self.loader.load_raw(self._find_file(s["pied_droit"]))[0],
                "LF": self.loader.load_raw(self._find_file(s["pied_gauche"]))[0],
            }"""

            t0 = xsens.find_t0(list(data.values()), mode=mode)

            return data, t0, mode

        return self._try_sensors(load)

    def _load_processed_data(self, disp=False):

        def load(s):
            data = {}
            freqs = []

            for key, label in zip(
                    ["HE", "LB", "RF", "LF"],
                    ["tete", "tronc", "pied_droit", "pied_gauche"],
            ):
                path = self._find_file(s[label])
                d, freq, default = xsens.import_xsens(
                    path,
                    t_delay=self.t0,
                    disp=disp,
                    mode=self.mode,
                )

                data[key] = d
                freqs.append(freq)

            if min(freqs) != max(freqs):
                raise ValueError(f"Sensors have different frequencies : {freqs}")
            else:
                self.freq = int(freqs[0])

            common_start = max(int(df["PacketCounter"].iloc[0]) for df in data.values())
            common_end = min(int(df["PacketCounter"].iloc[-1]) for df in data.values())

            if common_start > common_end:
                print(f"    [!] aucune plage de PacketCounter commune entre "
                      f"capteurs, epreuve ignoree")
                return None
            longest = max(len(df) for df in data.values())
            common_len = common_end - common_start + 1
            if longest - common_len > 10:
                print(f"    [!]  [i capteurs recadres sur la plage commune "
                      f"[{common_start}, {common_end}] ; rognage max {longest - common_len} ech.")

            for lab, df in data.items():
                mask = ((df["PacketCounter"] >= common_start)
                        & (df["PacketCounter"] <= common_end))
                data[lab] = df.loc[mask].reset_index(drop=True)

            return data

        return self._try_sensors(load)

    def _load_metadata(self):
        meta_files = [
            os.path.join(dp, f)
            for dp, _, files in os.walk(self.root)
            for f in files
            if "meta" in f and f.endswith(".json")
        ]

        if len(meta_files) != 1:
            if len(meta_files) == 0:
                print("Error: no file containing 'meta' and ending in .json found.")
            else:
                print("Error: multiple files containing 'meta' and ending in .json found.")
            return None

        with open(meta_files[0]) as f:
            return json.load(f)

    # =========================================================
    # ------------------- SEGMENTATION -------------------------
    # =========================================================

    def _run_segmentation(self):
        # ---- Gait events ----
        if self.metadata and "leftGaitEvents" in self.metadata:
            self.gait_events = events.ge_from_meta_to_events([
                self.metadata["leftGaitEvents"],
                self.metadata["rightGaitEvents"],
            ])
        else:
            self.gait_events, _ = ge.detect(
                self.data["RF"],
                self.data["LF"],
                freq=self.freq,
                method=self.segment_method["ge"],
            )

        # ---- U-turn ----
        self.uturn, self.uturn_val, self.r2, self.result, self.slope, self.uturn_side, self.n_uturns = u_sal.uturn_detection(
            self.data,
            freq=self.freq,
            show_plot=self.disp,
        )

        # ---- Segmentation ----
        self.seg = segmentation.uturn_to_seg(
            self.uturn,
            self.gait_events,
            method=self.segment_method["uturn"],
        )

        self.gait_events_corrected = ge.correct_detection(
            self.gait_events,
            self.seg,
        )

    # =========================================================
    # ---------------------- QUALITY ---------------------------
    # =========================================================

    def _compute_quality(self):
        self.quality_vec = quality.quality_vec(self)
        self.quality, self.quality_short = quality.quality_dict(
            self.quality_vec, self)

    # =========================================================
    # ----------------------- UTILS ----------------------------
    # =========================================================

    def _plot(self):
        import matplotlib.pyplot as plt

        viz.plot_segmentation_gait_events(
            self.gait_events_corrected,
            self.seg,
            self.data,
            self.freq,
            self.segment_method["ge"],
        )
        plt.show()

    # -----------------------------------------------------------------------------------------------
    # ---------------------------------- Semiogram --------------------------------------
    # -----------------------------------------------------------------------------------------------

    def get_semio_val(self, age=None):
        return semio_plot.compute_semio_val(self.data["LB"], self.seg, self.gait_events, age=age)

    def get_new_radar_semio(self, age=None):
        return semio_plot.new_radar(self.get_semio_val(), self.id_exp, age=age)

    # ---------------------------------- Rotagram ---------------------------------------------

    def plot_rotagram(self):
        rot.plot_rotagram(self.data["LB"], self.seg, self.gait_events, self.id_exp)

    def get_bornes(self):
        return segmentation.get_bornes(self.gait_events, self.seg)
