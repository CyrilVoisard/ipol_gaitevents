# Run the detection on the GAITRite validation database and produce the recall and precision rows.

# Signal pre-processing is left to the historical loader (xsens.import_XSens, through
# Experience_GR), so that the chain reproduces the one used for the published validation.

# Outputs, in RESULTS_DIR:
#   <tag>_recall.csv      one row per event annotated by the walkway
#   <tag>_precision.csv   one row per event returned by the algorithm
#   <tag>_runs.csv        one row per trial (status, duration, step count, ground truth source)
#   <tag>_config.json     the parameters used

import argparse
import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd
# --- paths to adapt -------------------------------------------------------------------

IPOL_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Historical modules used by the validation chain, shipped with this repository.
MODULES_DIR = os.path.join(IPOL_REPO, "validation_mpdtw")

# Folder holding input_validation_OK2 / _SEP / _EVF. Private data, not in the repository.
DATA_HOME = "/Users/cyril/Library/Mobile Documents/com~apple~CloudDocs/Borelli/6 - Publications/2024 - MPDTW/100 - Article final et suivi de soumission/0 - Données utilisées/Input/input"

RESULTS_DIR = os.path.join(IPOL_REPO, "results_sensitivity")

FREQ = 100

# --------------------------------------------------------------------------------------

sys.path.insert(0, IPOL_REPO)
sys.path.insert(0, MODULES_DIR)

from package import dtw_detection, config  # noqa: E402

import Patient as pt  # noqa: E402
import Patient_EVF as pt_evf  # noqa: E402
import Experience_GR as gr  # noqa: E402
import recall_precision  # noqa: E402


NAMES_RECALL = ["Patient", "Patient_Date", "Experience", "Pied", "Validate",
                "Event_tapis", "Value_tapis", "Event_detection", "Value_detection"]
NAMES_PRECISION = ["Patient", "Patient_Date", "Experience", "Pied", "Validate",
                   "Event_detection", "Value_detection", "Event_tapis", "Value_tapis"]
NAMES_RUNS = ["Cohort", "Patient", "Experience", "Status", "Ground_truth", "N_steps",
              "Duration_s", "Message"]


def to_legacy_format(steps_df):
    """Convert the step table to the layout expected by recall_precision.py.

    This repository returns [Foot, HO, TO, HS, FF, Score]. recall_precision.py, written for the
    original validation, reads index 3 as TO and index 4 as HS, which matches the former
    seven-column layout [Foot, Phase, HO, TO, HS, FF, Score]. The Phase column is reinserted.

    Arguments:
        steps_df {dataframe} -- step table returned by steps_detection_full

    Returns
    -------
    ndarray
        step table in the seven-column layout
    """

    a = steps_df.to_numpy()
    phase = np.full((len(a), 1), 100.0)

    return np.hstack([a[:, [0]], phase, a[:, 1:]])


def apply_overrides(overrides):
    """Apply parameter overrides of the form ALPHA_ACF=0.5 on the config module."""

    for item in overrides:
        if "=" not in item:
            raise ValueError("Invalid override: {} (expected NAME=value)".format(item))
        name, value = item.split("=", 1)
        name = name.strip()
        if not hasattr(config, name):
            raise ValueError("Unknown parameter in package/config.py: {}".format(name))
        setattr(config, name, float(value) if "." in value or "e" in value.lower() else int(value))
        print("    config.{} = {}".format(name, getattr(config, name)))


def list_experiments(cohort):
    """List the (patient, experiment name, folder) triplets holding a GAITRite recording.

    For the EVF cohort, Patient_EVF fills gr_exp from both the GR and the XSens folders, which
    duplicates every experiment. Only the GR entries are kept, since that folder holds both the
    sensor files and the CoP files.

    Arguments:
        cohort {str} -- OK2, SEP or EVF

    Returns
    -------
    list
        triplets (patient, experiment name, folder)
    """

    root = os.path.join(DATA_HOME, "input_validation_" + cohort)
    if not os.path.isdir(root):
        raise FileNotFoundError("Folder not found: {}".format(root))

    identifiers = sorted(d for d in os.listdir(root) if len(d) == 10)
    experiments = []
    seen = set()

    for id_patient in identifiers:
        if cohort == "EVF":
            patient = pt_evf.Patient(id_patient, root)
        else:
            patient = pt.Patient(id_patient, root)

        for id_exp in patient.gr_exp:
            if cohort == "EVF":
                rel_dir, exp_name = os.path.split(id_exp)
                if os.path.basename(rel_dir) != "GR":
                    continue
                path = os.path.join(patient.path, rel_dir)
            else:
                exp_name = id_exp
                path = patient.path

            key = (path, exp_name)
            if key in seen:
                continue
            seen.add(key)
            experiments.append((patient, exp_name, path))

    return experiments


def load_ground_truth(experience):
    """Load the signals and the ground truth, preferring the expert-reviewed walkway file.

    The _CoP_relu.txt file carries the two expert quality flags and the u-turn row that
    recall_precision.py relies on. When it is missing, the raw _CoP.txt is used instead and every
    event is marked as "Validated".

    Returns
    -------
    tuple
        right foot data, left foot data, ground truth array, source label
    """

    try:
        data_rf, data_lf, exo = experience.import_experience_pied()
        return data_rf, data_lf, exo, "reviewed"

    except (IndexError, FileNotFoundError):
        data_rf, data_lf, raw = experience.import_experience_pied(correction_gr=False)
        exo = np.array([list(row[:5]) + ["Validated", "Validated"] for row in raw], dtype=object)
        return data_rf, data_lf, exo, "raw"


def run_one(path, exp_name):
    """Detect the gait events of one trial.

    Returns
    -------
    tuple
        step table, ground truth, left foot data, right foot data, ground truth source
    """

    experience = gr.Experience_GR(path, exp_name, freq=FREQ)
    data_rf, data_lf, exo, source = load_ground_truth(experience)

    # output=0 skips the figures, which is required in batch mode
    steps = dtw_detection.steps_detection_full(data_rf, data_lf, FREQ, output=0)

    return to_legacy_format(steps), exo, data_lf, data_rf, source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True, help="name of the result series")
    parser.add_argument("--set", nargs="*", default=[], metavar="NAME=VALUE",
                        help="parameter override, for instance ALPHA_ACF=0.5")
    parser.add_argument("--cohorts", nargs="*", default=["OK2", "SEP", "EVF"])
    parser.add_argument("--max-per-patient", type=int, default=0,
                        help="limit the number of trials per patient (0 = all)")
    parser.add_argument("--require-relu", action="store_true",
                        help="skip the trials without a _CoP_relu.txt file")
    args = parser.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_recall = os.path.join(RESULTS_DIR, args.tag + "_recall.csv")
    if os.path.exists(out_recall):
        print("Tag '{}' already exists: {}".format(args.tag, out_recall))
        sys.exit(1)

    print("Parameters:")
    apply_overrides(args.set)

    with open(os.path.join(RESULTS_DIR, args.tag + "_config.json"), "w") as f:
        json.dump(config.as_dict(), f, indent=1)

    recall, precision, runs = [], [], []
    t_start = time.time()

    for cohort in args.cohorts:
        experiments = list_experiments(cohort)
        print("\n=== {}: {} trials ===".format(cohort, len(experiments)))

        seen = {}
        for patient, exp_name, path in experiments:
            if args.max_per_patient:
                seen[patient.id] = seen.get(patient.id, 0) + 1
                if seen[patient.id] > args.max_per_patient:
                    continue

            t0 = time.time()
            try:
                steps, exo, data_lf, data_rf, source = run_one(path, exp_name)

                if args.require_relu and source != "reviewed":
                    runs.append([cohort, patient.id, exp_name, "SKIPPED", source, 0,
                                 round(time.time() - t0, 1), "no _CoP_relu.txt"])
                    print("  {:<30} skipped (walkway not reviewed)".format(exp_name))
                    continue

                data = [data_lf, data_rf]
                recall = recall_precision.compute_recall(recall, data, steps, exo, exp_name)
                precision = recall_precision.compute_precision(precision, data, steps, exo, exp_name)

                runs.append([cohort, patient.id, exp_name, "OK", source, len(steps),
                             round(time.time() - t0, 1), ""])
                print("  {:<30} {:>3} steps   {:>5.1f} s   [{}]".format(
                    exp_name, len(steps), time.time() - t0, source))

            except Exception as exc:
                runs.append([cohort, patient.id, exp_name, "FAILED", "", 0,
                             round(time.time() - t0, 1), repr(exc)])
                print("  {:<30} FAILED: {}".format(exp_name, exc))
                traceback.print_exc(limit=2)

    pd.DataFrame(recall, columns=NAMES_RECALL).to_csv(out_recall, index=False)
    pd.DataFrame(precision, columns=NAMES_PRECISION).to_csv(
        os.path.join(RESULTS_DIR, args.tag + "_precision.csv"), index=False)
    pd.DataFrame(runs, columns=NAMES_RUNS).to_csv(
        os.path.join(RESULTS_DIR, args.tag + "_runs.csv"), index=False)

    n_ok = sum(1 for r in runs if r[3] == "OK")
    n_raw = sum(1 for r in runs if r[3] == "OK" and r[4] == "raw")
    print("\n{} trials visited, {} used, {} failed".format(
        len(runs), n_ok, sum(1 for r in runs if r[3] == "FAILED")))
    if n_raw:
        print("Warning: {} trials rely on a walkway file that was not reviewed.".format(n_raw))
    print("Total time: {:.0f} s. Results in {}".format(time.time() - t_start, RESULTS_DIR))


if __name__ == "__main__":
    main()