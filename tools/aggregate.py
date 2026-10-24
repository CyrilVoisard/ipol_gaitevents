# Aggregate the outputs of tools/run_validation.py into per-cohort metrics.
#
# Follows the methodology of the original validation study: an annotated event is correctly
# detected when the nearest detected event of the same kind lies within 20 % of the stride
# duration, this tolerance being applied by recall_precision.py through its Validate flag.
#
# Validate codes:
#    1    annotated "Validated" and correctly detected
#    0.5  annotated "Not_validated" and correctly detected
#    0    no detection within the tolerance
#   -1    event discarded by the expert ("Deleted")      -> excluded
#   -2    event too close to the trial boundaries        -> excluded
#   -3    event inside the u-turn                        -> excluded
#

import argparse
import os
import sys

import numpy as np
import pandas as pd

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_sensitivity")

COHORT_OF = {"OK2": "CS", "SEP": "MS", "PVE": "EVF", "EVF": "EVF"}

# Set to False for the strict variant, where "Not_validated" events do not count as correct.
COUNT_HALF_AS_CORRECT = True


def cohort_from_row(patient_date):
    """Infer the cohort from the experiment identifier.

    Arguments:
        patient_date {str} -- identifier of the form NAMEabc_COH_date

    Returns
    -------
    str
        CS, MS, EVF or ?
    """

    for key, label in COHORT_OF.items():
        if key in str(patient_date):
            return label
    return "?"


def load(tag):
    rec = pd.read_csv(os.path.join(RESULTS_DIR, tag + "_recall.csv"))
    pre = pd.read_csv(os.path.join(RESULTS_DIR, tag + "_precision.csv"))
    for df in (rec, pre):
        df["Cohort"] = df["Patient_Date"].map(cohort_from_row)
        df["Validate"] = pd.to_numeric(df["Validate"], errors="coerce")
        df["Value_tapis"] = pd.to_numeric(df["Value_tapis"], errors="coerce")
        df["Value_detection"] = pd.to_numeric(df["Value_detection"], errors="coerce")
    return rec, pre


def metrics(rec, pre, freq=100):
    """Compute recall, precision, F1 and median absolute errors for one cohort.

    Arguments:
        rec {dataframe} -- recall rows
        pre {dataframe} -- precision rows
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    dict or None
        metrics, None when the cohort holds no evaluable event
    """

    ok = {1.0} | ({0.5} if COUNT_HALF_AS_CORRECT else set())

    rec_eval = rec[rec["Validate"] >= 0]
    pre_eval = pre[pre["Validate"] >= 0]

    n_annotated = len(rec_eval)
    n_detected = len(pre_eval)
    if n_annotated == 0 or n_detected == 0:
        return None

    recall = rec_eval["Validate"].isin(ok).mean()
    precision = pre_eval["Validate"].isin(ok).mean()
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) else 0.0

    pairs = rec_eval[rec_eval["Validate"].isin(ok)].copy()
    pairs["delta_ms"] = (pairs["Value_tapis"] - pairs["Value_detection"]).abs() * 1000.0 / freq

    out = {"N_events": n_annotated,
           "Recall_%": 100 * recall,
           "Precision_%": 100 * precision,
           "F1_%": 100 * f1}

    for event in ("HS", "TO"):
        sub = pairs[pairs["Event_tapis"] == event]["delta_ms"]
        if len(sub):
            out["Delta{}_med_ms".format(event)] = np.median(sub)
            out["Delta{}_q1_ms".format(event)] = np.percentile(sub, 25)
            out["Delta{}_q3_ms".format(event)] = np.percentile(sub, 75)

    return out


def table_for_tag(tag):
    rec, pre = load(tag)
    rows = []
    for cohort in ["CS", "MS", "EVF"]:
        m = metrics(rec[rec["Cohort"] == cohort], pre[pre["Cohort"] == cohort])
        if m:
            rows.append(dict(Tag=tag, Cohort=cohort, **m))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tags", nargs="+")
    parser.add_argument("--latex", action="store_true", help="print a LaTeX tabular")
    args = parser.parse_args()

    tables = []
    for tag in args.tags:
        try:
            tables.append(table_for_tag(tag))
        except FileNotFoundError:
            print("Tag not found: {}".format(tag), file=sys.stderr)

    if not tables:
        sys.exit(1)

    full = pd.concat(tables, ignore_index=True)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print(full.round(1).to_string(index=False))

    out = os.path.join(RESULTS_DIR, "summary.csv")
    full.to_csv(out, index=False)
    print("\nWritten: {}".format(out))

    if args.latex:
        piv = full.pivot(index="Tag", columns="Cohort", values="F1_%").round(1)
        print("\n% --- for the sensitivity table ---")
        print(piv.to_latex())


if __name__ == "__main__":
    main()