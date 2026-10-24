# Snapshot of the algorithm outputs, to compare results before and after a change.
#
# Set the GAIT_DATA environment variable to the folder holding the *_rf.txt and *_lf.txt files.
# main.py is run as a subprocess, exactly as the demo does, so the snapshots stay valid
# whatever the internal signatures become.

import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.environ.get("GAIT_DATA", os.path.join(REPO, "example_datasets"))
SNAP_DIR = os.path.join(REPO, "snapshots")
TIMEOUT_S = 900


def git_state():
    """Current commit and working tree state, to trace where a snapshot comes from.

    Returns
    -------
    dict
        commit, branch and dirty flag
    """

    def run(cmd):
        try:
            return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                                  check=True).stdout.strip()
        except Exception:
            return "unknown"

    return {"commit": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "dirty": bool(run(["git", "status", "--porcelain"]))}


def quality_from_events(events):
    """Recompute the alternation quality index from the detected events.

    Arguments:
        events {dict} -- LeftFootEvents and RightFootEvents, each entry being [TO, HS]

    Returns
    -------
    int
        quality index between 0 and 100
    """

    ic = [(hs, 0) for _, hs in events.get("LeftFootEvents", [])]
    ic += [(hs, 1) for _, hs in events.get("RightFootEvents", [])]
    ic.sort()

    if len(ic) < 2:
        return 0

    feet = [foot for _, foot in ic]
    n_alt = sum(abs(feet[k + 1] - feet[k]) for k in range(len(feet) - 1))

    return round(100 * n_alt / (len(feet) - 1))


def run_one(rf_path, lf_path, freq):
    """Run main.py on one trial, in a temporary working folder.

    Arguments:
        rf_path {str} -- right foot file
        lf_path {str} -- left foot file
        freq {int} -- acquisition frequency (Hz)

    Returns
    -------
    dict
        detected events and quality index, or an error message
    """

    work = tempfile.mkdtemp(prefix="gait_snap_")
    env = dict(os.environ, MPLBACKEND="Agg")

    cmd = [sys.executable, os.path.join(REPO, "main.py"),
           "-i0", os.path.abspath(rf_path),
           "-i1", os.path.abspath(lf_path),
           "-freq", str(freq)]

    try:
        proc = subprocess.run(cmd, cwd=work, env=env, capture_output=True,
                              text=True, timeout=TIMEOUT_S)

        if proc.returncode != 0:
            return {"error": "returncode {}".format(proc.returncode),
                    "stderr": proc.stderr[-2000:]}

        json_path = os.path.join(work, "gait_events.json")
        if not os.path.exists(json_path):
            return {"error": "gait_events.json not found", "stderr": proc.stderr[-2000:]}

        with open(json_path) as f:
            events = json.load(f)

        events.pop("Detection date", None)

        left = sorted([[int(a), int(b)] for a, b in events.get("LeftFootEvents", [])])
        right = sorted([[int(a), int(b)] for a, b in events.get("RightFootEvents", [])])

        return {"n_left": len(left), "n_right": len(right),
                "quality": quality_from_events(events),
                "left": left, "right": right}

    except subprocess.TimeoutExpired:
        return {"error": "timeout after {} s".format(TIMEOUT_S)}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    tag = sys.argv[1]
    freq = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    rf_files = sorted(glob.glob(os.path.join(DATA_DIR, "*_rf.txt")))
    if not rf_files:
        print("No *_rf.txt file found in {}".format(DATA_DIR))
        print("Set the GAIT_DATA environment variable to your data folder.")
        sys.exit(1)

    os.makedirs(SNAP_DIR, exist_ok=True)

    snap = {"tag": tag, "freq": freq, "git": git_state(),
            "python": sys.version.split()[0], "trials": {}}

    for rf in rf_files:
        name = os.path.basename(rf)[:-len("_rf.txt")]
        lf = rf[:-len("_rf.txt")] + "_lf.txt"

        if not os.path.exists(lf):
            print("  {:<30} left foot file missing, trial skipped".format(name))
            continue

        print("  {:<30} ...".format(name), end=" ", flush=True)
        result = run_one(rf, lf, freq)
        snap["trials"][name] = result

        if "error" in result:
            print("FAILED ({})".format(result["error"]))
        else:
            print("{} left steps, {} right steps, quality {}".format(
                result["n_left"], result["n_right"], result["quality"]))

    out_path = os.path.join(SNAP_DIR, tag + ".json")
    with open(out_path, "w") as f:
        json.dump(snap, f, indent=1)

    print("\nWritten: {}".format(out_path))
    if snap["git"]["dirty"]:
        print("Warning: the git working tree holds uncommitted changes.")


if __name__ == "__main__":
    main()