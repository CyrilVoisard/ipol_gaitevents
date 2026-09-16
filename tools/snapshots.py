# New tests after review 
import sys, os, json, glob
import pandas as pd
from package import import_data, dtw_detection, quality

TAG = sys.argv[1]
FREQ = 100
DATASETS = sorted(glob.glob("example_datasets/*_rf.txt"))

os.makedirs("snapshots", exist_ok=True)
out = {}
for rf in DATASETS:
    name = os.path.basename(rf).replace("_rf.txt", "")
    lf = rf.replace("_rf.txt", "_lf.txt")
    data_rf = import_data.import_XSens(rf, FREQ)   # adapter après C7
    data_lf = import_data.import_XSens(lf, FREQ)
    steps = dtw_detection.steps_detection_full(data_rf, data_lf, FREQ, output="snapshots")
    steps = steps.sort_values(by=["Foot", "TO"]).reset_index(drop=True)
    out[name] = {
        "n_steps": int(len(steps)),
        "quality": int(quality.compute_quality(steps)),
        "events": steps[["Foot", "HO", "TO", "HS", "FF"]].astype(int).values.tolist(),
    }

with open(f"snapshots/{TAG}.json", "w") as f:
    json.dump(out, f, indent=1)
print(f"snapshots/{TAG}.json écrit")
