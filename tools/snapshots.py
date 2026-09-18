# New tests after review 

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
    """Commit courant et etat de l'arbre de travail, pour tracer d'ou vient l'instantane."""
    def run(cmd):
        try:
            return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                                  check=True).stdout.strip()
        except Exception:
            return "unknown"
 
    return {
        "commit": run(["git", "rev-parse", "--short", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(run(["git", "status", "--porcelain"])),
    }
 
 
def quality_from_events(events):
    """Recalcule l'indice de qualite (alternation droite / gauche) a partir des evenements.
 
    events : dict avec les cles LeftFootEvents et RightFootEvents, chaque entree etant [TO, HS].
    Reproduit quality.compute_quality : fusion des contacts initiaux (HS) des deux pieds,
    tri chronologique, proportion de paires consecutives qui alternent.
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
    """Lance main.py sur un essai, dans un dossier de travail temporaire."""
    work = tempfile.mkdtemp(prefix="gait_snap_")
    env = dict(os.environ, MPLBACKEND="Agg")  # pas d'affichage requis
 
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
            return {"error": "gait_events.json absent", "stderr": proc.stderr[-2000:]}
 
        with open(json_path) as f:
            events = json.load(f)
 
        # la date de detection change a chaque execution : elle doit sortir de l'instantane
        events.pop("Detection date", None)
 
        left = sorted([[int(a), int(b)] for a, b in events.get("LeftFootEvents", [])])
        right = sorted([[int(a), int(b)] for a, b in events.get("RightFootEvents", [])])
 
        return {
            "n_left": len(left),
            "n_right": len(right),
            "quality": quality_from_events(events),
            "left": left,
            "right": right,
        }
 
    except subprocess.TimeoutExpired:
        return {"error": "timeout apres {} s".format(TIMEOUT_S)}
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
        print("Aucun fichier *_rf.txt trouve dans {}".format(DATA_DIR))
        print("Definis la variable d'environnement GAIT_DATA vers ton dossier de donnees.")
        sys.exit(1)
 
    os.makedirs(SNAP_DIR, exist_ok=True)
 
    snap = {"tag": tag, "freq": freq, "git": git_state(),
            "python": sys.version.split()[0], "trials": {}}
 
    for rf in rf_files:
        name = os.path.basename(rf)[:-len("_rf.txt")]
        lf = rf[:-len("_rf.txt")] + "_lf.txt"
 
        if not os.path.exists(lf):
            print("  {:<30} fichier gauche manquant, essai ignore".format(name))
            continue
 
        print("  {:<30} ...".format(name), end=" ", flush=True)
        result = run_one(rf, lf, freq)
        snap["trials"][name] = result
 
        if "error" in result:
            print("ECHEC ({})".format(result["error"]))
        else:
            print("{} pas gauche, {} pas droit, qualite {}".format(
                result["n_left"], result["n_right"], result["quality"]))
 
    out_path = os.path.join(SNAP_DIR, tag + ".json")
    with open(out_path, "w") as f:
        json.dump(snap, f, indent=1)
 
    print("\nEcrit : {}".format(out_path))
    if snap["git"]["dirty"]:
        print("Attention : l'arbre de travail git contient des modifications non commitees.")
 
 
if __name__ == "__main__":
    main()
