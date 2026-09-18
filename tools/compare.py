
import json
import os
import sys
 
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP_DIR = os.path.join(REPO, "snapshots")
 
 
def load(tag):
    path = os.path.join(SNAP_DIR, tag + ".json")
    if not os.path.exists(path):
        print("Instantane introuvable : {}".format(path))
        sys.exit(2)
    with open(path) as f:
        return json.load(f)
 
 
def median(values):
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    return float(s[n // 2]) if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0
 
 
def diff_side(before, after):
    """Ecarts evenement par evenement, pour un pied, si les comptages concordent.
 
    Retourne (n_events, n_identical, max_diff, median_diff) ou None si les comptages different.
    """
    if len(before) != len(after):
        return None
 
    diffs = []
    for (to_b, hs_b), (to_a, hs_a) in zip(before, after):
        diffs.append(abs(to_a - to_b))
        diffs.append(abs(hs_a - hs_b))
 
    n_identical = sum(1 for d in diffs if d == 0)
 
    return len(diffs), n_identical, (max(diffs) if diffs else 0), median(diffs)
 
 
def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
 
    tag_a, tag_b = sys.argv[1], sys.argv[2]
    as_md = "--md" in sys.argv
 
    a, b = load(tag_a), load(tag_b)
 
    print("AVANT : {:<12} commit {} ({})".format(tag_a, a["git"]["commit"], a["git"]["branch"]))
    print("APRES : {:<12} commit {} ({})".format(tag_b, b["git"]["commit"], b["git"]["branch"]))
    if a.get("freq") != b.get("freq"):
        print("Attention : les deux instantanes n'ont pas la meme frequence "
              "({} contre {} Hz).".format(a.get("freq"), b.get("freq")))
    print()
 
    rows = []
    all_identical = True
 
    names = sorted(set(a["trials"]) | set(b["trials"]))
    for name in names:
        ta, tb = a["trials"].get(name), b["trials"].get(name)
 
        if ta is None or tb is None:
            print("{:<28} present dans un seul instantane".format(name))
            all_identical = False
            continue
 
        if "error" in ta or "error" in tb:
            print("{:<28} ECHEC  avant={}  apres={}".format(
                name, ta.get("error", "ok"), tb.get("error", "ok")))
            if "error" in tb and "error" not in ta:
                print("        stderr : {}".format(tb.get("stderr", "")[-400:]))
            all_identical = False
            continue
 
        n_before = ta["n_left"] + ta["n_right"]
        n_after = tb["n_left"] + tb["n_right"]
 
        print("{:<28} pas {} -> {}   qualite {} -> {}".format(
            name, n_before, n_after, ta["quality"], tb["quality"]))
 
        stats = []
        for side in ("left", "right"):
            d = diff_side(ta[side], tb[side])
            if d is None:
                print("        {:<6} nombre de pas different "
                      "({} -> {}), comparaison evenement par evenement impossible".format(
                          side, len(ta[side]), len(tb[side])))
                all_identical = False
                stats = None
                break
            stats.append(d)
 
        if stats is None:
            rows.append((name, n_before, n_after, ta["quality"], tb["quality"], "n/a", "n/a", "n/a"))
            continue
 
        n_tot = sum(s[0] for s in stats)
        n_id = sum(s[1] for s in stats)
        d_max = max(s[2] for s in stats)
        d_med = median([s[3] for s in stats])
        pct_id = 100.0 * n_id / n_tot if n_tot else 100.0
 
        if d_max == 0:
            print("        identique (tous les evenements au meme echantillon)")
        else:
            all_identical = False
            print("        ecart sur TO et HS : max {} ech., median {:.1f} ech., "
                  "{:.1f} % d'evenements inchanges".format(d_max, d_med, pct_id))
 
        rows.append((name, n_before, n_after, ta["quality"], tb["quality"],
                     d_max, "{:.1f}".format(d_med), "{:.1f} %".format(pct_id)))
 
    print()
    print("VERDICT : {}".format("IDENTIQUE" if all_identical else "DIFFERENT"))
 
    if as_md:
        print("\n| Essai | Pas avant | Pas apres | Qualite avant | Qualite apres | "
              "Ecart max (ech.) | Ecart median (ech.) | Evenements inchanges |")
        print("|---|---|---|---|---|---|---|---|")
        for r in rows:
            print("| {} | {} | {} | {} | {} | {} | {} | {} |".format(*r))
 
    sys.exit(0 if all_identical else 1)
 
 
if __name__ == "__main__":
    main()
 
