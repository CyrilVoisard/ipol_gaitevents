import sys, json, numpy as np
a = json.load(open(f"snapshots/{sys.argv[1]}.json"))
b = json.load(open(f"snapshots/{sys.argv[2]}.json"))
for k in a:
    ea, eb = np.array(a[k]["events"]), np.array(b[k]["events"])
    print(f"\n{k}: {a[k]['n_steps']} -> {b[k]['n_steps']} pas, "
          f"qualité {a[k]['quality']} -> {b[k]['quality']}")
    if ea.shape == eb.shape:
        d = np.abs(ea[:, 1:] - eb[:, 1:])
        print(f"  écart HO/TO/HS/FF (échantillons) : max={d.max()}, médian={np.median(d):.1f}, "
              f"identiques={100 * (d == 0).mean():.1f}%")
    else:
        print("  !! nombre de pas différent, comparaison ligne à ligne impossible")
