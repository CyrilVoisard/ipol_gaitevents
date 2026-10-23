"""Agrégation des sorties de run_validation.py en métriques par cohorte.

Reproduit la méthodologie de l'article JNER 2024 :
  - un événement annoté est correctement détecté si l'écart au plus proche événement détecté
    de même nature est inférieur à 20 % de la durée de foulée (c'est recall_precision.py qui
    applique cette tolérance et pose l'indicateur Validate) ;
  - recall    = proportion d'événements annotés correctement détectés ;
  - precision = proportion d'événements détectés correspondant à un événement annoté ;
  - F1        = moyenne harmonique des deux ;
  - accuracy  = médiane de l'écart absolu, sur les seules paires correctement appariées.

Codes de l'indicateur Validate produits par recall_precision.py :
     1   événement annoté "Validated" et correctement détecté
     0.5 événement annoté "Not_validated" et correctement détecté
     0   pas de détection dans la tolérance
    -1   événement supprimé par l'expert ("Deleted")     -> exclu
    -2   événement trop proche des bords de l'essai       -> exclu
    -3   événement situé dans le demi-tour                -> exclu

Usage
-----
    python aggregate.py baseline
    python aggregate.py baseline alpha_0.5 alpha_0.9 --latex
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_sensitivity")

# préfixe de l'identifiant patient -> cohorte
COHORT_OF = {"OK2": "CS", "SEP": "MS", "PVE": "EVF", "EVF": "EVF"}

#: si True, les événements annotés "Not_validated" mais détectés (Validate = 0.5) comptent
#: comme correctement détectés. Mets False pour la variante stricte.
COUNT_HALF_AS_CORRECT = True


def cohort_from_row(patient_date):
    """Déduit la cohorte de l'identifiant d'expérience (NOMabc_COH_date...)."""
    for key, label in COHORT_OF.items():
        if key in str(patient_date):
            return label
    return "?"


def load(tag):
    rec = pd.read_csv(os.path.join(RESULTS_DIR, tag + "_recall.csv"))
    pre = pd.read_csv(os.path.join(RESULTS_DIR, tag + "_precision.csv"))
    for df in (rec, pre):
        df["Cohorte"] = df["Patient_Date"].map(cohort_from_row)
        df["Validate"] = pd.to_numeric(df["Validate"], errors="coerce")
        df["Value_tapis"] = pd.to_numeric(df["Value_tapis"], errors="coerce")
        df["Value_detection"] = pd.to_numeric(df["Value_detection"], errors="coerce")
    return rec, pre


def metrics(rec, pre, freq=100):
    """Calcule les métriques pour une cohorte."""
    ok = {1.0} | ({0.5} if COUNT_HALF_AS_CORRECT else set())

    rec_eval = rec[rec["Validate"] >= 0]
    pre_eval = pre[pre["Validate"] >= 0]

    n_annotes = len(rec_eval)
    n_detectes = len(pre_eval)
    if n_annotes == 0 or n_detectes == 0:
        return None

    recall = rec_eval["Validate"].isin(ok).mean()
    precision = pre_eval["Validate"].isin(ok).mean()
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) else 0.0

    # écarts, en millisecondes, sur les seules paires correctement appariées
    paires = rec_eval[rec_eval["Validate"].isin(ok)].copy()
    paires["delta_ms"] = (paires["Value_tapis"] - paires["Value_detection"]).abs() * 1000.0 / freq

    out = {"N_evenements": n_annotes,
           "Recall_%": 100 * recall,
           "Precision_%": 100 * precision,
           "F1_%": 100 * f1}

    for event, label in (("HS", "HS"), ("TO", "TO")):
        sub = paires[paires["Event_tapis"] == event]["delta_ms"]
        if len(sub):
            out["Delta{}_med_ms".format(label)] = np.median(sub)
            out["Delta{}_q1_ms".format(label)] = np.percentile(sub, 25)
            out["Delta{}_q3_ms".format(label)] = np.percentile(sub, 75)

    return out


def table_for_tag(tag):
    rec, pre = load(tag)
    rows = []
    for cohorte in ["CS", "MS", "EVF"]:
        m = metrics(rec[rec["Cohorte"] == cohorte], pre[pre["Cohorte"] == cohorte])
        if m:
            rows.append(dict(Tag=tag, Cohorte=cohorte, **m))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tags", nargs="+")
    parser.add_argument("--latex", action="store_true", help="sortie au format tabular LaTeX")
    args = parser.parse_args()

    tables = []
    for tag in args.tags:
        try:
            tables.append(table_for_tag(tag))
        except FileNotFoundError:
            print("Tag introuvable : {}".format(tag), file=sys.stderr)

    if not tables:
        sys.exit(1)

    full = pd.concat(tables, ignore_index=True)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print(full.round(1).to_string(index=False))

    out = os.path.join(RESULTS_DIR, "synthese.csv")
    full.to_csv(out, index=False)
    print("\nEcrit : {}".format(out))

    if args.latex:
        piv = full.pivot(index="Tag", columns="Cohorte", values="F1_%").round(1)
        print("\n% --- à coller dans le tableau de sensibilité ---")
        print(piv.to_latex())


if __name__ == "__main__":
    main()
