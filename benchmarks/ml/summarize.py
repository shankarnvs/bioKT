"""Summarize held-out model parity and the fixed human protein panel."""

import csv
import json
import statistics
import sys
from pathlib import Path

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def summarize_iris(rows):
    summary = {}
    for algorithm in ("naive_bayes", "linear"):
        selected = [r for r in rows if r["model"].startswith(f"iris-{algorithm}-")]
        by_seed = {}
        for seed in sorted({r["seed"] for r in selected}):
            fold = [r for r in selected if r["seed"] == seed]
            expected = [r["expected"] for r in fold]
            observed = [r["biokt"] for r in fold]
            by_seed[seed] = {
                "cases": len(fold),
                "accuracy": accuracy_score(expected, observed),
                "macro_f1": f1_score(expected, observed, average="macro"),
            }
        summary[algorithm] = {
            "predictions": len(selected),
            "label_matches": sum(r["reference"] == r["biokt"] for r in selected),
            "max_probability_error": max(float(r["max_probability_error"]) for r in selected),
            "mean_accuracy": statistics.mean(row["accuracy"] for row in by_seed.values()),
            "mean_macro_f1": statistics.mean(row["macro_f1"] for row in by_seed.values()),
            "by_seed": by_seed,
        }
    return summary


def summarize_proteins(rows):
    expected = [r["expected"] == "membrane" for r in rows]
    observed = [r["biokt"] == "membrane" for r in rows]
    probabilities = [float(r["probability_membrane"]) for r in rows]
    tn, fp, fn, tp = confusion_matrix(expected, observed, labels=[False, True]).ravel()
    return {
        "cases": len(rows), "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        "accuracy": accuracy_score(expected, observed),
        "balanced_accuracy": balanced_accuracy_score(expected, observed),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "auroc": roc_auc_score(expected, probabilities),
        "brier": brier_score_loss(expected, probabilities),
        "balanced_chance_accuracy": 0.5,
        "balanced_constant_probability_brier": 0.25,
    }


def main():
    folder = Path(sys.argv[1])
    result = {
        "iris": summarize_iris(read_tsv(folder / "iris-results.tsv")),
        "proteins": summarize_proteins(read_tsv(folder / "protein-results.tsv")),
    }
    content = json.dumps(result, indent=2)
    (folder / "summary.json").write_text(content + "\n", encoding="utf-8")
    print(content)


if __name__ == "__main__":
    main()
