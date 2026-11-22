"""Export scikit-learn Iris models and prepare the fixed UniProt test panel."""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import sklearn
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB


SEEDS = (17, 23, 41, 59, 83)
PANEL_SIZE = 30
STANDARD_PROTEIN = re.compile("[ACDEFGHIKLMNPQRSTVWY]+")


def export_iris(output):
    iris = load_iris()
    labels = iris.target_names[iris.target]
    features = list(iris.feature_names)
    fields = ["model", "seed", "row", "expected", "reference", *features, *iris.target_names]
    with (output / "iris-reference.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(fields)
        for seed in SEEDS:
            train, test = train_test_split(
                range(len(labels)), test_size=0.30, random_state=seed, stratify=labels
            )
            for algorithm, model in (
                ("naive_bayes", GaussianNB()),
                ("linear", LogisticRegression(max_iter=1000)),
            ):
                model.fit(iris.data[train], labels[train])
                classes = model.classes_.tolist()
                weights = {"algorithm": algorithm, "name": f"Iris {algorithm}", "features": features}
                if algorithm == "naive_bayes":
                    weights.update(
                        priors=dict(zip(classes, model.class_prior_.tolist())),
                        means=dict(zip(classes, model.theta_.tolist())),
                        variances=dict(zip(classes, model.var_.tolist())),
                    )
                else:
                    weights.update(
                        weights=dict(zip(classes, model.coef_.tolist())),
                        biases=dict(zip(classes, model.intercept_.tolist())),
                    )
                filename = f"iris-{algorithm}-{seed}.json"
                (output / filename).write_text(json.dumps(weights, indent=2), encoding="utf-8")
                predicted = model.predict(iris.data[test])
                probabilities = model.predict_proba(iris.data[test])
                for i, row in enumerate(test):
                    writer.writerow([
                        filename, seed, row, labels[row], predicted[i],
                        *iris.data[row].tolist(), *probabilities[i].tolist(),
                    ])


def read_panel(path, label):
    with path.open(encoding="utf-8-sig", newline="") as infile:
        rows = list(csv.reader(infile, delimiter="\t"))
    assert rows[0][:2] == ["Entry", "Sequence"], path
    selected = []
    for accession, sequence, transmembrane, location in rows[1:]:
        if not STANDARD_PROTEIN.fullmatch(sequence):
            continue
        if label == "membrane" and "ECO:0000269" not in transmembrane:
            continue
        if label == "soluble" and (transmembrane or "membrane" in location.casefold()):
            continue
        selected.append((accession, label, sequence))
    selected.sort()
    if len(selected) < PANEL_SIZE:
        raise ValueError(f"{path}: need {PANEL_SIZE} eligible reviewed proteins, got {len(selected)}")
    return selected[:PANEL_SIZE]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--positive", type=Path, required=True)
    parser.add_argument("--negative", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    export_iris(args.output)
    panel = read_panel(args.positive, "membrane") + read_panel(args.negative, "soluble")
    if len({entry for entry, _, _ in panel}) != len(panel):
        raise ValueError("Overlapping UniProt accessions")
    with (args.output / "protein-panel.tsv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["accession", "expected", "sequence"])
        writer.writerows(panel)
    manifest = {
        "sklearn": sklearn.__version__,
        "iris_seeds": SEEDS,
        "panel_size": len(panel),
        "positive_sha256": hashlib.sha256(args.positive.read_bytes()).hexdigest(),
        "negative_sha256": hashlib.sha256(args.negative.read_bytes()).hexdigest(),
        "panel_sha256": hashlib.sha256((args.output / "protein-panel.tsv").read_bytes()).hexdigest(),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
