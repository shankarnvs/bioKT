"""Run from repository root after installDist. Arguments: archived ML input folder, output folder."""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from prepare import STANDARD_PROTEIN

inputs, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
with (inputs / "protein-panel.tsv").open() as stream:
    development = {r["accession"] for r in csv.DictReader(stream, delimiter="\t")}
records = []
for filename, label in [("uniprot-positive-experimental.tsv", 1), ("uniprot-negative.tsv", 0)]:
    with (inputs / filename).open(encoding="utf-8-sig") as stream:
        for r in csv.DictReader(stream, delimiter="\t"):
            tm = r["Transmembrane"]
            eligible = "ECO:0000269" in tm if label else not tm and "membrane" not in r["Subcellular location [CC]"].lower()
            if eligible and STANDARD_PROTEIN.fullmatch(r["Sequence"]):
                records.append((r["Entry"], label, r["Sequence"]))
train_seq = {s for a, _, s in records if a in development}
records = [r for r in records if r[0] in development or r[2] not in train_seq]
assert len({r[0] for r in records}) == len(records)
panel = output / "sequences.tsv"
with panel.open("w", newline="") as stream:
    csv.writer(stream, delimiter="\t").writerows(records)
run = subprocess.run(["java", "-cp", "biokt/build/install/biokt/lib/*",
    "benchmarks/ml/MembraneFeatures.java", str(panel), str(output / "regions.tsv")],
    check=True, text=True, capture_output=True)
rows = [line.split("\t") for line in run.stdout.splitlines()]
assert [r[0] for r in rows] == [a for a, _, _ in records]
x = np.array([r[2:] for r in rows], dtype=float)
y = np.array([label for _, label, _ in records])
train = np.array([a in development for a, _, _ in records])
assert train.sum() == len(development) == 60
test = ~train
print("Holdout positive/negative:", sum(y[test]), sum(1-y[test]))
def report(name, predicted):
    truth = y[test]
    counts = [sum((truth == a) & (predicted == b)) for a,b in [(1,1),(1,0),(0,1),(0,0)]]
    print(name, "TP FN FP TN:", *counts)
report("Shipped classifier", np.array([r[1]=="membrane" for r in rows])[test])
with (output / "heldout_predictions.tsv").open("w",newline="") as stream:
    writer = csv.writer(stream,delimiter="\t")
    writer.writerow(["model","accession","expected","score","threshold","prediction"])
    for name, n in [("Original six",6),("Plus peak",7),("Plus stretches and coverage",9),("All local features",x.shape[1])]:
        model = make_pipeline(StandardScaler(), LogisticRegression(C=1, max_iter=2000))
        model.fit(x[train,:n],y[train])
        # Highest training-negative score gives a common 100% training specificity target.
        threshold = model.predict_proba(x[train & (y==0),:n])[:,1].max()
        scores = model.predict_proba(x[test,:n])[:,1]
        predicted = scores > threshold
        report(name,predicted)
        for i,score,call in zip(np.flatnonzero(test),scores,predicted):
            writer.writerow([name,records[i][0],y[i],score,threshold,int(call)])
