"""Train the bundled expanded predictor. Run from repository root after installDist.
Arguments: sequence TSV (no header), development accession panel, output folder.
The packaged weights are biokt/membrane-model.tsv. Rebuild installDist after training.
"""
import csv
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

panel, development, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
command = ["java", "-cp", "biokt/build/install/biokt/lib/*",
           "benchmarks/ml/ExpandedMembraneBenchmark.java", str(panel)]
features = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="\t"))
names = list(features[0])[1:]
ids = [r["accession"] for r in features]
with panel.open() as stream:
    records = {a: (int(label), seq) for a, label, seq in csv.reader(stream, delimiter="\t")}
with development.open() as stream:
    dev = {r["accession"] for r in csv.DictReader(stream, delimiter="\t")}
assert dev and dev < set(ids) and len(ids) == len(set(ids))
assert not ({records[a][1] for a in dev} & {records[a][1] for a in set(ids)-dev})
x = np.array([[float(r[n]) for n in names] for r in features])
y = np.array([records[a][0] for a in ids])
train = np.array([a in dev for a in ids])
six = ["gravy_KyteDoolitle","hydrophobic_ratio","charged_fraction","aliphatic_index","aromatic_fraction","pi"]
with (output/"predictions.tsv").open("w",newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(["model","accession","split","expected","membrane_score","call"])
    for model_name, selected in [("six_refit",six),("expanded",names)]:
        data = x[:,[names.index(name) for name in selected]]
        scaler = StandardScaler().fit(data[train])
        # Fixed L2 regularization and 0.5 cutoff, shared by both models.
        model = LogisticRegression(C=0.1,max_iter=2000).fit(scaler.transform(data[train]),y[train])
        scores = model.predict_proba(scaler.transform(data))[:,1]
        calls = scores > 0.5
        for i, accession in enumerate(ids):
            writer.writerow([model_name,accession,"development" if train[i] else "comparison",
                             y[i],scores[i],int(calls[i])])
        print(model_name,len(selected),[(int(((calls==pred)&(y==truth)&~train).sum()))
              for truth,pred in [(1,1),(1,0),(0,1),(0,0)]])
        if model_name == "expanded":
            weights = model.coef_[0]/scaler.scale_
            bias = model.intercept_[0]-np.dot(weights,scaler.mean_)
            with Path("biokt/membrane-model.tsv").open("w",newline="") as file:
                table = csv.writer(file,delimiter="\t",lineterminator="\n")
                table.writerow(["intercept",repr(float(bias))])
                table.writerows((n,repr(float(w))) for n,w in zip(selected,weights))
