"""Build a balanced panel from UniProt TSVs and an earlier accession/label/sequence panel."""
import csv
import hashlib
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

folder, previous = map(Path, sys.argv[1:])
with previous.open() as stream:
    old = list(csv.reader(stream, delimiter="\t"))
old_ids = {a for a, _, _ in old}
old_sequences = {s for _, _, s in old}
old_domains = set()
groups = [[], []]
used = set(old_sequences)
for label, filename in [(1,"positive.tsv"), (0,"all-negatives.tsv")]:
    with (folder/filename).open(encoding="utf-8-sig") as stream:
        records = sorted(csv.DictReader(stream, delimiter="\t"),key=lambda r:r["Entry"])
    for r in records:
        a, seq = r["Entry"], r["Sequence"]
        domains = set(r["Pfam"].strip(";").split(";")) - {""}
        if a in old_ids or seq in old_sequences:
            old_domains.update(domains)
            continue
        tm = r["Transmembrane"]
        eligible = ("ECO:0000269" in tm if label else not tm and "membrane" not in r["Subcellular location [CC]"].lower())
        if not eligible or not domains or seq in used or not re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+",seq):
            continue
        groups[label].append((a,label,seq,domains))
        used.add(seq)
rng = random.Random(20261124)
n = min(map(len,groups))
records = sorted(rng.sample(groups[0],n) + rng.sample(groups[1],n))
# Connected components keep all shared Pfam annotations on one side of the split.
parent = {}
def root(x):
    while parent.setdefault(x,x) != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for a,_,_,domains in records:
    for domain in sorted(domains):
        parent[root(domain)] = root(a)
components = defaultdict(list)
for row in records:
    components[root(row[0])].append(row)
eligible = [g for g in components.values() if not any(r[3] & old_domains for r in g)]
target = round(n/4)
# Choose by class counts only. No predictions or fitted parameters enter this split.
test = set()
counts = [0,0]
rng.shuffle(eligible)
for component in eligible:
    proposed = [counts[i]+sum(r[1]==i for r in component) for i in (0,1)]
    if sum(abs(v-target) for v in proposed) < sum(abs(v-target) for v in counts):
        test.update(r[0] for r in component)
        counts = proposed
assert all(counts) and all(v<n for v in counts)
def domains_for(testing):
    return set().union(*(r[3] for r in records if (r[0] in test)==testing))
assert not domains_for(True) & (domains_for(False) | old_domains)
with (folder/"sequences.tsv").open("w",newline="") as seqfile, (folder/"split.tsv").open("w",newline="") as manifest, (folder/"development.tsv").open("w",newline="") as dev:
    sw,mw,dw = [csv.writer(f,delimiter="\t") for f in (seqfile,manifest,dev)]
    mw.writerow(["accession","label","split","pfam","sequence_sha256"])
    dw.writerow(["accession"])
    for a,label,seq,domains in records:
        sw.writerow([a,label,seq])
        mw.writerow([a,label,"comparison" if a in test else "development",";".join(sorted(domains)),hashlib.sha256(seq.encode()).hexdigest()])
        if a not in test: dw.writerow([a])
print("eligible negatives/positives",list(map(len,groups)),"selected per class",n,"test negatives/positives",counts)
