"""Extend the frozen human split with reviewed vertebrates; retain old assignments.
Arguments: new exports directory, 508-protein dataset directory, excluded 79-row TSV.
"""
import csv
import hashlib
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

out, previous, pilot = map(Path, sys.argv[1:])
rng = random.Random(20261125)

def read(path):
    with path.open(encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def plain(path):
    with path.open(encoding="utf-8-sig") as f:
        return list(csv.reader(f, delimiter="\t"))

prior = {r["accession"]:r for r in read(previous/"split.tsv")}
legacy = plain(pilot)
excluded = {a for a, _, _ in legacy}
used = {s for _, _, s in legacy + plain(previous/"sequences.tsv")}
records = {}
pool = defaultdict(lambda: [[], []])
domains = {"development":set(), "comparison":set()}
for label, group in [(1,"positive"),(0,"negative")]:
    for r in sorted(read(out/f"{group}.tsv"),key=lambda r:r["Entry"]):
        a,s,org = r["Entry"],r["Sequence"],r["Organism (ID)"]
        pf = set(r["Pfam"].strip(";").split(";"))-{""}
        if a in excluded:
            domains["development"].update(pf)
            continue
        if a in prior:
            split = prior[a]["split"]
            assert hashlib.sha256(s.encode()).hexdigest()==prior[a]["sequence_sha256"]
            records[a]=(label,s,org,pf,split,"previous")
            domains[split].update(pf)
            continue
        tm=r["Transmembrane"]
        valid="ECO:0000269" in tm if label else not tm and "membrane" not in r["Subcellular location [CC]"].lower()
        if org!="9606" and valid and pf and s not in used and re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+",s):
            pool[org][label].append((a,label,s,org,pf))
            used.add(s)
assert set(records)==set(prior)
assert domains["development"].isdisjoint(domains["comparison"])
new=[]
for org, pair in sorted(pool.items()):
    count=min(map(len,pair))
    for candidates in pair:
        new.extend(rng.sample(candidates,count))
rng.shuffle(new)
# Reserve fresh Pfam groups first, with equal positive and negative targets.
# No prior development, evaluation or pilot domains are eligible for this stratum.
seen_domains=domains["development"] | domains["comparison"]
target=max(1,len(new)//8)
counts=[0,0]
for a,label,s,org,pf in new:
    if counts[label]<target and not pf & seen_domains:
        records[a]=(label,s,org,pf,"comparison","new")
        domains["comparison"].update(pf)
        counts[label]+=1
for a,label,s,org,pf in new:
    if a not in records and not pf & domains["comparison"]:
        records[a]=(label,s,org,pf,"development","new")
        domains["development"].update(pf)
assert domains["development"].isdisjoint(domains["comparison"])
assert len({r[1] for r in records.values()})==len(records)
with (out/"sequences.tsv").open("w",newline="") as seq, (out/"split.tsv").open("w",newline="") as split, (out/"development.tsv").open("w",newline="") as dev:
    sw,mw,dw=[csv.writer(f,delimiter="\t") for f in (seq,split,dev)]
    mw.writerow(["accession","label","organism_id","split","cohort","pfam","sequence_sha256"])
    dw.writerow(["accession"])
    for a,(label,s,org,pf,part,cohort) in sorted(records.items()):
        sw.writerow([a,label,s])
        mw.writerow([a,label,org,part,cohort,";".join(sorted(pf)),hashlib.sha256(s.encode()).hexdigest()])
        if part=="development":dw.writerow([a])
from collections import Counter
print(len(records),Counter((r[4],r[5],r[0]) for r in records.values()),"species",len({r[2] for r in records.values()}))
