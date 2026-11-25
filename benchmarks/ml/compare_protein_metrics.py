"""Run from repository root with JDK 17, installDist, and Biopython 1.88."""
import csv
import gzip
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData as data

panel, output = map(Path, sys.argv[1:])
proteins = {r[0]: ProteinAnalysis(r[2]) for r in csv.reader(panel.open(), delimiter="\t")}
scales = list(data.gravy_scales)
args = ["java", "-cp", "biokt/build/install/biokt/lib/*", "benchmarks/ml/ProteinMetricValues.java", str(panel), ",".join(scales)]
run = subprocess.run(args, capture_output=True, text=True, check=True)
summary = Counter()
with output.open("w", newline="") as stream, gzip.open(output.with_suffix(".values.tsv.gz"), "wt", newline="") as details:
    values = csv.writer(details, delimiter="\t")
    values.writerow(["accession","metric","component_index","biokt","biopython"])
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(["accession","metric","values","matched","mismatched","max_abs_error"])
    for line in run.stdout.splitlines():
        accession, metric, raw = line.split("\t")
        p = proteins[accession]
        counts = p.count_amino_acids()
        frac = lambda aa: sum(counts[c] for c in aa)/p.length
        if metric.startswith("gravy_"):
            expected = [p.gravy(metric[6:])]
        elif metric.startswith("profile_"):
            _, w, e = metric.split("_")
            expected = p.protein_scale(data.kd,int(w),float(e))
        else:
            functions = {
                "counts": lambda: [counts[c] for c in sorted(counts)],
                "percentages": lambda: [p.amino_acids_percent[c] for c in sorted(counts)],
                "mass": lambda: [p.molecular_weight()],
                "mono_mass": lambda: [ProteinAnalysis(p.sequence, monoisotopic=True).molecular_weight()],
                "aromaticity": lambda: [p.aromaticity()],
                "instability": lambda: [p.instability_index()],
                "flexibility": p.flexibility,
                "charge": lambda: [p.charge_at_pH(ph) for ph in (4,7,10)],
                "pI": lambda: [p.isoelectric_point()],
                "secondary": p.secondary_structure_fraction,
                "extinction": p.molar_extinction_coefficient,
                "legacy_six": lambda: [p.gravy(),frac("AVILMFWP"),frac("RKHDE"),100*(frac("A")+2.9*frac("V")+3.9*frac("IL")),p.aromaticity(),p.isoelectric_point()]
            }
            expected = functions[metric]()
        actual = [float(v) for v in raw.split(",")] if raw else []
        assert len(actual)==len(expected)
        values.writerows((accession, metric, i, a, b) for i,(a,b) in enumerate(zip(actual,expected)))
        errors = [abs(a-b) for a,b in zip(actual,expected)]
        misses = sum(e>1e-8 for e in errors)
        writer.writerow([accession,metric,len(errors),len(errors)-misses,misses,max(errors,default=0)])
        if metric=="legacy_six":
            for i,e in enumerate(errors): summary["legacy_misses_"+str(i)] += int(e>1e-8)
        else:
            summary["comparisons"] += len(errors)
            summary["mismatches"] += misses
            if misses: summary["miss_"+metric] += misses
print(len(proteins),"proteins;",len(scales),"GRAVY scales;",dict(summary))
