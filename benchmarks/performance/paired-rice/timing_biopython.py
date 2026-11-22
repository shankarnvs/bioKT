import gzip, sys, time
from collections import defaultdict
from pathlib import Path
deps, root, batches = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
if deps:
    sys.path.insert(0, deps)
from Bio import SeqIO
from Bio.Seq import Seq

with gzip.open(root/'rice-chr9.fa.gz','rt') as handle:
    reference = SeqIO.read(handle, 'fasta')
genes, transcripts = {}, defaultdict(list)
coding = defaultdict(list)
with gzip.open(root/'rice-exercise-region.gff3.gz','rt') as handle:
    for row in handle:
        if not row.strip() or row.startswith('#'):
            continue
        chrom, _, kind, start, end, _, strand, phase, attributes = row.rstrip().split('\t')
        if chrom != reference.id:
            continue
        attrs = dict(part.split('=',1) for part in attributes.split(';') if '=' in part)
        if kind in ('gene','pseudogene'):
            genes[attrs['ID']] = (int(start),int(end))
        elif kind in ('mRNA','transcript','processed_transcript'):
            transcripts[attrs['Parent'].split(',')[0]].append((attrs['ID'],strand))
        elif kind == 'CDS':
            coding[attrs['Parent'].split(',')[0]].append((int(start),int(end),int(phase)))
models = []
for gene, bounds in genes.items():
    if transcripts[gene]:
        tid, strand = max(transcripts[gene],key=lambda t:sum(e-s+1 for s,e,_ in coding[t[0]]))
        segments = sorted(coding[tid],reverse=strand=='-')
        if segments:
            models.append((*bounds,tid,strand,segments))

variants = [(int(pos),ref,alt) for _,pos,_,ref,alt in
            (line.split() for line in (root/'rice_vep.vcf').read_text().splitlines() if line.strip())]
assert len(variants) == 50

def annotate():
    rows = []
    for pos, ref, alt in variants:
        for begin,end,tid,strand,segments in models:
            if not begin <= pos <= end or not any(s<=pos<=e for s,e,_ in segments):
                continue
            phase=segments[0][2]
            preceding = 0
            for s,e,_ in segments:
                if s <= pos <= e:
                    offset = preceding + (e-pos if strand=='-' else pos-s) - phase
                    break
                preceding += e-s+1
            if offset<0 or reference.seq[pos-1]!=ref:
                continue
            codon_start = offset//3*3 + phase
            window = []
            preceding = 0
            for s,e,_ in segments:
                length = e-s+1
                for i in range(max(0,codon_start-preceding), min(length,codon_start+3-preceding)):
                    window.append(e-i if strand=='-' else s+i)
                if len(window)==3:
                    break
                preceding += length
            if len(window) != 3:
                continue
            old = Seq(''.join(reference.seq[p-1] for p in window))
            base = Seq(alt)
            if strand=='-':
                old,base = old.complement(),base.complement()
            offset_in_codon = offset%3
            new = old[:offset_in_codon]+base+old[offset_in_codon+1:]
            rows.append(f'{pos} {tid} {old}>{new} {old.translate()}{offset//3+1}{new.translate()}')
    return rows

print('['+', '.join(annotate())+']')
if batches:
    count = sum(len(annotate()) for _ in range(1000))
    start = time.perf_counter()
    count += sum(len(annotate()) for _ in range(batches))
    print((time.perf_counter()-start)*1000/batches)
    assert count == batches+1000