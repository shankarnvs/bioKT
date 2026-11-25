import gzip
import sys
from pathlib import Path
work = Path(sys.argv[1])
source = Path(sys.argv[2])
sys.path.insert(0, str(work/'deps'))
from Bio import SeqIO, __version__
from Bio.Seq import Seq

rows = [line.split() for line in (work/'biokt-codons.tsv').read_text(encoding='utf-8-sig').splitlines()]
assert len(rows) == len({row[0] for row in rows}) == 64
for codon, actual in rows:
    assert str(Seq(codon).translate(table=1)) == actual, (codon, actual)
print(f'Biopython {__version__}: 64/64 BioKt standard codons matched')

for name, before, after, expected in [
    ('missense','GCT','GAT','AD'), ('synonymous','GCT','GCA','AA'),
    ('stop gain','TGG','TAG','W*'), ('stop loss','TAA','TCA','*S')
]:
    result = str(Seq(before+after).translate(table=1))
    assert result == expected, (name, result)
    print(f'{name}: {before}>{after} {result[0]}>{result[1]}')

with gzip.open(source/'rice-chr9.fa.gz','rt') as f:
    ref = SeqIO.read(f,'fasta')
transcript = 'transcript:Os09t0103500-01'
segments = []
with gzip.open(source/'rice-exercise-region.gff3.gz','rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        chrom, _, feature, start, end, _, strand, phase, attrs = line.rstrip().split('\t')
        fields = dict(part.split('=',1) for part in attrs.split(';') if '=' in part)
        if chrom == ref.id and feature == 'CDS' and transcript in fields.get('Parent','').split(','):
            assert strand == '+'
            segments.append((int(start),int(end),int(phase)))
segments.sort()
assert segments and segments[0][2] == 0
position = 437255
assert any(row.split() == [ref.id,str(position),'.','T','G'] for row in (source/'rice_vep.vcf').read_text().splitlines())
positions = [pos for start,end,_ in segments for pos in range(start,end+1)]
cds = Seq(''.join(ref.seq[pos-1] for pos in positions))
offset = positions.index(position)
assert cds[offset] == 'T'
start = offset // 3 * 3
old = cds[start:start+3]
alt = cds[:offset] + Seq('G') + cds[offset+1:]
new = alt[start:start+3]
result = (str(old),str(new),str(old.translate()),str(new.translate()),offset//3+1)
assert result == ('CTC','CGC','L','R',16), result
print(f'IRGSP-1.0 chr9:{position} T>G: {old}>{new}, Leu{result[4]}Arg')