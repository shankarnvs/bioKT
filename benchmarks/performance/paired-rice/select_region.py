import gzip, pathlib, sys
work = pathlib.Path(sys.argv[1])
positions = [int(line.split()[1]) for line in (work/'rice_vep.vcf').read_text().splitlines() if line.strip() and not line.startswith('#')]
with gzip.open(work/'rice-release62.gff3.gz', 'rt') as source:
    rows = [l for l in source if l.startswith('9\t')]
start=min(positions)-2000
end=max(positions)+2000
bounds=[l.split('\t') for l in rows if l.split('\t')[2] in {'gene','pseudogene'} and int(l.split('\t')[3]) <= end and int(l.split('\t')[4]) >= start]
start=min([start]+[int(r[3]) for r in bounds])
end=max([end]+[int(r[4]) for r in bounds])
with gzip.open(work/'rice-exercise-region.gff3.gz','wt') as out:
    out.writelines(l for l in rows if int(l.split('\t')[3]) <= end and int(l.split('\t')[4]) >= start)
print(f'Chromosome 9 window {start}-{end}, {len(positions)} exercise variants')