# BioKt JVM performance benchmarks

The JMH suite measures three operations on a fixed synthetic rice workload:

- Parse VCF text already in memory.
- Annotate parsed variants using the generic coding SNP labels.
- Annotate the same variants with a reference sequence and per-allele codon effects.

Each operation uses 128 or 2,048 variants across 128 genes on one contig.
Fixture creation, GFF3 parsing, FASTA loading, and file I/O are outside the
timed operations. Annotation uses pre-parsed VCF data. Setup checks that every
variant resolves to one coding SNP effect before timing begins.

## Run the baseline

Use JDK 17 and run these commands from the repository root.

```powershell
New-Item -ItemType Directory -Force biokt/build/reports/jmh | Out-Null
.\biokt\gradlew.bat -p biokt jmh '--args=-foe true -prof gc -rf json -rff build/reports/jmh/rice-baseline.json'
```

On Linux or macOS, use `mkdir -p biokt/build/reports/jmh` and replace
`.\biokt\gradlew.bat` with `./biokt/gradlew`.

JMH 1.37 runs one thread with two JVM forks, three 1-second warmup
iterations, and five 1-second measurement iterations per fork. The forked JVM
uses a 256 MiB initial heap and a 512 MiB maximum heap. The primary metric is
milliseconds per full operation. The GC profiler adds allocated bytes per
operation, which measures allocation volume rather than peak memory.

Keep the fixture sizes, JVM settings, and machine fixed when comparing commits.
Record the commit, JDK, operating system, CPU, and result JSON. The synthetic
workload provides a repeatable local baseline. Android runtime performance,
large real datasets, and cold file-loading costs require separate measurements.

## Human serum albumin pilot

`AlbuminBindingBenchmark` uses the mature chain of human serum albumin
(UniProt P02768, residues 25–609) and three ligands from one published
binding study. The JMH benchmark reports time and allocation per ligand
for scoring pre-parsed molecules and for parsing plus scoring. It reuses
the protein in both paths and excludes file I/O and fixture setup.

```powershell
New-Item -ItemType Directory -Force biokt/build/reports/jmh | Out-Null
.\biokt\gradlew.bat -p biokt molecularBindingCheck jmh '--args=AlbuminBindingBenchmark -foe true -prof gc -rf json -rff build/reports/jmh/albumin-baseline.json'
```

The comparison converts the reported association constants to apparent
dissociation constants in micromolar. It holds out each ligand in turn,
fits the calibration on the other two, and prints the held-out prediction,
fold error, and whether its score falls outside the training range. It also
reports the fold error of predicting the geometric mean affinity of the
training ligands. BioKt's score uses sequence and molecular descriptors.
The pilot evaluates one target and one assay, with three related ligands.

Sources: [UniProt P02768](https://www.uniprot.org/uniprot/P02768) and
[Kober and Sjöholm, 1980](https://pubmed.ncbi.nlm.nih.gov/7464807/).

## Interaction pairs pilot

The fixtures use synthetic 120-residue proteins, synthetic 20-base
DNA and RNA sequences, and pre-parsed aspirin and acetaminophen. The
benchmarks cover protein–protein scoring, DNA duplex thermodynamics,
molecule–DNA profiling, molecule–RNA scoring, and Morgan fingerprint
similarity between molecules. The fingerprint benchmarks time the
similarity calculation both with and without fingerprint generation.
The separate albumin pilot covers molecule–protein scoring.

```powershell
New-Item -ItemType Directory -Force biokt/build/reports/jmh | Out-Null
.\biokt\gradlew.bat -p biokt interactionChecks jmh '--args=InteractionBenchmarks -foe true -f 2 -wi 2 -i 3 -prof gc -rf json -rff build/reports/jmh/interaction-baseline.json'
```

The fixture checks verify protein-score symmetry, a perfect duplex,
finite molecule–sequence outputs, and fingerprint identity and symmetry.
A terminal mismatch is also reported as a diagnostic. Fingerprint
similarity measures structural overlap. Physical association between
two small molecules and experimental binding accuracy are outside
the scope of these timings.

### Measured DNA melting temperatures

Run the four unmodified, self-complementary DNA duplexes from Table 1 of
[Rozners and Moulder](https://pmc.ncbi.nlm.nih.gov/articles/PMC373285/).
The check uses the reported 2 µM total strand concentration and 300 mM
NaCl. The study buffer also contained sodium cacodylate and EDTA, which
the BioKt API cannot represent separately.

```powershell
.\biokt\gradlew.bat -p biokt installDist
java -cp "biokt/build/install/biokt/lib/*" benchmarks/performance/DuplexMeltingCheck.java
```

The check prints predicted and measured melting temperatures, individual
errors, mean absolute error, and mean signed error.
