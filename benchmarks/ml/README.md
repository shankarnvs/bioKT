# BioKt machine learning benchmark

The first check trains scikit-learn 1.7.2 Gaussian Naive Bayes and multinomial
logistic regression on Iris and exports the fitted parameters to BioKt JSON.
Each of five fixed seeds produces a stratified 105-row training split and a
45-row test split. BioKt loads each model and the report compares its labels
and class probabilities with scikit-learn on the held-out rows. Accuracy and
macro F1 compare predictions with the Iris labels. The same observations can
appear in the test sets of multiple splits.

The second check evaluates the existing membrane predictor without tuning.
Use UniProtKB release 2026_04 exports, sorted by accession, with the fields
`accession,sequence,ft_transmem,cc_subcellular_location` and 40 results per
query. The common filter is `(reviewed:true) AND (organism_id:9606) AND
(length:[80 TO 1200])`. For the positive candidates add
`AND (ft_transmem_exp:*)`. For the cytoplasm-annotated comparison candidates
add `AND (cc_scl_term:SL-0086) AND NOT (ft_transmem:*) AND NOT
(cc_scl_term:SL-0162)`. The preparation script selects 30 per class after
filtering on sequence alphabet and annotation evidence. These labels are
annotation proxies for a convenience sample of reviewed human proteins.

Run from the repository root in a Python environment with scikit-learn 1.7.2
and with JDK 17 available. Replace the two input TSV paths with the saved
UniProt exports.

```powershell
.\biokt\gradlew.bat -p biokt installDist
python benchmarks/ml/prepare.py biokt/build/ml-benchmark --positive positives.tsv --negative negatives.tsv
java -cp "biokt/build/install/biokt/lib/*" benchmarks/ml/MlBenchmark.java biokt/build/ml-benchmark
python benchmarks/ml/summarize.py biokt/build/ml-benchmark
```

Preserve the original input exports alongside the output manifest and summary
to reproduce this exact panel. This benchmark measures prediction quality and
model export fidelity. It does not measure inference latency.

### Expanded membrane model

The current `predictMembrane` uses the bundled `biokt/membrane-model.tsv`
for standard protein sequences. Inputs cover the ProtParam descriptor families,
the hydrophobic, charged and aliphatic composition measures, and the ten local
candidate-region features. Variable-length flexibility and hydropathy profiles
use their minimum, mean and maximum; profiles shorter than one window use zeros.
The original six-feature model remains accessible via `PretrainedModels`.

From the repository root, build the distribution and supply the three-column
accession/label/sequence TSV (integer labels, no header) and the original
development panel. The training command writes the weights into the repository.

```powershell
python benchmarks/ml/train_membrane.py sequences.tsv biokt/build/ml-benchmark/protein-panel.tsv biokt/build/expanded-membrane
.\biokt\gradlew.bat -p biokt test installDist
java -cp "biokt/build/install/biokt/lib/*" benchmarks/ml/ExpandedMembraneBenchmark.java sequences.tsv predict
```

The fit standardizes features using only development proteins and uses L2
logistic regression (C = 0.1). The class with the higher score wins.
The original 60/19 split was exploratory. The human-only revision used 380 new
development proteins. Another 128 proteins are reserved for evaluation.

To reconstruct the larger panel, use release 2026_04 exports from the same
queries above, adding the field `xref_pfam`. Save all positive results as
`positive.tsv` and all negative results as `all-negatives.tsv`. Run
`python benchmarks/ml/expand_panel.py OUTPUT_FOLDER ORIGINAL_79_SEQUENCE_TSV`.
Use the generated sequences.tsv and development.tsv as the training inputs.
The split manifest lists the Pfam identifiers and sequence hashes. The split
excludes proteins without Pfam annotations, rejects earlier accessions and exact
sequences, and keeps connected groups sharing a Pfam identifier together.

The current vertebrate revision retains the human split and adds 229 proteins.
Export the same two UniProt queries with `taxonomy_id:7742` in place of
`organism_id:9606`, include `organism_id` among the returned fields, and save
the complete results as positive.tsv and negative.tsv. Run
`python benchmarks/ml/extend_vertebrate_panel.py VERTEBRATE_FOLDER HUMAN_FOLDER ORIGINAL_79_SEQUENCE_TSV`.
Use the generated sequences.tsv and development.tsv with train_membrane.py.
The vertebrate split.tsv identifies earlier and new entries, organism, Pfam
identifiers and sequence checksum. Training uses 545 entries. Report the 128
previously examined human cases separately from the 64 new comparison cases.
