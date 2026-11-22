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
