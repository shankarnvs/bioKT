# Rice variant annotation correctness benchmark

This benchmark compares BioKt's packaged annotator with manually labeled
synthetic VCF and GFF3 cases. The expected labels use BioKt's terminology.

From the repository root with JDK 17 on your path, build and run on Windows PowerShell:

```powershell
.\biokt\gradlew.bat -p biokt installDist
java -cp "biokt/build/install/biokt/lib/*" benchmarks/rice-annotation/RiceAnnotationBenchmark.java benchmarks/rice-annotation
```

On Linux or macOS, use `./biokt/gradlew -p biokt installDist` followed by the
same Java command. Redirect stdout to a file to save the TSV. The summary goes
to stderr. Exit code 1 means there was a mismatch, missing result, or unexpected
result. Exit code 0 means all results matched.

The baseline at commit `9e1124a` has five matches and two mismatches. The 5′ UTR
SNP receives a generic exon label, and the deletion anchored at position 114
receives an exon-indel label even though it removes the first coding base at
position 115. The expected labels record the intended behavior so both cases
serve as regression checks. With both fixes applied, all seven cases should pass.
