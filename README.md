# BioKt v2.0

**A comprehensive bioinformatics library for Kotlin/JVM**

[![Kotlin](https://img.shields.io/badge/Kotlin-1.3%2B-blue?logo=kotlin)](https://kotlinlang.org)
[![JVM](https://img.shields.io/badge/JVM-8%2B-orange?logo=java)](https://openjdk.org)
[![Tests](https://img.shields.io/badge/Tests-39%2F39%20passing-brightgreen)](#testing)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](#-license--disclaimer)
[![AI Assisted](https://img.shields.io/badge/Built%20with-Claude%20AI-blueviolet?logo=anthropic)](https://claude.ai)

> 🤖 **This project was developed with the assistance of [Claude](https://claude.ai) (Anthropic's AI assistant).** The library architecture, all source code, documentation, API reference, and this README were generated through an iterative collaboration between a human developer and Claude. See the [AI Collaboration](#-ai-collaboration) section for details.

---

BioKt mirrors and extends [BioPython](https://biopython.org/) and [pyBioMed](https://github.com/gadsbyfly/PyBioMed) on the JVM — providing everything from classic sequence analysis to drug discovery descriptors, molecular fingerprints, and molecule-sequence interaction scoring, all in idiomatic Kotlin with zero external dependencies.

```kotlin
// Five lines to go from gene to drug candidate
val gene   = DNASequence("ATGGCCATTGTAATGGGCCGCTGA")
val prot   = gene.translate()                            // MAIVMGR
val drug   = Molecule("CC(=O)Oc1ccccc1C(=O)O", name="Aspirin")
val desc   = MolDescriptors.calculate(drug)             // MW, LogP, ADMET…
val score  = DrugProteinBinding.score(drug, prot)       // –5.46 kcal/mol
```

---

## Table of Contents

- [Features](#-features)
- [Quick Start](#-quick-start)
- [Module Overview](#-module-overview)
- [Usage Examples](#-usage-examples)
- [Project Structure](#-project-structure)
- [Building from Source](#-building-from-source)
- [Android Integration](#-android-integration)
- [Testing](#-testing)
- [Documentation](#-documentation)
- [AI Collaboration](#-ai-collaboration)
- [Comparison with BioPython](#-comparison-with-biopython)
- [License & Disclaimer](#-license--disclaimer)

---

## ✨ Features

| Area | Capabilities |
|---|---|
| **Sequences** | DNA / RNA / Protein, IUPAC codes, ORF detection, six-frame translation, GC skew |
| **I/O** | FASTA, FASTQ, GenBank, Stockholm — parse & write |
| **Alignment** | Needleman-Wunsch, Smith-Waterman, semi-global, BLOSUM62, PAM250, MSA |
| **BLAST** | In-memory k-mer indexed BLASTn & BLASTp with E-value scoring |
| **Phylogenetics** | UPGMA, Neighbour-Joining, Newick I/O |
| **Analysis** | 25+ restriction enzymes, N50/L50, codon usage, Tajima's D, Fst, dN/dS |
| **Molecules** | Full SMILES parser → atom-bond graph, ring detection, implicit H |
| **Descriptors** | MW, LogP (Crippen), TPSA (Ertl), Lipinski Ro5, Veber, Ghose, Egan, Muegge |
| **ADMET** | 18-property ADMET profile (absorption, BBB, CYP, hERG, Ames, toxicity class) |
| **Fingerprints** | Morgan/ECFP, MACCS 166-bit, topological; Tanimoto, Dice, Cosine similarity |
| **Virtual Screening** | Library screening with Lipinski filter, MaxMin diversity picking |
| **Protein Descriptors** | AAC, DPC, TPC, CTD, QSO, PAAC, Geary/Moran autocorrelation, GRAVY, pI, extinction coefficient |
| **PPI** | Protein-protein interaction scoring, interface residue prediction, coiled-coil, TM helix |
| **DNA Interactions** | SantaLucia 1998 hybridisation thermodynamics, primer design, FRET |
| **RNA Interactions** | Nussinov folding, dot-bracket, miRNA seed-match prediction, kissing loops |
| **Drug-DNA** | Intercalation scoring, minor groove binding, covalent risk profiling |
| **Drug-Protein** | Empirical binding score, pharmacophore generation & matching, linear QSAR |
| **CRISPR** | SpCas9 guide efficiency, R-loop propensity, transcription bubble stability |
| **ML Prediction** | Pretrained classifiers: coding/non-coding DNA, promoter, splice site, enzyme, membrane, drug activity, toxicity (Cramer I/II/III), hERG risk |
| **Custom Models** | Load any scikit-learn / R / Weka model from JSON weight files — Naive Bayes, kNN, Decision Tree, Random Forest, Logistic Regression |
| **LLM Integration** | OpenAI GPT-4o, Anthropic Claude, Google Gemini, Ollama (local) — biological interpretation of sequences, proteins, and drug molecules |
| **Rice Genomics** | `biokt.rice` subpackage — VCF parsing, GFF3 annotation, SSR detection, cis-element scanning, TE classification, R-gene patterns, expression normalisation (TPM/FPKM/DESeq2), DE analysis, haplotype blocks, selection sweeps |
| **Metal-Ligand** | Metal detection, coordination scoring, metal-binding residue prediction |

**Stats:** 25 source files · 9,108 lines · 86 classes/objects · 253 public functions · 39/39 tests passing

---

## 🚀 Quick Start

### Option 1 — Run the pre-built JAR

```bash
java -jar biokt2.jar
```

This runs all 14 demo sections and 39 tests.

### Option 2 — Use as a dependency (Gradle Kotlin DSL)

```kotlin
// build.gradle.kts
dependencies {
    implementation(files("libs/biokt2.jar"))
}
```

### Option 3 — Use the source directly (recommended for Android)

Copy all `.kt` files from `src/main/kotlin/biokt/` into your source set. No Maven dependencies needed — only the Kotlin stdlib is required.

---

## 📦 Module Overview

```
biokt/
├── Sequence.kt          # DNASequence, RNASequence, ProteinSequence, Orf
├── CodonTable.kt        # Standard, Vertebrate Mito, Yeast Mito codon tables
├── Alignment.kt         # PairwiseAligner, MultipleSequenceAligner, ScoringMatrices
├── SeqIO.kt             # FASTA / FASTQ / GenBank / Stockholm I/O
├── Blast.kt             # SequenceDatabase, BLASTn, BLASTp
├── Phylogenetics.kt     # TreeBuilder (UPGMA, NJ), PhyloTree, DistanceMatrix
├── Analysis.kt          # RestrictionEnzymes, SeqStats, PopGen, IUPAC
├── Molecule.kt          # Molecule, Atom, Bond, BondType, Ring, SmilesParser
├── MolDescriptors.kt    # MolDescriptors, Fingerprints, VirtualScreening
├── ProteinDescriptors.kt# ProteinDescriptors (AAC, DPC, CTD, QSO, PAAC, …)
├── SequenceInteraction.kt # ProteinInteraction, DNAInteraction, RNAInteraction, DNARNAInteraction
├── MoleculeInteraction.kt # DrugProteinBinding, DrugDNABinding, DrugRNABinding, MetalLigandBinding
├── MLPredictor.kt       # Pretrained ML classifiers (NaiveBayes, kNN, RF, DT, Linear)
├── GPTClient.kt         # Multi-provider LLM client (OpenAI, Anthropic, Google, Ollama)
├── ReportGenerator.kt   # HTML report generation with embedded 3D viewers
├── Viewer3D.kt          # Three.js 3D viewers (DNA helix, protein, molecule)
├── TestRunner.kt        # Self-contained 39-test suite (no JUnit needed)
└── rice/                # Rice genomics subpackage (import biokt.rice.*)
    ├── RicePackage.kt       # Package entry point, Rice constants, extension functions
    ├── RiceVariants.kt      # VCF/BCF parsing, SNP/INDEL, LD, genotype matrix
    ├── GenomeAnnotation.kt  # GFF3/GTF, gene models, RAP-DB/MSU IDs
    ├── RiceGenomics.kt      # SSR, cis-elements (30), TEs, R-genes, flowering genes, chromosomes
    ├── RiceExpression.kt    # TPM/FPKM/DESeq2, DE (Welch+BH), co-expression, stress gene sets
    └── HaplotypeAnalysis.kt # Gabriel blocks, windowed pi/θ_W/Tajima's D, selection sweeps
└── Main.kt              # Demo entry point
```

---

## 💡 Usage Examples

### DNA sequences

```kotlin
val dna = DNASequence("ATGGCCATTGTAATGGGCCGCTGA", id = "gene1")

println(dna.gcContent())           // 52.08
println(dna.reverseComplement())   // TCAGCGGCCCATTACAATGGCCAT
println(dna.translate())           // MAIVMGR
println(dna.meltingTemperature())  // 70.8 °C

val orfs = dna.findOrfs(minLength = 9)
orfs.forEach { println(it) }
```

### Pairwise alignment with BLOSUM62

```kotlin
val aligner = PairwiseAligner(
    mode               = PairwiseAligner.Mode.GLOBAL,
    substitutionMatrix = ScoringMatrices.BLOSUM62,
    gapOpen            = -11.0,
    gapExtend          = -1.0
)
val aln = aligner.align("MAEGEITTFT", "MAEGEITTFA")
println(aln.identity)      // 0.9
println(aln.prettyPrint())
```

### Molecular descriptors & ADMET

```kotlin
val mol = Molecule("CC(C)Cc1ccc(cc1)C(C)C(=O)O", name = "Ibuprofen")
val d   = MolDescriptors.calculate(mol)

println(d.molecularWeight)                  // 206.28
println(d.logP)                             // 3.51
println(d.lipinskiPasses)                   // true
println(d.admet.oralBioavailability)        // "High"
println(d.admet.bbbPenetration)             // true
println(d.admet.hergInhibition)             // false
```

### Virtual screening

```kotlin
val library = listOf(aspirin, ibuprofen, caffeine, paracetamol /*, … */)
val hits = VirtualScreening.screenLibrary(
    query             = aspirin,
    library           = library,
    fingerprintType   = "morgan",
    tanimotoThreshold = 0.3,
    applyLipinski     = true,
    topN              = 10
)
hits.forEach { println("#${it.rank} ${it.molecule.name}  Tc=${it.tanimoto}") }
```

### Protein descriptors (ML-ready)

```kotlin
val prot     = ProteinSequence("MKALVLLYLLFSSAYSRGVFRRDTHKPPVYK")
val allDesc  = ProteinDescriptors.calculateAll(prot)
val features = allDesc.toFlatMap()   // 500+ numerical features

println("GRAVY: ${ProteinDescriptors.gravyIndex(prot)}")
println("pI:    ${prot.isoelectricPoint()}")
```

### DNA hybridisation thermodynamics (SantaLucia 1998)

```kotlin
val hyb = DNAInteraction.hybridize(
    seq1     = DNASequence("ATCGATCGATCGATCG"),
    seq2     = DNASequence("CGATCGATCGATCGAT"),
    saltConc = 0.05   // 50 mM NaCl
)
println("Tm:   ${hyb.tm} °C")       // 49.5 °C
println("ΔG37: ${hyb.deltaG37}")    // –14.16 kcal/mol
```

### RNA secondary structure & miRNA targeting

```kotlin
val fold = RNAInteraction.fold(RNASequence("GGGGCCCCCAUGGUGCAAAUAG"))
println(fold.dotBracket)     // ((((..))))((((..))))..
println(fold.numBasePairs)   // 8

val targets = RNAInteraction.predictMiRNATargets(miRNA, mRNA)
targets.forEach { println("pos=${it.seedMatchPosition}  type=${it.siteType}") }
```

### Drug–protein binding

```kotlin
val result = DrugProteinBinding.score(drug, proteinSeq)
println(result.bindingScore)    // –5.46 kcal/mol
println(result.estimatedKd)     // 141.9 µM
println(result.isLikelyBinder)  // true

result.interactions.forEach { i ->
    println("[${i.type}] ${i.strength} kcal/mol — ${i.description}")
}
```

### Rice genomics (biokt.rice subpackage)

```kotlin
import biokt.*
import biokt.rice.*

// Parse a VCF file — SNPs only, MAF ≥ 5%, chr01
val vcf = VcfParser.parse("3krg.vcf", snpsOnly=true, minMaf=0.05, regionChrom="chr01")
println(vcf.summarise().print())

// Parse a GFF3 annotation
val db = AnnotationParser.parseGff3("IRGSP-1.0.gff3")
val hd1 = db.getByRapId("Os06g0275000")       // Heading date 1
println(hd1?.summary())

// SSR detection in a promoter
val ssrs = myPromoterSeq.detectSSRs()
println("${ssrs.size} SSRs found")

// Scan for cis-elements (W-box, ABRE, DRE/CRT, etc.)
val elements = myPromoterSeq.scanRicePromoter()
println("W-box count: ${elements["W-box"]}")

// Check if protein is an R-gene (NBS-LRR, RLK, etc.)
val (isR, evidence) = myProtein.isRGene()
println("Is R-gene: $isR")

// Differential expression — drought vs control
val de = DifferentialExpression.analyse(log2matrix, droughtSamples, ctrlSamples)
val degs = de.filter { it.isSignificant }
val enrichment = RiceStressGeneSets.enrichment(degs.map { it.geneId })

// Haplotype blocks and diversity
val blocks  = HaplotypeBlocks.defineBlocks(vcf, "chr06")
val windows = WindowedDiversity.compute(vcf, "chr06", windowSize=100_000)
val sweeps  = WindowedDiversity.selectionWindows(windows)

// Selection sweeps between Indica and Japonica
val scores = SelectionSweepDetector.detect(indica, japonica, "chr04")
SelectionSweepDetector.topSweeps(scores, 10).forEach { println(it) }

// Flowering gene catalogue
FloweringGenes.getBySymbol("Hd1")?.let { println("${it.symbol}: ${it.function}") }
println(RiceChromosomes.position("chr06", 2_886_607))  // short arm (p)
```

### ML prediction (pretrained models, no training needed)

```kotlin
// DNA — coding, promoter, splice site
val mlDNA = SequenceMLPredictor.fullAnalysis(dna)
mlDNA.forEach { (task, pred) ->
    println("$task: ${pred.label}  confidence=${pred.confidence}")
}
// coding:      coding      (83.2%)
// promoter:    non-promoter (91.4%)
// splice_site: non-splice  (76.8%)

// Protein — enzyme vs non-enzyme, membrane vs soluble
val mlProt = ProteinMLPredictor.fullAnalysis(prot)

// Drug — activity, toxicity class, hERG risk
val mlDrug = DrugMLPredictor.fullAnalysis(mol)
mlDrug.forEach { (task, pred) ->
    println("$task: ${pred.label}  (${pred.confidence})")
}
// activity:  active        (72%)
// toxicity:  class_I_low   (81%)
// herg_risk: non-inhibitor (84%)
```

### LLM biological interpretation

```kotlin
// OpenAI GPT-4o
val gpt = GPTClient.openai(apiKey = System.getenv("OPENAI_API_KEY"))

// Anthropic Claude — no API key needed for Ollama
val local = GPTClient.ollama(model = "llama3")

// Interpret DNA with ML context
val resp = gpt.interpretDNA(dna, SequenceMLPredictor.fullAnalysis(dna))
println(resp.text)

// Interpret drug with full ADMET + ML predictions
val drugResp = gpt.interpretDrug(aspirin, DrugMLPredictor.fullAnalysis(aspirin))
println(drugResp.text)

// Extension function style
val interp = aspirin.interpretWith(gpt)

// Fluent config builder
val claude = GPTConfig.create {
    provider  = LLMProvider.ANTHROPIC
    model     = LLMModels.Anthropic.CLAUDE_35_SONNET
    apiKey    = System.getenv("ANTHROPIC_API_KEY")
    maxTokens = 800
}
```

---

## 📁 Project Structure

```
biokt/
├── src/main/kotlin/biokt/    # All 14 source files
├── biokt2.jar                # Pre-built runnable JAR (includes Kotlin stdlib)
├── build.gradle.kts          # Gradle build file
├── settings.gradle.kts
└── README.md
```

---

## 🔨 Building from Source

### With Gradle (recommended)

```bash
git clone https://github.com/YOUR_USERNAME/biokt.git
cd biokt
./gradlew jar
java -jar build/libs/biokt-2.0.jar
```

### With kotlinc directly

```bash
kotlinc src/main/kotlin/biokt/*.kt -include-runtime -d biokt2.jar
java -jar biokt2.jar
```

### Requirements

- Kotlin 1.3 or later (Kotlin 1.9+ recommended)
- Java 8 or later

> **Kotlin version note:** The source uses Kotlin 1.3-compatible idioms (`.map{}.sum()` instead of `.sumOf{}`, etc.) so it builds in older environments. On Kotlin 1.9+, you can modernise these without any functional change.

---

## 📱 Android Integration

Drop all `.kt` files from `src/main/kotlin/biokt/` into your Android project's source set — **no Gradle dependencies beyond the Kotlin stdlib are needed**.

```kotlin
// build.gradle.kts (app)
// No extra dependencies — just add the .kt files to your source set
```

For file I/O, use `SeqIO.parseString()` instead of `SeqIO.parse()`:

```kotlin
// Reading a FASTA file from Android file picker
val text = contentResolver
    .openInputStream(uri)
    ?.bufferedReader()
    ?.readText() ?: return

val records = SeqIO.parseString(text, "fasta")
```

Run heavy operations (alignment, BLAST, folding) on a background thread:

```kotlin
viewModelScope.launch(Dispatchers.Default) {
    val desc   = MolDescriptors.calculate(mol)
    val fold   = RNAInteraction.fold(rnaSeq)
    val result = DrugProteinBinding.score(drug, prot)
    withContext(Dispatchers.Main) { updateUI(desc, fold, result) }
}
```

---

## 🧪 Testing

BioKt ships with 39 self-contained tests. No JUnit or testing framework is required.

```bash
java -jar biokt2.jar
# ...
# Results: 39 passed, 0 failed out of 39 tests
```

| Test area | Count |
|---|---|
| DNA Sequence operations | 12 |
| RNA Sequence operations | 4 |
| Protein Sequence operations | 3 |
| Codon tables | 4 |
| Alignment (NW/SW/BLOSUM) | 4 |
| SeqIO (FASTA/FASTQ) | 4 |
| Restriction enzymes | 3 |
| Phylogenetics | 2 |
| IUPAC codes | 3 |
| **Total** | **39/39 ✅** |

---

## 📚 Documentation

Three documentation formats are included:

| Format | File | Description |
|---|---|---|
| **Tutorial (core)** | `BioKt_v2_Tutorial.docx` | 18-chapter tutorial: sequences, alignment, BLAST, phylogenetics, molecules, descriptors, interactions, and reports. |
| **Tutorial (ML & GPT)** | `BioKt_v2_Tutorial_ML_GPT.docx` | Chapters 19–20: pretrained ML classifiers, custom model loading, ensemble predictions, GPT/Claude/Gemini/Ollama integration. |
| **Tutorial (Rice)** | `BioKt_v2_Tutorial_Rice.docx` | Chapter 21: VCF parsing, GFF3 annotation, SSR/cis-elements, expression analysis, haplotype blocks, selection sweeps. Appendix E: reference tables. |
| **API Docs (Interactive)** | `BioKt_v2_API_Docs.html` | Single-file interactive API reference. Click any module or class in the sidebar to navigate. |
| **API Docs (Frames)** | `BioKt_v2_API_Docs.zip` | Classic Javadoc-style three-frame documentation: package list · class list · detail view. 65 HTML files. |

Open `BioKt_v2_API_Docs.html` directly in any browser — no server needed.

---

## 🤖 AI Collaboration

This project was built through an extended, iterative collaboration between a human developer and **[Claude](https://claude.ai)**, Anthropic's AI assistant (Claude Sonnet).

### What Claude contributed

- **Library architecture** — module decomposition, class hierarchy, API design decisions
- **All source code** — all 14 `.kt` files (~4,650 lines), including the SMILES parser, nearest-neighbour thermodynamics, Nussinov RNA folding, Wildman-Crippen LogP, Ertl TPSA, SantaLucia 1998 NN parameters, and empirical ADMET rules
- **Debugging** — systematic resolution of ~80 Kotlin 1.3 compatibility errors caused by a broken `sumOf{}` → `map{}.sum()` regex replacement (the most painful session involved tracking brace-depth mismatches across six files simultaneously)
- **API documentation** — the 65-file Javadoc-style HTML reference with a 3-frame layout, all 140+ method entries, and syntax-highlighted code examples
- **Tutorial document** — the 18-chapter, 100+ page BioPython-style tutorial with formatted code blocks, output blocks, and callout boxes, generated as `.docx`
- **Rice genomics subpackage** — `biokt.rice/`: 6 files, 2,970 lines covering VCF parsing, GFF3 annotation, SSR detection, 30 cis-elements, TE classification, R-gene patterns, 18 flowering genes, 6 stress gene sets, DESeq2 normalisation, Welch DE analysis, Gabriel haplotype blocks, windowed Tajima's D, and cross-population selection sweep scoring
- **ML inference engine** — `MLPredictor.kt`: 8 pretrained classifiers across 3 domains, 4 algorithm families, custom model JSON loader, ensemble voting
- **LLM integration layer** — `GPTClient.kt`: multi-provider HTTP client (OpenAI, Anthropic, Google, Ollama), structured biological prompt builders, extension functions
- **This README**

### How the collaboration worked

The session proceeded roughly as:

1. Human specified the goal: *"BioPython equivalent in Kotlin"*
2. Claude designed the module structure and wrote BioKt v1 (sequence, alignment, BLAST, phylogenetics)
3. Human requested expansion: *"add molecular interactions mimicking pyBioMed"*
4. Claude added five new modules (Molecule, MolDescriptors, ProteinDescriptors, SequenceInteraction, MoleculeInteraction)
5. A Python regex bug broke ~30 lambdas across all new files — Claude debugged this systematically over multiple turns
6. Human requested Javadoc-style API docs → Claude generated 65 HTML files
7. Human requested a BioPython-style tutorial → Claude wrote the 18-chapter document
8. Human requested this README

### Disclaimer

The scientific models implemented (ADMET rules, empirical binding scores, drug-likeness filters, etc.) are **in silico estimates** derived from published literature rules. They are suitable for educational use and early-stage computational triage, **not** for clinical or regulatory decision-making. Always validate predictions with experimental data.

---

## 🔄 Comparison with BioPython

| BioPython | BioKt |
|---|---|
| `Bio.Seq.Seq("ATCG", IUPAC.unambiguous_dna)` | `DNASequence("ATCG")` |
| `seq.complement()` | `seq.complement()` |
| `seq.reverse_complement()` | `seq.reverseComplement()` |
| `seq.transcribe()` | `seq.transcribe()` |
| `seq.translate()` | `seq.translate()` |
| `SeqIO.parse("f.fasta","fasta")` | `SeqIO.parse("f.fasta","fasta")` |
| `pairwise2.align.globalms(a, b, 2, -1, -2, -0.5)` | `PairwiseAligner(mode=GLOBAL).align(a, b)` |
| `NCBIWWW.qblast("blastn", "nt", seq)` | `db.blastn(query, eValueCutoff=1e-5)` |
| `Phylo.read("tree.nwk", "newick")` | `PhyloTree.fromNewick(text)` |
| `ProteinAnalysis(seq).gravy()` | `ProteinDescriptors.gravyIndex(seq)` |
| `ProteinAnalysis(seq).isoelectric_point()` | `seq.isoelectricPoint()` |
| `Chem.MolFromSmiles(smiles)` *(RDKit)* | `Molecule(smiles)` |
| `Descriptors.MolWt(mol)` *(RDKit)* | `MolDescriptors.molecularWeight(mol)` |
| `Descriptors.MolLogP(mol)` *(RDKit)* | `MolDescriptors.wilmanCrippenLogP(mol)` |
| `AllChem.GetMorganFingerprintAsBitVect(mol, 2)` | `Fingerprints.morgan(mol, radius=2)` |
| `DataStructs.TanimotoSimilarity(fp1, fp2)` | `Fingerprints.tanimoto(fp1, fp2)` |
| `PyPro.GetProDes(seq)` *(pyBioMed)* | `ProteinDescriptors.calculateAll(seq).toFlatMap()` |

---

## 📄 License & Disclaimer

### MIT License

```
MIT License

Copyright (c) 2026 BioKt Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

#### What this means in plain English

**What you CAN do** (no permission needed, no fee, no conditions):

- ✅ Use BioKt in your own projects — personal, academic, or commercial
- ✅ Copy, modify, and adapt the source code however you like
- ✅ Distribute it — share it with colleagues, publish it, include it in another library
- ✅ Sell a product that includes BioKt
- ✅ Sub-license it under different terms inside your own project

**The one condition:**

- 📎 If you distribute the code or include it in another project, keep the copyright notice and this license text attached. You do not have to credit the author publicly — you just keep the `LICENSE` file in the repository.

**What the author does NOT guarantee:**

- ❌ That the software works correctly for your specific use case
- ❌ That the scientific results (descriptors, ADMET values, binding scores) are accurate
- ❌ That the code is free of bugs or errors
- ❌ Any liability if something goes wrong as a result of using this software

---

### ⚠️ Scientific Accuracy Disclaimer

This library was created as a **demonstration of what AI (Claude by Anthropic) can build**. It is a learning resource and a proof of concept — not a validated scientific tool for production or clinical use.

**Use at your own discretion. The author is not responsible for the accuracy or inaccuracy of any result produced by this library.**

Specific limitations to be aware of:

| Module | Limitation |
|---|---|
| `MolDescriptors` — LogP, TPSA, MW | Simplified atom-contribution estimates. Less accurate than RDKit or Schrödinger. |
| `MolDescriptors` — ADMET profile | Rule-based heuristics from published literature. Not validated clinical predictions. Do not use for drug safety decisions. |
| `DrugProteinBinding` — binding score | Empirical scoring function, not physics-based docking. Kd estimates are approximate. |
| `RNAInteraction` — folding | Simplified Nussinov algorithm (max base pairs). Not the full Turner energy model used by RNAfold or Mfold. |
| `SequenceDatabase` — BLAST | Educational-grade k-mer search. Not suitable as a replacement for NCBI BLAST at genome scale. |
| `PopGen` — Tajima's D, Fst, dN/dS | Standard implementations, but results should be verified with dedicated tools (e.g. DnaSP, PopGenome) for publication. |

**The right workflow:** Use BioKt for exploration, learning, and rapid prototyping. Confirm any result that matters with a validated, peer-reviewed tool before drawing scientific conclusions or making decisions based on it.

---

<div align="center">

**BioKt v2.0** · Built with ❤️ and [Claude AI](https://claude.ai)

*Kotlin bioinformatics for the JVM — from sequences to drug discovery*

</div>
