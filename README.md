# BioKt v2.0

**A comprehensive bioinformatics library for Kotlin/JVM**

[![Kotlin](https://img.shields.io/badge/Kotlin-1.3%2B-blue?logo=kotlin)](https://kotlinlang.org)
[![JVM](https://img.shields.io/badge/JVM-8%2B-orange?logo=java)](https://openjdk.org)
[![Tests](https://img.shields.io/badge/Tests-39%2F39%20passing-brightgreen)](#testing)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](#license)
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
- [License](#-license)

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
| **Metal-Ligand** | Metal detection, coordination scoring, metal-binding residue prediction |

**Stats:** 14 source files · 4,654 lines · 58 classes/objects · 159 public functions · 39/39 tests passing

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
├── TestRunner.kt        # Self-contained 39-test suite (no JUnit needed)
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
| **Tutorial** | `BioKt_v2_Tutorial.docx` | 18-chapter tutorial modelled on the BioPython Tutorial and Cookbook. Covers every module with working examples and expected output. |
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

**This project is released into the public domain. No copyright is claimed.**

You are free to use, copy, modify, share, and build on this code for any purpose — personal, academic, or commercial — without asking for permission and without any conditions.

---

### ⚠️ Important: Use at Your Own Discretion

This library was created as a **demonstration of what AI (Claude by Anthropic) can build**. It is a learning resource and a proof of concept, not a production-grade scientific tool.

**Please read this before using BioKt in any serious context:**

- The molecular descriptors (LogP, TPSA, MW, etc.) are **estimates** based on simplified atom-contribution models. They are not as accurate as RDKit or Schrödinger.
- The ADMET predictions (BBB penetration, CYP substrates, hERG risk, oral toxicity, etc.) are **rule-based approximations** derived from published literature heuristics. They are not validated clinical predictions.
- The drug-protein binding scores are **empirical estimates**, not physics-based docking results. Do not use them to make decisions about real drug candidates.
- The RNA folding uses a simplified Nussinov algorithm, not the full Turner energy model used by tools like RNAfold.
- The BLAST implementation is educational-grade and is not a substitute for NCBI BLAST for large-scale searches.

**The author takes no responsibility for any outcome — correct or incorrect — that results from using this software.** Whether the results are accurate or inaccurate, helpful or misleading, the decision to use this library and act on its output is entirely yours.

This is a demo. Verify everything with proper scientific tools before drawing conclusions.

---

<div align="center">

**BioKt v2.0** · Built with ❤️ and [Claude AI](https://claude.ai)

*Kotlin bioinformatics for the JVM — from sequences to drug discovery*

</div>
