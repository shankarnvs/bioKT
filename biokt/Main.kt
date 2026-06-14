package biokt

fun main() {
    println("=" .repeat(60))
    println("  BioKt — Bioinformatics Library for Kotlin")
    println("=" .repeat(60))
    println()

    // ─── 1. DNA Sequence Operations ─────────────────────────
    println("▶ 1. DNA SEQUENCE OPERATIONS")
    println("-".repeat(40))

    val dna = DNASequence("ATGCGATACGCTTACGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCATGCGATACGCTTACGATCGATCG", id = "gene1", description = "Example gene")
    println("Sequence   : ${dna.sequence}")
    println("Length     : ${dna.length} bp")
    println("GC Content : ${"%.1f".format(dna.gcContent())}%")
    println("Complement : ${dna.complement()}")
    println("Rev Comp   : ${dna.reverseComplement()}")
    println("Tm (°C)    : ${"%.1f".format(dna.meltingTemperature())}")
    println("A: ${dna.count('A')}  T: ${dna.count('T')}  G: ${dna.count('G')}  C: ${dna.count('C')}")
    println()

    // ─── 2. Transcription & Translation ─────────────────────
    println("▶ 2. TRANSCRIPTION & TRANSLATION")
    println("-".repeat(40))

    val testDna = DNASequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    val rna = testDna.transcribe()
    println("DNA    : ${testDna.sequence}")
    println("mRNA   : $rna")

    val protein = testDna.translate()
    println("Protein: $protein")

    println("\nSix-frame translation:")
    testDna.sixFrameTranslation().forEach { (frame, prot) ->
        println("  Frame $frame: ${prot.sequence}")
    }
    println()

    // ─── 3. ORF Finding ─────────────────────────────────────
    println("▶ 3. ORF FINDING")
    println("-".repeat(40))

    val genomicSeq = DNASequence("NNNATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAGNNNATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGTAA")
    val orfs = genomicSeq.findOrfs(minLength = 30)
    println("Found ${orfs.size} ORF(s) in sequence:")
    orfs.forEach { println("  $it") }
    println()

    // ─── 4. Protein Analysis ────────────────────────────────
    println("▶ 4. PROTEIN ANALYSIS")
    println("-".repeat(40))

    val prot = ProteinSequence("MKALVLLYLLFSSAYSRGVFRRDTHKPPVYKKVEHPLSPSTVQQLLQHLESLRLFHTKFSSAQALKNLQKLQESTFDDIPNLSSIPNLTQESKSASNEFSFLNKPQIQPSEQEVYVDTFHRQLLKFIPGLEEDFQDVPHFLQNLQYVFTKAQIQLMEMRQAYKDLAQEFEQKLLSEWRQKAAQQSKPAAQQQEAQAGQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPPQQEAQKPP", id = "prot1")
    println("Length        : ${prot.length} aa")
    println("MW            : ${"%.1f".format(prot.molecularWeight())} Da")
    println("pI            : ${"%.2f".format(prot.isoelectricPoint())}")
    println("Aromaticity   : ${"%.3f".format(prot.aromaticity())}")
    println("Instability   : ${"%.2f".format(prot.instabilityIndex())} (>40 = unstable)")
    println()

    // ─── 5. Sequence Alignment ──────────────────────────────
    println("▶ 5. PAIRWISE ALIGNMENT")
    println("-".repeat(40))

    val seqA = "ATCGATCGATCGATCGATCG"
    val seqB = "ATCGATCGATCGATCG"

    val globalAligner = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL)
    val globalAln = globalAligner.align(seqA, seqB)
    println("Global Alignment (Needleman-Wunsch):")
    println("  Score   : ${"%.1f".format(globalAln.score)}")
    println("  Identity: ${"%.1f".format(globalAln.identity * 100)}%")
    println("  Aligned A: ${globalAln.alignedA}")
    println("  Aligned B: ${globalAln.alignedB}")
    println()

    val localAligner = PairwiseAligner(mode = PairwiseAligner.Mode.LOCAL)
    val localAln = localAligner.align("XXXATCGATCGXXX", "ATCGATCGATCG")
    println("Local Alignment (Smith-Waterman):")
    println("  Score   : ${"%.1f".format(localAln.score)}")
    println("  Identity: ${"%.1f".format(localAln.identity * 100)}%")
    println("  Aligned A: ${localAln.alignedA}")
    println("  Aligned B: ${localAln.alignedB}")
    println()

    // Protein alignment with BLOSUM62
    val protAligner = PairwiseAligner(
        mode = PairwiseAligner.Mode.GLOBAL,
        substitutionMatrix = ScoringMatrices.BLOSUM62
    )
    val protAln = protAligner.align("MKALVLLY", "MAALVLLY")
    println("Protein Alignment (BLOSUM62):")
    println("  Score   : ${"%.1f".format(protAln.score)}")
    println(protAln.prettyPrint())

    // ─── 6. FASTA I/O ────────────────────────────────────────
    println("▶ 6. FASTA I/O")
    println("-".repeat(40))

    val fastaInput = """
        >seq1 First sequence
        ATCGATCGATCGATCGATCG
        ATCGATCG
        >seq2 Second sequence
        GCTAGCTAGCTAGCTAGCTA
        >seq3 Third sequence
        TTTTAAAACCCCGGGG
    """.trimIndent()

    val records = SeqIO.parseString(fastaInput, "fasta")
    println("Parsed ${records.size} FASTA records:")
    records.forEach { rec ->
        println("  ${rec.id}: ${rec.length} bp — ${rec.description}")
        println("    GC: ${"%.1f".format(rec.sequence.gcContent())}%")
    }

    // Write back to FASTA string
    val fastaOut = SeqIO.toString(records, "fasta")
    println("\nRe-serialized FASTA:\n$fastaOut")

    // ─── 7. Restriction Digest ──────────────────────────────
    println("▶ 7. RESTRICTION ENZYME ANALYSIS")
    println("-".repeat(40))

    val vector = DNASequence("AAGCTTGAATTCATGCGATCGATCGATCGATCGGATCCATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCATGCGATACGCTTACGATCGATCGTCTAGA")
    val sites = RestrictionEnzymes.digest(vector)
    println("Restriction sites in vector:")
    sites.forEach { println("  ${it.enzyme} at position ${it.position} (${it.recognition})") }

    val fragments = RestrictionEnzymes.virtualDigest(vector, listOf("EcoRI", "BamHI"))
    println("\nEcoRI + BamHI digest fragments:")
    fragments.forEachIndexed { i, frag -> println("  Fragment ${i+1}: ${frag.length} bp") }
    println()

    // ─── 8. Multiple Sequence Alignment ─────────────────────
    println("▶ 8. MULTIPLE SEQUENCE ALIGNMENT")
    println("-".repeat(40))

    val sequences = listOf(
        DNASequence("ATCGATCGATCG", id = "seq1"),
        DNASequence("ATCGAATCGATCG", id = "seq2"),
        DNASequence("ATCGATCGTTCG", id = "seq3"),
        DNASequence("ATCGATCGATGG", id = "seq4")
    )
    val msaAligner = MultipleSequenceAligner()
    val msa = msaAligner.align(sequences)
    println("MSA (${msa.numSeqs} seqs, ${msa.length} cols):")
    println(msa.prettyPrint())
    println("Consensus: ${msa.consensusSequence()}")
    println("% Identity: ${"%.1f".format(msa.percentIdentity())}%")
    println()

    // ─── 9. Phylogenetics ────────────────────────────────────
    println("▶ 9. PHYLOGENETICS")
    println("-".repeat(40))

    val distMatrix = DistanceMatrix(
        taxa = listOf("Human", "Chimp", "Gorilla", "Orangutan", "Gibbon"),
        distances = arrayOf(
            doubleArrayOf(0.0, 0.01, 0.03, 0.12, 0.18),
            doubleArrayOf(0.01, 0.0, 0.04, 0.13, 0.19),
            doubleArrayOf(0.03, 0.04, 0.0, 0.11, 0.17),
            doubleArrayOf(0.12, 0.13, 0.11, 0.0, 0.16),
            doubleArrayOf(0.18, 0.19, 0.17, 0.16, 0.0)
        )
    )

    println("Distance Matrix:\n$distMatrix")
    val tree = TreeBuilder.upgma(distMatrix)
    println("UPGMA Tree (Newick):\n  ${tree.toNewick()}")
    println("\nTree structure:")
    print(tree.printAscii())
    println()

    // ─── 10. Statistics ─────────────────────────────────────
    println("▶ 10. SEQUENCE STATISTICS")
    println("-".repeat(40))

    val seqs = listOf(
        DNASequence("ATGCGATACGCTTACGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGC"),
        DNASequence("ATCGATCG"),
        DNASequence("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
    )
    val nStats = SeqStats.nStats(seqs)
    println("Assembly stats:")
    println("  Total: ${nStats.totalLength} bp in ${nStats.numSeqs} seqs")
    println("  N50: ${nStats.n50}  L50: ${nStats.l50}")
    println("  Max: ${nStats.maxLen}  Min: ${nStats.minLen}  Mean: ${"%.1f".format(nStats.meanLen)}")

    println("\nEntropy of first seq: ${"%.3f".format(SeqStats.entropy(seqs[0]))}")
    println("Linguistic complexity: ${"%.3f".format(SeqStats.linguisticComplexity(seqs[0]))}")

    val codons = SeqStats.codonUsageTable(seqs[0])
    println("\nCodon usage (first 5):")
    codons.entries.take(5).forEach { (c, n) -> println("  $c: $n") }
    println()

    // ─── 11. IUPAC Ambiguity ────────────────────────────────
    println("▶ 11. IUPAC AMBIGUITY CODES")
    println("-".repeat(40))

    val ambig = "ATR"   // R = A or G
    val expanded = IUPAC.expandAmbiguous(ambig)
    println("Expanding '$ambig': ${expanded.joinToString(", ")}")
    println("Ambiguity code for {A,G,T}: ${IUPAC.ambiguityCode(setOf('A','G','T'))}")
    println()

    // ─── 12. Population Genetics ────────────────────────────
    println("▶ 12. POPULATION GENETICS (dN/dS)")
    println("-".repeat(40))

    val geneA = DNASequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
    val geneB = DNASequence("ATGGCCATCGTAATGGGACGCTGAAAGGGCGCCCGA")
    val (ka, ks) = PopGen.dNdS(geneA, geneB)
    println("dN (Ka) = ${"%.4f".format(ka)}")
    println("dS (Ks) = ${"%.4f".format(ks)}")
    println("dN/dS   = ${if (ks > 0) "%.4f".format(ka / ks) else "N/A"}")
    println()

    println("=".repeat(60))
    println("  BioKt demo complete!")
    println("=".repeat(60))
    println()
    runTests()
}
