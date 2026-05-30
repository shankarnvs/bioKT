package biokt

// Simple test runner — works without JUnit platform
fun runTests() {
    var passed = 0; var failed = 0; val failures = mutableListOf<String>()

    fun test(name: String, block: () -> Unit) {
        try { block(); passed++; println("  ✓ $name") }
        catch (e: Exception) { failed++; failures.add(name + ": " + e.message); println("  ✗ $name: ${e.message}") }
    }

    fun assertEquals(expected: Any?, actual: Any?, msg: String = "") {
        if (expected != actual) throw AssertionError("Expected $expected but was $actual $msg")
    }
    fun assertTrue(v: Boolean, msg: String = "") { if (!v) throw AssertionError("Expected true $msg") }
    fun assertFails(block: () -> Unit) {
        try { block(); throw AssertionError("Expected an exception but none was thrown") }
        catch (e: IllegalArgumentException) { /* expected */ }
    }

    println("\n=== BioKt Test Suite ===\n")

    println("DNA Sequence Tests:")
    test("basic construction") { assertEquals("ATGC", DNASequence("ATGC").sequence) }
    test("lowercase normalization") { assertEquals("ATGC", DNASequence("atgc").sequence) }
    test("invalid char throws") { assertFails { DNASequence("ATGX") } }
    test("complement") { assertEquals("TACG", DNASequence("ATGC").complement().sequence) }
    test("reverse complement") { assertEquals("GCAT", DNASequence("ATGC").reverseComplement().sequence) }
    test("gc content") { assertTrue(Math.abs(DNASequence("ATGC").gcContent() - 50.0) < 0.01) }
    test("count base") { assertEquals(3, DNASequence("AAATTTGC").count('A')) }
    test("transcription") { assertEquals("AUGCUU", DNASequence("ATGCTT").transcribe().sequence) }
    test("translation MAIVMGR") {
        val p = DNASequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG").translate()
        assertEquals("MAIVMGR", p.sequence)
    }
    test("find motif") { assertEquals(listOf(0, 3, 6), DNASequence("ATGATGATG").findMotif("ATG")) }
    test("melting temperature > 0") { assertTrue(DNASequence("ATGCAT").meltingTemperature() > 0) }
    test("find ORFs") {
        val seq = DNASequence("NNNATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAGNNNATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGTAA")
        val orfs = seq.findOrfs(minLength = 30)
        assertTrue(orfs.isNotEmpty(), "Should find ORFs")
    }

    println("\nRNA Sequence Tests:")
    test("RNA construction") { assertEquals("AUGCUU", RNASequence("AUGCUU").sequence) }
    test("back transcription") { assertEquals("ATGCTT", RNASequence("AUGCUU").backTranscribe().sequence) }
    test("RNA translation") { assertEquals("MAI", RNASequence("AUGGCCAUUUAA").translate().sequence) }
    test("RNA invalid base") { assertFails { RNASequence("AUGXT") } }

    println("\nProtein Sequence Tests:")
    test("molecular weight > 100") { assertTrue(ProteinSequence("M").molecularWeight() > 100) }
    test("amino acid composition") {
        val comp = ProteinSequence("AACCGG").aminoAcidComposition()
        assertTrue(Math.abs((comp['A'] ?: 0.0) - 33.33) < 0.1)
    }
    test("isoelectric point 0-14") {
        val pi = ProteinSequence("MKALVLLY").isoelectricPoint()
        assertTrue(pi in 0.0..14.0)
    }

    println("\nCodon Table Tests:")
    test("standard ATG=M") { assertEquals('M', CodonTables.STANDARD.translate("ATG")) }
    test("standard TAA=stop") { assertEquals('*', CodonTables.STANDARD.translate("TAA")) }
    test("all 64 codons valid") {
        for (a in "ATGC") for (b in "ATGC") for (c in "ATGC") {
            val aa = CodonTables.STANDARD.translate("$a$b$c")
            assertTrue(aa in Alphabets.AMINO_ACIDS, "Invalid AA '$aa' for $a$b$c")
        }
    }
    test("mito TGA=Trp, standard TGA=stop") {
        assertEquals('W', CodonTables.VERTEBRATE_MITOCHONDRIAL.translate("TGA"))
        assertEquals('*', CodonTables.STANDARD.translate("TGA"))
    }

    println("\nAlignment Tests:")
    test("global identical score > 0") {
        val r = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL).align("ATGC", "ATGC")
        assertTrue(r.score > 0)
    }
    test("global identity 100% for identical") {
        val r = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL).align("ATGC", "ATGC")
        assertTrue(Math.abs(r.identity - 1.0) < 0.01)
    }
    test("local alignment score > 0") {
        val r = PairwiseAligner(mode = PairwiseAligner.Mode.LOCAL).align("XXXATGCXXX", "ATGC")
        assertTrue(r.score > 0)
    }
    test("global alignment has gap for different-length seqs") {
        val r = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL).align("ATCGATCG", "ATCG")
        assertTrue(r.alignedA.contains('-') || r.alignedB.contains('-'))
    }

    println("\nSeqIO Tests:")
    test("parse FASTA 2 records") {
        val r = SeqIO.parseString(">seq1 desc\nATGCATGC\n>seq2\nGCTAGCTA", "fasta")
        assertEquals(2, r.size)
        assertEquals("seq1", r[0].id)
        assertEquals("ATGCATGC", r[0].sequence.sequence)
    }
    test("parse FASTA multiline") {
        val r = SeqIO.parseString(">seq1\nATGC\nATGC\nATGC", "fasta")
        assertEquals("ATGCATGCATGC", r[0].sequence.sequence)
    }
    test("FASTA roundtrip") {
        val records = listOf(SeqRecord(DNASequence("ATGCATGC"), "seq1", description = "test"))
        val parsed = SeqIO.parseString(SeqIO.toString(records, "fasta"), "fasta")
        assertEquals("ATGCATGC", parsed[0].sequence.sequence)
    }
    test("parse FASTQ") {
        val r = SeqIO.parseString("@seq1\nATGCATGC\n+\nIIIIIIII", "fastq")
        assertEquals("ATGCATGC", r[0].sequence.sequence)
    }

    println("\nRestriction Enzyme Tests:")
    test("EcoRI detection") {
        val sites = RestrictionEnzymes.digest(DNASequence("AAAGAATTCAAA"))
        assertTrue(sites.any { it.enzyme == "EcoRI" })
    }
    test("EcoRI position 0") {
        val sites = RestrictionEnzymes.digest(DNASequence("GAATTCATGCGGATCC"))
        assertEquals(0, sites.first { it.enzyme == "EcoRI" }.position)
    }
    test("virtual digest preserves length") {
        val dna = DNASequence("GAATTCATGCGGATCC")
        val frags = RestrictionEnzymes.virtualDigest(dna, listOf("EcoRI", "BamHI"))
        assertEquals(dna.length, frags.map { it.length }.sum())
    }

    println("\nPhylogenetics Tests:")
    test("UPGMA produces 3 leaves") {
        val m = DistanceMatrix(listOf("A","B","C"), arrayOf(
            doubleArrayOf(0.0,0.1,0.2), doubleArrayOf(0.1,0.0,0.2), doubleArrayOf(0.2,0.2,0.0)
        ))
        assertEquals(3, TreeBuilder.upgma(m).leaves().size)
    }
    test("Newick roundtrip contains taxa") {
        val t = PhyloTree.fromNewick("((A:0.1,B:0.1):0.05,C:0.2):0.0;")
        val nw = t.toNewick()
        assertTrue(nw.contains("A") && nw.contains("B") && nw.contains("C"))
    }

    println("\nIUPAC Tests:")
    test("expand R -> A,G") {
        val e = IUPAC.expandAmbiguous("R")
        assertTrue("A" in e && "G" in e && e.size == 2)
    }
    test("expand N -> 4 bases") { assertEquals(4, IUPAC.expandAmbiguous("N").size) }
    test("ambiguity code A,G -> R") { assertEquals('R', IUPAC.ambiguityCode(setOf('A','G'))) }

    println("\n==============================")
    println("Results: $passed passed, $failed failed out of ${passed+failed} tests")
    if (failures.isNotEmpty()) { println("\nFailures:"); failures.forEach { println("  - $it") } }
    println("==============================\n")
}
