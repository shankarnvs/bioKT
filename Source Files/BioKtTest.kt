package biokt

import kotlin.test.*
import org.junit.jupiter.api.Test

class DNASequenceTest {

    @Test fun testBasicConstruction() {
        val seq = DNASequence("ATGC")
        assertEquals(4, seq.length)
        assertEquals("ATGC", seq.sequence)
    }

    @Test fun testLowerCaseNormalization() {
        val seq = DNASequence("atgc")
        assertEquals("ATGC", seq.sequence)
    }

    @Test fun testInvalidCharThrows() {
        assertFailsWith<IllegalArgumentException> { DNASequence("ATGX") }
    }

    @Test fun testComplement() {
        val seq = DNASequence("ATGC")
        assertEquals("TACG", seq.complement().sequence)
    }

    @Test fun testReverseComplement() {
        val seq = DNASequence("ATGC")
        assertEquals("GCAT", seq.reverseComplement().sequence)
    }

    @Test fun testGcContent() {
        val seq = DNASequence("ATGC")  // 2/4 = 50%
        assertTrue(Math.abs(seq.gcContent() - 50.0) < 0.01)
    }

    @Test fun testCountBase() {
        val seq = DNASequence("AAATTTGC")
        assertEquals(3, seq.count('A'))
        assertEquals(3, seq.count('T'))
        assertEquals(1, seq.count('G'))
    }

    @Test fun testTranscription() {
        val dna = DNASequence("ATGCTT")
        val rna = dna.transcribe()
        assertEquals("AUGCUU", rna.sequence)
    }

    @Test fun testTranslation() {
        // ATG=M, GCC=A, ATT=I, GTA=V, ATG=M, GGC=G, CGC=R, TGA=*
        val dna = DNASequence("ATGGCCATTGTAATGGGCCGCTGA")
        val protein = dna.translate()
        assertEquals("MAIMMGR", protein.sequence)
    }

    @Test fun testTranslateFromStartCodon() {
        val dna = DNASequence("NNNATGGCCATTTGA")
        val protein = dna.translate(startCodon = true)
        assertEquals("MA", protein.sequence)
    }

    @Test fun testFindMotif() {
        val seq = DNASequence("ATGATGATG")
        val positions = seq.findMotif("ATG")
        assertEquals(listOf(0, 3, 6), positions)
    }

    @Test fun testMeltingTemperature() {
        val shortSeq = DNASequence("ATGCAT")  // short, uses 2AT + 4GC formula
        val tm = shortSeq.meltingTemperature()
        assertTrue(tm > 0)
    }

    @Test fun testFindOrfs() {
        val seq = DNASequence("XXXATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGTAAXXX".replace("X","N"))
        val orfs = seq.findOrfs(minLength = 10)
        assertTrue(orfs.isNotEmpty())
    }
}

class RNASequenceTest {

    @Test fun testRNAConstruction() {
        val rna = RNASequence("AUGCUU")
        assertEquals("AUGCUU", rna.sequence)
    }

    @Test fun testBackTranscription() {
        val rna = RNASequence("AUGCUU")
        assertEquals("ATGCTT", rna.backTranscribe().sequence)
    }

    @Test fun testRNATranslation() {
        val rna = RNASequence("AUGGCCAUUUAA")
        val prot = rna.translate()
        assertEquals("MA", prot.sequence)
    }

    @Test fun testRNAInvalidBase() {
        assertFailsWith<IllegalArgumentException> { RNASequence("AUGXT") }
    }
}

class ProteinSequenceTest {

    @Test fun testMolecularWeight() {
        val prot = ProteinSequence("M")  // Met = 149.21
        assertTrue(prot.molecularWeight() > 100)
    }

    @Test fun testAminoAcidComposition() {
        val prot = ProteinSequence("AACCGG")
        val comp = prot.aminoAcidComposition()
        assertTrue(Math.abs((comp['A'] ?: 0.0) - 33.33) < 0.1)
    }

    @Test fun testIsoelectricPoint() {
        val prot = ProteinSequence("MKALVLLY")
        val pi = prot.isoelectricPoint()
        assertTrue(pi in 0.0..14.0)
    }
}

class CodonTableTest {

    @Test fun testStandardTranslation() {
        assertEquals('M', CodonTables.STANDARD.translate("ATG"))
        assertEquals('*', CodonTables.STANDARD.translate("TAA"))
        assertEquals('*', CodonTables.STANDARD.translate("TAG"))
        assertEquals('*', CodonTables.STANDARD.translate("TGA"))
    }

    @Test fun testAllSixtyFourCodons() {
        val bases = listOf('A','T','G','C')
        var count = 0
        for (a in bases) for (b in bases) for (c in bases) {
            val codon = "$a$b$c"
            val aa = CodonTables.STANDARD.translate(codon)
            assertTrue(aa in Alphabets.AMINO_ACIDS, "Invalid AA '$aa' for codon $codon")
            count++
        }
        assertEquals(64, count)
    }

    @Test fun testMitochondrialOverrides() {
        // TGA = Trp in vertebrate mito, Stop in standard
        assertEquals('W', CodonTables.VERTEBRATE_MITOCHONDRIAL.translate("TGA"))
        assertEquals('*', CodonTables.STANDARD.translate("TGA"))
    }
}

class AlignmentTest {

    @Test fun testGlobalAlignmentIdentical() {
        val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL)
        val result = aligner.align("ATGC", "ATGC")
        assertTrue(Math.abs(result.identity - 1.0) < 0.01)
    }

    @Test fun testGlobalAlignmentScore() {
        val aligner = PairwiseAligner(
            mode = PairwiseAligner.Mode.GLOBAL, matchScore = 2.0, mismatchScore = -1.0, gapOpen = -2.0
        )
        val result = aligner.align("ATGC", "ATGC")
        assertTrue(result.score > 0)
    }

    @Test fun testLocalAlignment() {
        val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.LOCAL)
        val result = aligner.align("XXXATGCXXX", "ATGC")
        assertTrue(result.score > 0)
    }

    @Test fun testAlignmentWithGap() {
        val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL)
        val result = aligner.align("ATGC", "AGC")
        assertTrue(result.gaps > 0)
    }

    @Test fun testGlobalAlignmentContainsGap() {
        val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL)
        val result = aligner.align("ATCGATCG", "ATCG")
        assertTrue(result.alignedA.contains('-') || result.alignedB.contains('-'))
    }
}

class SeqIOTest {

    @Test fun testParseFasta() {
        val fasta = ">seq1 description\nATGCATGC\n>seq2\nGCTAGCTA"
        val records = SeqIO.parseString(fasta, "fasta")
        assertEquals(2, records.size)
        assertEquals("seq1", records[0].id)
        assertEquals("description", records[0].description)
        assertEquals("ATGCATGC", records[0].sequence.sequence)
    }

    @Test fun testParseFastaMultiline() {
        val fasta = ">seq1\nATGC\nATGC\nATGC"
        val records = SeqIO.parseString(fasta, "fasta")
        assertEquals("ATGCATGCATGC", records[0].sequence.sequence)
    }

    @Test fun testWriteAndParseFasta() {
        val records = listOf(
            SeqRecord(DNASequence("ATGCATGC"), "seq1", description = "test")
        )
        val fastaStr = SeqIO.toString(records, "fasta")
        val parsed = SeqIO.parseString(fastaStr, "fasta")
        assertEquals(1, parsed.size)
        assertEquals("ATGCATGC", parsed[0].sequence.sequence)
    }

    @Test fun testParseFastq() {
        val fastq = "@seq1\nATGCATGC\n+\nIIIIIIII"
        val records = SeqIO.parseString(fastq, "fastq")
        assertEquals(1, records.size)
        assertEquals("ATGCATGC", records[0].sequence.sequence)
    }
}

class RestrictionTest {

    @Test fun testEcoRIRecognition() {
        val dna = DNASequence("AAAGAATTCAAA")
        val sites = RestrictionEnzymes.digest(dna)
        assertTrue(sites.any { it.enzyme == "EcoRI" })
    }

    @Test fun testDigestPosition() {
        val dna = DNASequence("GAATTCATGCGGATCC")
        val sites = RestrictionEnzymes.digest(dna)
        val ecoRI = sites.first { it.enzyme == "EcoRI" }
        assertEquals(0, ecoRI.position)
    }

    @Test fun testVirtualDigest() {
        val dna = DNASequence("GAATTCATGCGGATCC")
        val fragments = RestrictionEnzymes.virtualDigest(dna, listOf("EcoRI", "BamHI"))
        assertTrue(fragments.size >= 2)
        assertEquals(dna.length, fragments.map { it.length }.sum())
    }
}

class PhylogeneticsTest {

    @Test fun testUpgmaTree() {
        val matrix = DistanceMatrix(
            taxa = listOf("A", "B", "C"),
            distances = arrayOf(
                doubleArrayOf(0.0, 0.1, 0.2),
                doubleArrayOf(0.1, 0.0, 0.2),
                doubleArrayOf(0.2, 0.2, 0.0)
            )
        )
        val tree = TreeBuilder.upgma(matrix)
        assertNotNull(tree.root)
        assertEquals(3, tree.leaves().size)
    }

    @Test fun testNewickRoundtrip() {
        val newick = "((A:0.1,B:0.1):0.05,C:0.2):0.0;"
        val tree = PhyloTree.fromNewick(newick)
        val output = tree.toNewick()
        assertNotNull(output)
        assertTrue(output.contains("A"))
        assertTrue(output.contains("B"))
        assertTrue(output.contains("C"))
    }
}

class IUPACTest {

    @Test fun testExpandR() {
        val expanded = IUPAC.expandAmbiguous("R")
        assertTrue("A" in expanded && "G" in expanded)
        assertEquals(2, expanded.size)
    }

    @Test fun testExpandN() {
        val expanded = IUPAC.expandAmbiguous("N")
        assertEquals(4, expanded.size)
    }

    @Test fun testAmbiguityCode() {
        assertEquals('R', IUPAC.ambiguityCode(setOf('A', 'G')))
        assertEquals('N', IUPAC.ambiguityCode(setOf('A', 'T', 'G', 'C')))
    }
}
