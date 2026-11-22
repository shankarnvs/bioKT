package biokt.rice

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CodingSnpAnnotationTest {
    private val forward = model('+', "1 4 0", "10 17 2")
    private val forwardDna = "ATGGNNNNNCTTGGTAA"

    @Test
    fun annotatesEachAlternateAlleleInACodonSplitAcrossExons() {
        val annotation = annotate(forward, forwardDna, 10, "C", "A,T")
        assertEquals("missense", annotation.consequence)
        val partial = annotate(forward, forwardDna, 11, "T", "A,N")
        assertEquals("coding_variant", partial.consequence)
        val synonymous = partial.codingSnpEffects.single()
        assertEquals(
            AnnotationOverlap.CodingSnpEffect("tx", "A", "missense", "GCT", "GAT", 2, 'A', 'D'),
            annotation.codingSnpEffects.first()
        )
        assertEquals("synonymous", synonymous.consequence)
        assertEquals("GCA", synonymous.alternateCodon)
        assertEquals(2, annotation.codingSnpEffects.size)
        assertEquals("GTT", annotation.codingSnpEffects[1].alternateCodon)
    }

    @Test
    fun translatesAReverseStrandCodonSplitAcrossExons() {
        val reverse = model('-', "4 11 2", "20 23 0")
        val dna = "NNNTTACCAAGNNNNNNNNCCAT"
        val effect = annotate(reverse, dna, 11, "G", "T").codingSnpEffects.single()
        assertEquals(
            AnnotationOverlap.CodingSnpEffect("tx", "T", "missense", "GCT", "GAT", 2, 'A', 'D'),
            effect
        )
    }

    @Test
    fun distinguishesStopGainAndLoss() {
        val gain = annotate(forward, forwardDna, 13, "G", "A").codingSnpEffects.single()
        assertEquals("stop_gained", gain.consequence)
        assertEquals('*', gain.alternateAminoAcid)
        val loss = annotate(forward, forwardDna, 16, "A", "C").codingSnpEffects.single()
        assertEquals("stop_lost", loss.consequence)
        assertEquals('S', loss.alternateAminoAcid)
    }

    @Test
    fun fallsBackWhenReferenceOrPhaseCannotSupportACall() {
        val vcf = variant(10, "C", "A")
        val invalidPhase = model('+', "1 4 0", "10 17 1")
        val unresolved = listOf(
            AnnotationOverlap.annotateVariants(vcf, forward).single(),
            annotate(forward, forwardDna, 10, "G", "A"),
            annotate(invalidPhase, forwardDna, 10, "C", "A"),
            annotate(forward, forwardDna.replaceRange(3, 4, "N"), 10, "C", "A")
        )
        for (result in unresolved) {
            assertEquals("missense_or_synonymous", result.consequence)
            assertTrue(result.codingSnpEffects.isEmpty())
        }
    }

    private fun model(strand: Char, vararg segments: String): AnnotationDb {
        val rows = mutableListOf(
            "chr1 test gene 1 24 . $strand . ID=g",
            "chr1 test mRNA 1 24 . $strand . ID=tx;Parent=g"
        )
        for (segment in segments) {
            val (start, end, phase) = segment.split(' ')
            rows.add("chr1 test CDS $start $end . $strand $phase Parent=tx")
        }
        return AnnotationParser.parseGff3String(rows.joinToString("\n") { it.replace(' ', '\t') })
    }

    private fun annotate(db: AnnotationDb, dna: String, pos: Int, ref: String, alt: String) =
        AnnotationOverlap.annotateVariants(variant(pos, ref, alt), db, mapOf("chr1" to dna)).single()

    private fun variant(pos: Int, ref: String, alt: String): VcfFile = VcfParser.parseString(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t$pos\t.\t$ref\t$alt\t.\tPASS\t."
    )
}
