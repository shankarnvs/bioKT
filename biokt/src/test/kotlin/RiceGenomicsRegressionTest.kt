package biokt.rice

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

class RiceGenomicsRegressionTest {
    @Test
    fun parsesDiploidCToTTransition() {
        val fixture = """
            ##fileformat=VCFv4.2
            #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT cc ct tt phased missing
            chr1 10 c_to_t C T 60 PASS . GT 0/0 0/1 1/1 1|0 ./.
        """.trimIndent().replace(' ', '\t')

        val snp = VcfParser.parseString(fixture).records.single()
        assertEquals("C", snp.ref)
        assertEquals(listOf("T"), snp.alt)
        assertTrue(snp.isTransition)
        assertEquals(listOf(0, 1, 2, 1, null),
            listOf("cc", "ct", "tt", "phased", "missing").map { snp.genotypes.getValue(it).dosage() })
        assertTrue(snp.genotypes.getValue("ct").isHet())
        assertTrue(snp.genotypes.getValue("phased").isPhased)
        assertEquals(listOf("1", "0"), snp.genotypes.getValue("phased").alleles())
        assertEquals(0.5, snp.hetFrequency(), 1e-12)
    }

    @Test
    fun parsesGenotypesAndFiltersByRegionAndMissingness() {
        val fixture = """
            ##fileformat=VCFv4.2
            #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 sample3
            chr1 10 snp A G 60 PASS . GT 0/1 1|1 ./.
            chr1 20 indel A AT 50 PASS . GT 0/0 0/1 0/1
            chr1 30 filtered A C 10 LowQual . GT 0/0 0/1 1/1
        """.trimIndent().replace(' ', '\t')

        val all = VcfParser.parseString(fixture)
        assertEquals(listOf("sample1", "sample2", "sample3"), all.header.samples)
        val snp = all.records.first()
        assertTrue(snp.genotypes.getValue("sample1").isHet())
        assertTrue(snp.genotypes.getValue("sample2").isHomAlt())
        assertNull(snp.genotypes.getValue("sample3").dosage())
        assertEquals(0.25, snp.maf(), 1e-12)
        assertEquals(0.5, snp.hetFrequency(), 1e-12)

        val passing = VcfParser.parseString(
            fixture, regionChrom = "chr1", regionStart = 10, regionEnd = 30,
            passOnly = true, maxMissing = 0.25
        )
        assertEquals(listOf("indel"), passing.records.map { it.id })
    }

    @Test
    fun computesLengthAdjustedTpmAndKeepsZeroCountSamplesAtZero() {
        val counts = arrayOf(doubleArrayOf(100.0, 120.0, 0.0), doubleArrayOf(100.0, 240.0, 0.0))
        val lengths = doubleArrayOf(1000.0, 2000.0)
        val expected = arrayOf(
            doubleArrayOf(2_000_000.0 / 3, 500_000.0, 0.0),
            doubleArrayOf(1_000_000.0 / 3, 500_000.0, 0.0)
        )

        val tpm = ExpressionNormaliser.toTPM(counts, lengths)
        expected.forEachIndexed { gene, row ->
            row.forEachIndexed { sample, value ->
                assertEquals(value, tpm[gene][sample], 1e-8, "Gene $gene, sample $sample")
            }
        }
    }
}
