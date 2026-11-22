package biokt.rice

import kotlin.test.Test
import kotlin.test.assertEquals

class AnnotationOverlapTest {
    private val db = AnnotationParser.parseGff3String(
        """
        chr1 test gene 100 240 . + . ID=gene1
        chr1 test mRNA 100 240 . + . ID=tx1;Parent=gene1
        chr1 test exon 110 130 . + . Parent=tx1
        chr1 test CDS 110 130 . + 0 Parent=tx1
        chr1 test exon 200 230 . + . Parent=tx1
        """.trimIndent().replace(' ', '\t')
    )

    @Test
    fun distinguishesCodingIndelsAndOtherVariants() {
        assertEquals("inframe_indel", consequence(115, "A", "ATCG"))
        assertEquals("inframe_indel", consequence(115, "ATCG", "A"))
        assertEquals("frameshift", consequence(115, "A", "AT"))
        assertEquals("coding_variant", consequence(115, "AT", "CG"))
        assertEquals("coding_variant", consequence(115, "A", "<DEL>"))
        assertEquals("coding_variant", consequence(115, "A", "AT,ATCG"))
    }

    @Test
    fun usesBothInternalExonBoundaries() {
        for (position in listOf(131, 132, 198, 199)) {
            assertEquals("splice_region", consequence(position, "A", "G"))
        }
        for (position in listOf(109, 133, 197, 231)) {
            assertEquals("intronic", consequence(position, "A", "G"))
        }
        assertEquals("splice_region", consequence(131, "A", "AT"))
        assertEquals("exon_SNP", consequence(200, "A", "G"))
    }

    private fun consequence(position: Int, ref: String, alt: String): String {
        val vcf = VcfParser.parseString(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                "chr1\t$position\t.\t$ref\t$alt\t.\tPASS\t."
        )
        return AnnotationOverlap.annotateVariants(vcf, db).single().consequence
    }
}
