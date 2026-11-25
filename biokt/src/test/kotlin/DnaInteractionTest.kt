package biokt

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class DnaInteractionTest {
    @Test
    fun countsMismatchesAtBothEndsOfTheAlignedDuplex() {
        val strand = DNASequence("ACTGATGCTAGTCCAGTACG")
        val complement = strand.reverseComplement().sequence
        val perfect = DNAInteraction.hybridize(strand, DNASequence(complement))
        assertEquals(0, perfect.mismatches)
        for (position in listOf(complement.indices.first, complement.indices.last)) {
            val replacement = if (complement[position] == 'A') 'C' else 'A'
            val mismatched = DNASequence(complement.replaceRange(position, position + 1, replacement.toString()))
            val result = DNAInteraction.hybridize(strand, mismatched)
            assertEquals(1, result.mismatches, "Mismatch on complementary strand at index $position")
            assertTrue(result.tm < perfect.tm)
        }
    }
}
