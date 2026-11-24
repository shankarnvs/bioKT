package biokt

import kotlin.test.Test
import kotlin.test.assertEquals

class LocalMembraneFeaturesTest {
    @Test
    fun countsDistinctStretchesAndTheirCoverage() {
        val segment = "I".repeat(19)
        val sequence = ProteinSequence(segment + "D".repeat(40) + segment, "test", "")
        val features = FeatureExtractor.localMembraneFeatures(sequence)
        assertEquals(4.5, features[8], 1e-10)
        assertEquals(2.0, features[1])
        assertEquals(listOf(0..27, 50..77), FeatureExtractor.membraneCandidateRegions(sequence))
        assertEquals(56.0 / 78, features[2], 1e-10)
        assertEquals(9.0 / 28, features[3], 1e-10)
        assertEquals(1.0, features[4])
        assertEquals(0.0, features[5])
        assertEquals(1.0, features[6])
    }
}
