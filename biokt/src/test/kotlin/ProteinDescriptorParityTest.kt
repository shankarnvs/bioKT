package biokt

import kotlin.test.Test
import kotlin.test.assertEquals

class ProteinDescriptorParityTest {
    @Test
    fun membraneInputsMatchProtParam() {
        val cases = listOf(
            Triple("MHFYWDE", 3.0 / 7, 4.710667991638183),
            Triple("AGK", 0.0, 8.795051002502444),
            Triple("DDDE", 0.0, 4.0500284194946286)
        )
        for ((sequence, aromaticity, pI) in cases) {
            val features = FeatureExtractor.membraneFeatures(ProteinSequence(sequence))
            assertEquals(aromaticity, features[4], 1e-8)
            assertEquals(pI, features[5], 1e-8)
        }
    }
}
