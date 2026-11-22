package biokt

import biokt.DrugProteinBinding.MeasuredAffinity
import kotlin.math.exp
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class DrugProteinCalibrationTest {
    private val target = ProteinSequence("ACDEFGHIKLMNPQRSTVWY", "reference")
    private val scorer = DrugProteinBinding
    private val ligands = listOf(Molecule("CCO"), Molecule("c1ccccc1"), Molecule("CC(=O)O"))
        .sortedBy { scorer.score(it, target).bindingScore }

    @Test
    fun interpolatesMeasuredKdAndPreservesTheHeuristicResult() {
        val training = listOf(ligands.first(), ligands.last()).map { ligand ->
            MeasuredAffinity(ligand, expectedKd(ligand))
        }
        val fitted = scorer.fitAffinityCalibration(target, training)
        val predicted = scorer.score(ligands[1], target, fitted)
        assertEquals(expectedKd(ligands[1]), predicted.kdMicromolar, 1e-10)
        assertEquals(scorer.score(ligands[1], target), predicted.heuristic)
        assertFalse(predicted.isExtrapolation)
    }

    @Test
    fun flagsExtrapolationAndRejectsOtherProteinSequences() {
        val training = ligands.drop(1).map { MeasuredAffinity(it, expectedKd(it)) }
        val fitted = scorer.fitAffinityCalibration(target, training)
        assertTrue(scorer.score(ligands.first(), target, fitted).isExtrapolation)
        assertFailsWith<IllegalArgumentException> {
            scorer.score(ligands.first(), ProteinSequence("AAAA", target.id), fitted)
        }
    }

    @Test
    fun rejectsInvalidKdAndDegenerateCalibrationInputs() {
        val reference = MeasuredAffinity(ligands.first(), 1.0)
        for (invalid in listOf(0.0, -1.0, Double.NaN, Double.POSITIVE_INFINITY)) {
            assertFailsWith<IllegalArgumentException> {
                scorer.fitAffinityCalibration(target, listOf(reference, MeasuredAffinity(ligands.last(), invalid)))
            }
        }
        assertFailsWith<IllegalArgumentException> {
            scorer.fitAffinityCalibration(target, listOf(reference, reference.copy(kdMicromolar = 2.0)))
        }
    }

    private fun expectedKd(ligand: Molecule): Double =
        exp(0.3 * scorer.score(ligand, target).bindingScore + 0.4)
}
