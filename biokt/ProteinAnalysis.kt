package biokt

import kotlin.math.pow

/** ProtParam-compatible descriptors for a non-empty sequence of standard residues. */
class ProteinAnalysis(protein: ProteinSequence) {
    private val data = ProteinAnalysisData
    private val sequence = protein.sequence
    private val n = sequence.length
    private val indices = sequence.map { data.letters.indexOf(it) }
    init { require(n > 0 && indices.all { it >= 0 }) { "Expected standard amino acids" } }
    val counts = data.letters.associateWith { aa -> sequence.count { it == aa } }
    val percentages = counts.mapValues { 100.0 * it.value / n }
    private fun fraction(group: String) = group.sumOf { counts.getValue(it) }.toDouble() / n
    fun aromaticity() = fraction("FWY")
    fun secondaryStructureFractions() = listOf(fraction("EMALK"), fraction("NPGSD"), fraction("VIYFWLT"))
    fun extinctionCoefficients(): List<Int> {
        val reduced = 5500 * counts.getValue('W') + 1490 * counts.getValue('Y')
        return listOf(reduced, reduced + 125 * (counts.getValue('C') / 2))
    }
    fun molecularWeight(monoisotopic: Boolean = false): Double {
        val masses = if (monoisotopic) data.mono else data.mass
        val water = if (monoisotopic) 18.010565 else 18.0153
        return indices.sumOf { masses[it] } - (n - 1) * water
    }
    fun instabilityIndex() = (1 until n).sumOf { data.instability[indices[it-1]][indices[it]] } * (10.0 / n)
    fun gravy(scale: String = "KyteDoolitle") = indices.map { data.scales.getValue(scale)[it] }.average()
    fun flexibility(): List<Double> {
        // Preserve Biopython's profile length and its offset for the central contribution.
        val weights = listOf(0.25, 0.4375, 0.625, 0.8125)
        return (0 until n-9).map { start ->
            val sides = weights.indices.sumOf { j ->
                weights[j] * (data.flex[indices[start+j]] + data.flex[indices[start+8-j]])
            }
            (sides + data.flex[indices[start+5]]) / 5.25
        }
    }
    fun proteinScale(scale: Map<Char, Double>, window: Int, edge: Double = 1.0): List<Double> {
        require(window >= 2 && edge in 0.0..1.0)
        val half = window / 2
        val weights = (0 until half).map { edge + 2 * (1-edge) * it / (window-1) }
        val total = 1 + 2 * weights.sum()
        return (0..n-window).map { start ->
            val sides = weights.indices.sumOf { j ->
                val left = scale[sequence[start+j]]
                val right = scale[sequence[start+window-1-j]]
                if (left == null || right == null) 0.0 else weights[j] * (left+right)
            }
            (sides + (scale[sequence[start+half]] ?: 0.0)) / total
        }
    }
    fun proteinScale(scale: String, window: Int, edge: Double = 1.0): List<Double> {
        val values = data.scales.getValue(scale)
        return proteinScale(data.letters.mapIndexed { i, aa -> aa to values[i] }.toMap(), window, edge)
    }
    private val ionization = ProteinIonization(sequence, counts)
    fun chargeAtPH(ph: Double) = ionization.chargeAtPH(ph)
    fun isoelectricPoint() = ionization.isoelectricPoint()
}

internal class ProteinIonization(private val sequence: String, private val counts: Map<Char, Int>) {
    fun chargeAtPH(ph: Double): Double {
        val nTerm = when(sequence.firstOrNull()) {
            'A' -> 7.59; 'M' -> 7.0; 'S' -> 6.93; 'P' -> 8.36
            'T' -> 6.82; 'V' -> 7.44; 'E' -> 7.7; else -> 7.5
        }
        val cTerm = when(sequence.lastOrNull()) { 'D' -> 4.55; 'E' -> 4.75; else -> 3.55 }
        fun basic(pk: Double) = 1 / (1 + 10.0.pow(ph-pk))
        fun acidic(pk: Double) = 1 / (1 + 10.0.pow(pk-ph))
        return basic(nTerm) + counts.getOrDefault('K', 0) * basic(10.0) + counts.getOrDefault('R', 0) * basic(12.0) +
            counts.getOrDefault('H', 0) * basic(5.98) - acidic(cTerm) - counts.getOrDefault('D', 0) * acidic(4.05) -
            counts.getOrDefault('E', 0) * acidic(4.45) - counts.getOrDefault('C', 0) * acidic(9.0) - counts.getOrDefault('Y', 0) * acidic(10.0)
    }
    fun isoelectricPoint(): Double {
        var low = 4.05
        var high = 12.0
        var mid = 7.775
        while (high-low > 0.0001) {

            if (chargeAtPH(mid)>0) low=mid else high=mid
            mid = (low+high)/2
        }
        return mid
    }
}
