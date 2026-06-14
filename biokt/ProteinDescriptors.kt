package biokt

import kotlin.math.*

// ============================================================
// PROTEIN DESCRIPTORS
// Mirrors & extends pyBioMed's PyPro module
// AAC, DPC, TPC, CTD, QSO, PAAC, APAAC, SOCN, Geary, Moran
// ============================================================

object ProteinDescriptors {

    // ── Amino Acid Composition (AAC) ─────────────────────────

    fun aminoAcidComposition(seq: ProteinSequence): Map<String, Double> {
        val aas = "ACDEFGHIKLMNPQRSTVWY"
        val counts = seq.sequence.groupingBy { it }.eachCount()
        return aas.associate { aa ->
            aa.toString() to (counts[aa] ?: 0).toDouble() / seq.length
        }
    }

    // ── Dipeptide Composition (DPC) ──────────────────────────

    fun dipeptideComposition(seq: ProteinSequence): Map<String, Double> {
        val aas = "ACDEFGHIKLMNPQRSTVWY"
        val counts = mutableMapOf<String, Int>()
        for (i in 0 until seq.length - 1) {
            val dp = "${seq[i]}${seq[i + 1]}"
            if (seq[i] in aas && seq[i + 1] in aas)
                counts[dp] = (counts[dp] ?: 0) + 1
        }
        val total = (seq.length - 1).toDouble()
        val result = mutableMapOf<String, Double>()
        for (a in aas) for (b in aas) result["$a$b"] = (counts["$a$b"] ?: 0) / total
        return result
    }

    // ── Tripeptide Composition (TPC) ─────────────────────────

    fun tripeptideComposition(seq: ProteinSequence): Map<String, Double> {
        val aas = "ACDEFGHIKLMNPQRSTVWY"
        val counts = mutableMapOf<String, Int>()
        for (i in 0 until seq.length - 2) {
            val tp = "${seq[i]}${seq[i+1]}${seq[i+2]}"
            if (tp.all { it in aas }) counts[tp] = (counts[tp] ?: 0) + 1
        }
        val total = (seq.length - 2).toDouble()
        return counts.mapValues { (_, v) -> v / total }
    }

    // ── CTD (Composition, Transition, Distribution) ──────────

    object CTDProperties {
        // Hydrophobicity (Eisenberg 1984)
        val HYDROPHOBICITY = mapOf(
            'A' to 1, 'R' to 2, 'N' to 2, 'D' to 2, 'C' to 1,
            'Q' to 2, 'E' to 2, 'G' to 2, 'H' to 3, 'I' to 1,
            'L' to 1, 'K' to 2, 'M' to 1, 'F' to 1, 'P' to 2,
            'S' to 2, 'T' to 2, 'W' to 1, 'Y' to 3, 'V' to 1
        )
        // Normalized van der Waals volume
        val VDW_VOLUME = mapOf(
            'A' to 1, 'R' to 3, 'N' to 2, 'D' to 2, 'C' to 2,
            'Q' to 3, 'E' to 3, 'G' to 1, 'H' to 3, 'I' to 3,
            'L' to 3, 'K' to 3, 'M' to 3, 'F' to 3, 'P' to 2,
            'S' to 1, 'T' to 2, 'W' to 3, 'Y' to 3, 'V' to 2
        )
        // Polarity
        val POLARITY = mapOf(
            'A' to 1, 'R' to 3, 'N' to 3, 'D' to 3, 'C' to 2,
            'Q' to 3, 'E' to 3, 'G' to 1, 'H' to 2, 'I' to 1,
            'L' to 1, 'K' to 3, 'M' to 1, 'F' to 1, 'P' to 2,
            'S' to 2, 'T' to 2, 'W' to 2, 'Y' to 2, 'V' to 1
        )
        // Secondary structure (helix, sheet, coil tendencies)
        val SEC_STRUCTURE = mapOf(
            'A' to 1, 'R' to 1, 'N' to 3, 'D' to 3, 'C' to 2,
            'Q' to 1, 'E' to 1, 'G' to 3, 'H' to 2, 'I' to 2,
            'L' to 1, 'K' to 1, 'M' to 1, 'F' to 2, 'P' to 3,
            'S' to 3, 'T' to 2, 'W' to 2, 'Y' to 2, 'V' to 2
        )
        // Charge at pH 7
        val CHARGE = mapOf(
            'A' to 2, 'R' to 3, 'N' to 2, 'D' to 1, 'C' to 2,
            'Q' to 2, 'E' to 1, 'G' to 2, 'H' to 3, 'I' to 2,
            'L' to 2, 'K' to 3, 'M' to 2, 'F' to 2, 'P' to 2,
            'S' to 2, 'T' to 2, 'W' to 2, 'Y' to 2, 'V' to 2
        )

        val ALL = mapOf(
            "hydrophobicity" to HYDROPHOBICITY,
            "vdwVolume" to VDW_VOLUME,
            "polarity" to POLARITY,
            "secStructure" to SEC_STRUCTURE,
            "charge" to CHARGE
        )
    }

    fun ctdDescriptors(seq: ProteinSequence): Map<String, Double> {
        val result = mutableMapOf<String, Double>()
        val s = seq.sequence

        for ((propName, propMap) in CTDProperties.ALL) {
            val encoded = s.mapNotNull { propMap[it] }
            val len = encoded.size.toDouble()

            // Composition: fraction in each class
            for (cls in 1..3) {
                result["${propName}_C$cls"] = encoded.count { it == cls } / len
            }

            // Transition: fraction of adjacent pairs that differ in class
            var transitions = 0
            for (i in 0 until encoded.size - 1)
                if (encoded[i] != encoded[i + 1]) transitions++
            result["${propName}_T"] = transitions / (len - 1)

            // Distribution: position of 1st, 25%, 50%, 75%, 100% of each class
            for (cls in 1..3) {
                val positions = encoded.indices.filter { encoded[it] == cls }
                if (positions.isEmpty()) {
                    listOf("D1","D25","D50","D75","D100").forEach { result["${propName}_${cls}_$it"] = 0.0 }
                    continue
                }
                val n = positions.size
                result["${propName}_${cls}_D1"]   = positions.first() / len
                result["${propName}_${cls}_D25"]  = positions[maxOf(0, n / 4)] / len
                result["${propName}_${cls}_D50"]  = positions[maxOf(0, n / 2)] / len
                result["${propName}_${cls}_D75"]  = positions[maxOf(0, 3 * n / 4)] / len
                result["${propName}_${cls}_D100"] = positions.last() / len
            }
        }
        return result
    }

    // ── Quasi-Sequence-Order (QSO) ───────────────────────────

    private val SCHNEIDER_WREDE_MATRIX: Map<Pair<Char, Char>, Double> by lazy {
        // Physicochemical distance matrix (Schneider & Wrede 1994)
        val aas = "ACDEFGHIKLMNPQRSTVWY"
        val data = arrayOf(
            doubleArrayOf(0.00,0.11,0.83,0.79,0.11,0.56,0.79,0.58,0.83,0.11,0.11,0.83,0.49,0.79,0.79,0.18,0.11,0.11,0.79,0.11),
            doubleArrayOf(0.11,0.00,0.94,0.90,0.00,0.67,0.90,0.69,0.94,0.00,0.00,0.94,0.60,0.90,0.90,0.29,0.22,0.00,0.90,0.22),
            doubleArrayOf(0.83,0.94,0.00,0.04,0.94,0.27,0.04,0.25,0.00,0.94,0.94,0.00,0.34,0.04,0.04,0.65,0.72,0.94,0.04,0.72),
            doubleArrayOf(0.79,0.90,0.04,0.00,0.90,0.23,0.00,0.21,0.04,0.90,0.90,0.04,0.30,0.00,0.00,0.61,0.68,0.90,0.00,0.68),
            doubleArrayOf(0.11,0.00,0.94,0.90,0.00,0.67,0.90,0.69,0.94,0.00,0.00,0.94,0.60,0.90,0.90,0.29,0.22,0.00,0.90,0.22),
            doubleArrayOf(0.56,0.67,0.27,0.23,0.67,0.00,0.23,0.02,0.27,0.67,0.67,0.27,0.07,0.23,0.23,0.38,0.45,0.67,0.23,0.45),
            doubleArrayOf(0.79,0.90,0.04,0.00,0.90,0.23,0.00,0.21,0.04,0.90,0.90,0.04,0.30,0.00,0.00,0.61,0.68,0.90,0.00,0.68),
            doubleArrayOf(0.58,0.69,0.25,0.21,0.69,0.02,0.21,0.00,0.25,0.69,0.69,0.25,0.09,0.21,0.21,0.40,0.47,0.69,0.21,0.47),
            doubleArrayOf(0.83,0.94,0.00,0.04,0.94,0.27,0.04,0.25,0.00,0.94,0.94,0.00,0.34,0.04,0.04,0.65,0.72,0.94,0.04,0.72),
            doubleArrayOf(0.11,0.00,0.94,0.90,0.00,0.67,0.90,0.69,0.94,0.00,0.00,0.94,0.60,0.90,0.90,0.29,0.22,0.00,0.90,0.22),
            doubleArrayOf(0.11,0.00,0.94,0.90,0.00,0.67,0.90,0.69,0.94,0.00,0.00,0.94,0.60,0.90,0.90,0.29,0.22,0.00,0.90,0.22),
            doubleArrayOf(0.83,0.94,0.00,0.04,0.94,0.27,0.04,0.25,0.00,0.94,0.94,0.00,0.34,0.04,0.04,0.65,0.72,0.94,0.04,0.72),
            doubleArrayOf(0.49,0.60,0.34,0.30,0.60,0.07,0.30,0.09,0.34,0.60,0.60,0.34,0.00,0.30,0.30,0.31,0.38,0.60,0.30,0.38),
            doubleArrayOf(0.79,0.90,0.04,0.00,0.90,0.23,0.00,0.21,0.04,0.90,0.90,0.04,0.30,0.00,0.00,0.61,0.68,0.90,0.00,0.68),
            doubleArrayOf(0.79,0.90,0.04,0.00,0.90,0.23,0.00,0.21,0.04,0.90,0.90,0.04,0.30,0.00,0.00,0.61,0.68,0.90,0.00,0.68),
            doubleArrayOf(0.18,0.29,0.65,0.61,0.29,0.38,0.61,0.40,0.65,0.29,0.29,0.65,0.31,0.61,0.61,0.00,0.07,0.29,0.61,0.07),
            doubleArrayOf(0.11,0.22,0.72,0.68,0.22,0.45,0.68,0.47,0.72,0.22,0.22,0.72,0.38,0.68,0.68,0.07,0.00,0.22,0.68,0.00),
            doubleArrayOf(0.11,0.00,0.94,0.90,0.00,0.67,0.90,0.69,0.94,0.00,0.00,0.94,0.60,0.90,0.90,0.29,0.22,0.00,0.90,0.22),
            doubleArrayOf(0.79,0.90,0.04,0.00,0.90,0.23,0.00,0.21,0.04,0.90,0.90,0.04,0.30,0.00,0.00,0.61,0.68,0.90,0.00,0.68),
            doubleArrayOf(0.11,0.22,0.72,0.68,0.22,0.45,0.68,0.47,0.72,0.22,0.22,0.72,0.38,0.68,0.68,0.07,0.00,0.22,0.68,0.00)
        )
        val m = mutableMapOf<Pair<Char, Char>, Double>()
        aas.forEachIndexed { i, a -> aas.forEachIndexed { j, b -> m[a to b] = data[i][j] } }
        m
    }

    private fun qsoCorrelation(seq: String, lag: Int, matrix: Map<Pair<Char, Char>, Double>): Double {
        if (lag >= seq.length) return 0.0
        var sum = 0.0; var count = 0
        for (i in 0 until seq.length - lag) {
            val d = matrix[seq[i] to seq[i + lag]]
            if (d != null) { sum += d * d; count++ }
        }
        return if (count == 0) 0.0 else sum / count
    }

    fun quasiSequenceOrder(seq: ProteinSequence, maxLag: Int = 30, weight: Double = 0.1): Map<String, Double> {
        val result = mutableMapOf<String, Double>()
        val s = seq.sequence

        val correlations = (1..minOf(maxLag, s.length - 1))
            .map { lag -> qsoCorrelation(s, lag, SCHNEIDER_WREDE_MATRIX) }

        val tau = correlations.sum()
        val aas = "ACDEFGHIKLMNPQRSTVWY"
        val counts = s.groupingBy { it }.eachCount()
        val totalAAC = aas.map { aa -> (counts[aa] ?: 0).toDouble() }.sum()
        val denom = totalAAC + weight * tau

        for (aa in aas) {
            result["QSO_${aa}"] = (counts[aa] ?: 0) / denom
        }
        for ((i, corr) in correlations.withIndex()) {
            result["QSO_lag${i + 1}"] = weight * corr / denom
        }
        return result
    }

    // ── Pseudo Amino Acid Composition (PAAC) ─────────────────

    fun pseudoAAC(seq: ProteinSequence, lambdaMax: Int = 30, weight: Double = 0.05): Map<String, Double> {
        val result = mutableMapOf<String, Double>()
        val s = seq.sequence
        val aas = "ACDEFGHIKLMNPQRSTVWY"

        // Physicochemical properties (Kuo-Chen Chou)
        val hydrophobicity = mapOf('A' to 0.62,'C' to 0.29,'D' to -0.90,'E' to -0.74,'F' to 1.19,
            'G' to 0.48,'H' to -0.40,'I' to 1.38,'K' to -1.50,'L' to 1.06,'M' to 0.64,
            'N' to -0.78,'P' to 0.12,'Q' to -0.85,'R' to -2.53,'S' to -0.18,'T' to -0.05,
            'V' to 1.08,'W' to 0.81,'Y' to 0.26)

        val hydrophilicity = mapOf('A' to -0.5,'C' to -1.0,'D' to 3.0,'E' to 3.0,'F' to -2.5,
            'G' to 0.0,'H' to -0.5,'I' to -1.8,'K' to 3.0,'L' to -1.8,'M' to -1.3,
            'N' to 0.2,'P' to 0.0,'Q' to 0.2,'R' to 3.0,'S' to 0.3,'T' to -0.4,
            'V' to -1.5,'W' to -3.4,'Y' to -2.3)

        val sideChainMass = mapOf('A' to 15.0,'C' to 47.0,'D' to 59.0,'E' to 73.0,'F' to 91.0,
            'G' to 1.0,'H' to 82.0,'I' to 57.0,'K' to 73.0,'L' to 57.0,'M' to 75.0,
            'N' to 58.0,'P' to 42.0,'Q' to 72.0,'R' to 100.0,'S' to 31.0,'T' to 45.0,
            'V' to 43.0,'W' to 130.0,'Y' to 107.0)

        fun normalize(prop: Map<Char, Double>): Map<Char, Double> {
            val mean = prop.values.average()
            val std = sqrt(prop.values.map { (it - mean).pow(2) }.average())
            return prop.mapValues { (_, v) -> if (std == 0.0) 0.0 else (v - mean) / std }
        }

        val nHydro = normalize(hydrophobicity)
        val nHydroPhil = normalize(hydrophilicity)
        val nMass = normalize(sideChainMass)

        fun theta(lag: Int): Double {
            var sum = 0.0; var n = 0
            for (i in 0 until s.length - lag) {
                val a = s[i]; val b = s[i + lag]
                val h1 = nHydro[a] ?: 0.0; val h2 = nHydro[b] ?: 0.0
                val p1 = nHydroPhil[a] ?: 0.0; val p2 = nHydroPhil[b] ?: 0.0
                val m1 = nMass[a] ?: 0.0; val m2 = nMass[b] ?: 0.0
                sum += ((h1 - h2).pow(2) + (p1 - p2).pow(2) + (m1 - m2).pow(2)) / 3.0
                n++
            }
            return if (n == 0) 0.0 else sum / n
        }

        val maxL = minOf(lambdaMax, s.length - 1)
        val thetas = (1..maxL).map { theta(it) }
        val thetaSum = thetas.sum()

        val counts = s.groupingBy { it }.eachCount()
        val denom = s.length + weight * thetaSum

        for (aa in aas) result["PAAC_${aa}"] = (counts[aa] ?: 0) / denom
        for ((i, t) in thetas.withIndex()) result["PAAC_lambda${i + 1}"] = weight * t / denom

        return result
    }

    // ── Geary Autocorrelation ────────────────────────────────

    fun gearyAutocorrelation(seq: ProteinSequence, maxLag: Int = 30): Map<String, Double> {
        val props = mapOf(
            "MW"   to mapOf('A' to 71.03,'C' to 103.01,'D' to 115.03,'E' to 129.04,'F' to 147.07,'G' to 57.02,'H' to 137.06,'I' to 113.08,'K' to 128.09,'L' to 113.08,'M' to 131.04,'N' to 114.04,'P' to 97.05,'Q' to 128.06,'R' to 156.10,'S' to 87.03,'T' to 101.05,'V' to 99.07,'W' to 186.08,'Y' to 163.06),
            "pK1"  to mapOf('A' to 2.34,'C' to 1.96,'D' to 1.88,'E' to 2.19,'F' to 1.83,'G' to 2.34,'H' to 1.82,'I' to 2.36,'K' to 2.18,'L' to 2.36,'M' to 2.28,'N' to 2.02,'P' to 1.99,'Q' to 2.17,'R' to 2.17,'S' to 2.21,'T' to 2.09,'V' to 2.32,'W' to 2.83,'Y' to 2.20)
        )
        val result = mutableMapOf<String, Double>()
        val s = seq.sequence

        for ((propName, propMap) in props) {
            val values = s.mapNotNull { propMap[it] }
            val mean = values.average()
            val n = values.size.toDouble()
            val variance = values.map { (it - mean).pow(2) }.sum() / n

            for (lag in 1..minOf(maxLag, values.size - 1)) {
                var num = 0.0
                for (i in 0 until values.size - lag)
                    num += (values[i] - values[i + lag]).pow(2)
                val geary = if (variance == 0.0) 0.0 else
                    (n / (2 * (n - lag))) * (num / ((n - 1) * variance))
                result["Geary_${propName}_lag$lag"] = geary
            }
        }
        return result
    }

    // ── Moran Autocorrelation ────────────────────────────────

    fun moranAutocorrelation(seq: ProteinSequence, maxLag: Int = 30): Map<String, Double> {
        val hydrophobicity = mapOf('A' to 0.62,'C' to 0.29,'D' to -0.90,'E' to -0.74,'F' to 1.19,
            'G' to 0.48,'H' to -0.40,'I' to 1.38,'K' to -1.50,'L' to 1.06,'M' to 0.64,
            'N' to -0.78,'P' to 0.12,'Q' to -0.85,'R' to -2.53,'S' to -0.18,'T' to -0.05,
            'V' to 1.08,'W' to 0.81,'Y' to 0.26)
        val result = mutableMapOf<String, Double>()
        val s = seq.sequence
        val values = s.mapNotNull { hydrophobicity[it] }
        val mean = values.average()
        val n = values.size.toDouble()
        val variance = values.map { (it - mean).pow(2) }.sum() / n

        for (lag in 1..minOf(maxLag, values.size - 1)) {
            var num = 0.0
            for (i in 0 until values.size - lag)
                num += (values[i] - mean) * (values[i + lag] - mean)
            val moran = if (variance == 0.0) 0.0 else
                (num / (n - lag)) / variance
            result["Moran_hydro_lag$lag"] = moran
        }
        return result
    }

    // ── Full descriptor set ──────────────────────────────────

    data class ProteinDescriptorSet(
        val aac: Map<String, Double>,
        val dpc: Map<String, Double>,
        val ctd: Map<String, Double>,
        val qso: Map<String, Double>,
        val paac: Map<String, Double>,
        val geary: Map<String, Double>,
        val moran: Map<String, Double>,
        val physicochemical: Map<String, Double>
    ) {
        fun toFlatMap(): Map<String, Double> {
            val all = mutableMapOf<String, Double>()
            aac.forEach { (k, v) -> all["AAC_$k"] = v }
            dpc.entries.take(100).forEach { (k, v) -> all["DPC_$k"] = v }
            ctd.forEach { (k, v) -> all["CTD_$k"] = v }
            qso.forEach { (k, v) -> all["QSO_$k"] = v }
            paac.forEach { (k, v) -> all["$k"] = v }
            geary.forEach { (k, v) -> all["$k"] = v }
            moran.forEach { (k, v) -> all["$k"] = v }
            physicochemical.forEach { (k, v) -> all["PC_$k"] = v }
            return all
        }
    }

    fun calculateAll(seq: ProteinSequence): ProteinDescriptorSet {
        return ProteinDescriptorSet(
            aac = aminoAcidComposition(seq),
            dpc = dipeptideComposition(seq),
            ctd = ctdDescriptors(seq),
            qso = quasiSequenceOrder(seq),
            paac = pseudoAAC(seq),
            geary = gearyAutocorrelation(seq),
            moran = moranAutocorrelation(seq),
            physicochemical = physicochemicalProfile(seq)
        )
    }

    // GRAVY (Grand Average of Hydropathicity — Kyte & Doolittle)
    fun gravyIndex(seq: ProteinSequence): Double {
        val kd = mapOf('A' to 1.8,'R' to -4.5,'N' to -3.5,'D' to -3.5,'C' to 2.5,'Q' to -3.5,
            'E' to -3.5,'G' to -0.4,'H' to -3.2,'I' to 4.5,'L' to 3.8,'K' to -3.9,'M' to 1.9,
            'F' to 2.8,'P' to -1.6,'S' to -0.8,'T' to -0.7,'V' to 4.2,'W' to -0.9,'Y' to -1.3)
        return seq.sequence.mapNotNull { kd[it] }.average()
    }

    // Aliphatic index (Ikai 1980)
    fun aliphaticIndex(seq: ProteinSequence): Double {
        if (seq.length == 0) return 0.0
        val a = seq.count('A').toDouble() / seq.length
        val v = seq.count('V').toDouble() / seq.length
        val i = seq.count('I').toDouble() / seq.length
        val l = seq.count('L').toDouble() / seq.length
        return 100 * (a + 2.9 * v + 3.9 * (i + l))
    }

    // Boman index (interaction potential with membranes)
    fun bomanIndex(seq: ProteinSequence): Double {
        val boman = mapOf('L' to -4.92,'I' to -4.44,'V' to -4.04,'F' to -2.98,'M' to -2.35,
            'W' to -2.33,'A' to -1.81,'C' to -1.28,'G' to -0.94,'Y' to 0.14,'T' to 2.57,
            'S' to 3.40,'H' to 4.66,'Q' to 5.54,'K' to 5.55,'N' to 6.64,'E' to 6.81,
            'D' to 8.72,'R' to 14.92,'P' to 0.0)
        return seq.sequence.mapNotNull { boman[it] }.average()
    }

    fun chargeAtPH(seq: ProteinSequence, pH: Double): Double {
        val pKa = mapOf('D' to 3.9,'E' to 4.1,'H' to 6.0,'C' to 8.3,'Y' to 10.1,'K' to 10.5,'R' to 12.5)
        var charge = 1.0 / (1.0 + 10.0.pow(pH - 8.0))   // N-terminus
        charge -= 1.0 / (1.0 + 10.0.pow(3.1 - pH))       // C-terminus
        for (aa in seq.sequence) {
            val pk = pKa[aa] ?: continue
            charge += when (aa) {
                'D','E','C','Y' -> -1.0 / (1.0 + 10.0.pow(pk - pH))
                'H','K','R'     ->  1.0 / (1.0 + 10.0.pow(pH - pk))
                else -> 0.0
            }
        }
        return charge
    }

    fun hydrophobicRatio(seq: ProteinSequence): Double {
        val hydrophobic = setOf('A','V','I','L','M','F','W','P')
        return seq.sequence.count { it in hydrophobic }.toDouble() / seq.length
    }

    fun extinctionCoefficient(seq: ProteinSequence): Double {
        // Pace 1995: ε = nW×5500 + nY×1490 + nC×125
        val nW = seq.count('W'); val nY = seq.count('Y'); val nC = seq.count('C')
        return nW * 5500.0 + nY * 1490.0 + nC * 125.0
    }
    fun physicochemicalProfile(seq: ProteinSequence): Map<String, Double> {
        return mapOf(
            "length"           to seq.length.toDouble(),
            "molecularWeight"  to seq.molecularWeight(),
            "isoelectricPoint" to seq.isoelectricPoint(),
            "aromaticity"      to seq.aromaticity(),
            "instabilityIndex" to seq.instabilityIndex(),
            "gravy"            to gravyIndex(seq),
            "aliphaticIndex"   to aliphaticIndex(seq),
            "bomanIndex"       to bomanIndex(seq),
            "chargeAtPH7"      to chargeAtPH(seq, 7.0),
            "hydrophobicRatio" to hydrophobicRatio(seq),
            "extinctionCoeff"  to extinctionCoefficient(seq)
        )
    }

}
