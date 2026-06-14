package biokt

import kotlin.math.*

// ============================================================
// SEQUENCE–SEQUENCE INTERACTIONS
// Protein-Protein, DNA-DNA, RNA-RNA, DNA-RNA interactions
// ============================================================

// ─── Interaction Score ───────────────────────────────────────

data class InteractionScore(
    val score: Double,
    val confidence: Double,       // 0-1
    val evidence: List<String>,   // supporting features
    val details: Map<String, Double>
) {
    val isLikelyInteracting: Boolean get() = confidence >= 0.5
    override fun toString() =
        "InteractionScore(score=%.3f, confidence=%.2f, likely=$isLikelyInteracting)".format(score, confidence)
}

// ─── Protein–Protein Interaction ────────────────────────────

object ProteinInteraction {

    // ── Amino-acid-based PPI features ────────────────────────

    // Chou & Shen (2006) — conjunction-trinucleotide composition for PPI
    fun conjunctionFeatures(seqA: ProteinSequence, seqB: ProteinSequence): Map<String, Double> {
        val result = mutableMapOf<String, Double>()
        val pdA = ProteinDescriptors.calculateAll(seqA)
        val pdB = ProteinDescriptors.calculateAll(seqB)

        // Pairwise feature products (key interaction signal)
        val aacA = pdA.aac; val aacB = pdB.aac
        for ((aa, va) in aacA) {
            val vb = aacB[aa] ?: 0.0
            result["AAC_prod_$aa"] = va * vb
            result["AAC_diff_$aa"] = abs(va - vb)
        }

        result["MW_ratio"]      = if (seqB.molecularWeight() == 0.0) 0.0 else seqA.molecularWeight() / seqB.molecularWeight()
        result["pI_diff"]       = abs(seqA.isoelectricPoint() - seqB.isoelectricPoint())
        result["GRAVY_sum"]     = ProteinDescriptors.gravyIndex(seqA) + ProteinDescriptors.gravyIndex(seqB)
        result["GRAVY_product"] = ProteinDescriptors.gravyIndex(seqA) * ProteinDescriptors.gravyIndex(seqB)
        result["charge_product"]= ProteinDescriptors.chargeAtPH(seqA, 7.0) * ProteinDescriptors.chargeAtPH(seqB, 7.0)
        result["length_ratio"]  = if (seqB.length == 0) 0.0 else seqA.length.toDouble() / seqB.length
        result["hydro_ratio"]   = if (ProteinDescriptors.hydrophobicRatio(seqB) == 0.0) 0.0
                                   else ProteinDescriptors.hydrophobicRatio(seqA) / ProteinDescriptors.hydrophobicRatio(seqB)
        return result
    }

    // ── Interface residue prediction ─────────────────────────
    // Propensity-based prediction of interface residues

    fun predictInterfaceResidues(seq: ProteinSequence): List<InterfaceResidue> {
        // Interface propensity scores (Ofran & Rost 2003)
        val propensity = mapOf(
            'A' to 0.87,'C' to 1.23,'D' to 0.96,'E' to 1.05,'F' to 1.16,'G' to 0.88,
            'H' to 1.29,'I' to 1.00,'K' to 1.06,'L' to 0.98,'M' to 1.33,'N' to 1.08,
            'P' to 0.65,'Q' to 1.10,'R' to 1.29,'S' to 0.96,'T' to 0.94,'V' to 0.88,
            'W' to 1.37,'Y' to 1.35
        )

        return seq.sequence.mapIndexed { i, aa ->
            val score = propensity[aa] ?: 1.0
            // Window-averaged score (context matters)
            val window = (maxOf(0, i - 3)..minOf(seq.length - 1, i + 3))
                .map { propensity[seq[it]] ?: 1.0 }.average()
            InterfaceResidue(
                position = i,
                aminoAcid = aa,
                propensityScore = window,
                isLikelyInterface = window >= 1.05
            )
        }
    }

    data class InterfaceResidue(
        val position: Int,
        val aminoAcid: Char,
        val propensityScore: Double,
        val isLikelyInterface: Boolean
    )

    // ── PPI scoring ──────────────────────────────────────────

    fun score(seqA: ProteinSequence, seqB: ProteinSequence): InteractionScore {
        val evidence = mutableListOf<String>()
        val details = mutableMapOf<String, Double>()

        // Charge complementarity (opposite charges favor interaction)
        val chargeA = ProteinDescriptors.chargeAtPH(seqA, 7.0)
        val chargeB = ProteinDescriptors.chargeAtPH(seqB, 7.0)
        val chargeComp = if (chargeA * chargeB < 0) 1.0 else 0.3
        details["charge_complementarity"] = chargeComp
        if (chargeComp > 0.7) evidence.add("Charge complementarity")

        // Hydrophobic patch complementarity
        val hydroA = ProteinDescriptors.hydrophobicRatio(seqA)
        val hydroB = ProteinDescriptors.hydrophobicRatio(seqB)
        val hydroScore = 1.0 - abs(hydroA - hydroB)
        details["hydrophobic_complementarity"] = hydroScore
        if (hydroScore > 0.7) evidence.add("Hydrophobic patch match")

        // GRAVY complementarity (one hydrophilic, one hydrophobic often interact)
        val gravyA = ProteinDescriptors.gravyIndex(seqA)
        val gravyB = ProteinDescriptors.gravyIndex(seqB)
        val gravyScore = if (gravyA * gravyB < 0) 0.8 else 0.4
        details["gravy_compatibility"] = gravyScore
        if (gravyScore > 0.7) evidence.add("Surface complementarity (GRAVY)")

        // Interface residue density
        val ifA = predictInterfaceResidues(seqA).count { it.isLikelyInterface }.toDouble() / seqA.length
        val ifB = predictInterfaceResidues(seqB).count { it.isLikelyInterface }.toDouble() / seqB.length
        val ifScore = (ifA + ifB) / 2.0
        details["interface_density"] = ifScore
        if (ifScore > 0.15) evidence.add("High interface propensity")

        // Size complementarity (roughly similar size proteins often interact)
        val sizeDiff = abs(log2(maxOf(seqA.length, 1).toDouble() / maxOf(seqB.length, 1).toDouble()))
        val sizeScore = exp(-sizeDiff / 3.0)
        details["size_compatibility"] = sizeScore

        // Dipeptide similarity (coevolution proxy)
        val dpcA = ProteinDescriptors.dipeptideComposition(seqA)
        val dpcB = ProteinDescriptors.dipeptideComposition(seqB)
        val keys = dpcA.keys.intersect(dpcB.keys)
        val dpcSim = if (keys.isEmpty()) 0.0 else
            keys.map { minOf(dpcA[it]!!, dpcB[it]!!) }.sum()
        details["dipeptide_similarity"] = dpcSim
        if (dpcSim > 0.3) evidence.add("Dipeptide composition similarity")

        val raw = chargeComp * 0.25 + hydroScore * 0.20 + gravyScore * 0.15 +
                  ifScore * 0.25 + sizeScore * 0.10 + dpcSim * 0.05

        val confidence = minOf(1.0, raw * 1.2)
        return InteractionScore(raw, confidence, evidence, details)
    }

    // ── Coiled-coil prediction (heptad repeat) ───────────────

    fun predictCoiledCoil(seq: ProteinSequence, window: Int = 21): List<CoiledCoilSegment> {
        val cc3Score = mapOf('A' to 1.44,'E' to 1.44,'I' to 1.40,'L' to 1.40,'M' to 1.36,
            'Q' to 1.33,'R' to 1.37,'S' to 0.99,'T' to 0.87,'V' to 1.31,
            'K' to 1.23,'N' to 1.00,'H' to 1.05,'D' to 0.98,'F' to 1.19,
            'W' to 1.08,'C' to 1.04,'Y' to 0.87,'G' to 0.52,'P' to 0.0)
        val segments = mutableListOf<CoiledCoilSegment>()
        val s = seq.sequence

        for (i in 0..s.length - window) {
            val win = s.substring(i, i + window)
            val score = win.mapNotNull { cc3Score[it] }.average()
            if (score >= 1.1) {
                segments.add(CoiledCoilSegment(i, i + window, score,
                    probability = minOf(1.0, (score - 1.0) * 3.0)))
            }
        }
        // Merge overlapping
        return mergeSegments(segments)
    }

    data class CoiledCoilSegment(val start: Int, val end: Int, val score: Double, val probability: Double)

    private fun mergeSegments(segs: List<CoiledCoilSegment>): List<CoiledCoilSegment> {
        if (segs.isEmpty()) return emptyList()
        val sorted = segs.sortedBy { it.start }
        val merged = mutableListOf(sorted.first())
        for (s in sorted.drop(1)) {
            val last = merged.last()
            if (s.start <= last.end) {
                merged[merged.size - 1] = last.copy(end = maxOf(last.end, s.end),
                    score = maxOf(last.score, s.score))
            } else merged.add(s)
        }
        return merged
    }

    // ── Transmembrane helix prediction (simple Kyte-Doolittle) ─

    fun predictTransmembraneHelices(seq: ProteinSequence, window: Int = 19, threshold: Double = 1.6): List<IntRange> {
        val kd = mapOf('A' to 1.8,'R' to -4.5,'N' to -3.5,'D' to -3.5,'C' to 2.5,'Q' to -3.5,
            'E' to -3.5,'G' to -0.4,'H' to -3.2,'I' to 4.5,'L' to 3.8,'K' to -3.9,'M' to 1.9,
            'F' to 2.8,'P' to -1.6,'S' to -0.8,'T' to -0.7,'V' to 4.2,'W' to -0.9,'Y' to -1.3)
        val s = seq.sequence
        val helices = mutableListOf<IntRange>()
        var i = 0
        while (i <= s.length - window) {
            val score = s.substring(i, i + window).mapNotNull { kd[it] }.average()
            if (score >= threshold) {
                val start = i
                while (i <= s.length - window &&
                    s.substring(i, i + window).mapNotNull { kd[it] }.average() >= threshold) i++
                helices.add(start..i + window)
            }
            i++
        }
        return helices
    }
}

// ─── DNA–DNA Interaction ─────────────────────────────────────

object DNAInteraction {

    // ── Hybridization thermodynamics (nearest-neighbor) ──────

    data class HybridizationResult(
        val deltaH: Double,    // kcal/mol
        val deltaS: Double,    // cal/mol·K
        val deltaG37: Double,  // kcal/mol at 37°C
        val tm: Double,        // melting temperature (°C)
        val numGC: Int,
        val numAT: Int,
        val mismatches: Int
    ) {
        override fun toString() =
            "Hybridization(ΔH=%.1f kcal/mol, ΔG37=%.1f kcal/mol, Tm=%.1f°C)".format(deltaH, deltaG37, tm)
    }

    // SantaLucia 1998 nearest-neighbor parameters
    private val NN_DH = mapOf(
        "AA" to -7.9, "AT" to -7.2, "TA" to -7.2, "CA" to -8.5,
        "GT" to -8.4, "CT" to -7.8, "GA" to -8.2, "CG" to -10.6,
        "GC" to -9.8, "GG" to -8.0, "TT" to -7.9, "TG" to -8.5,
        "AG" to -7.8, "TC" to -8.2, "AC" to -8.4, "CC" to -8.0
    )
    private val NN_DS = mapOf(
        "AA" to -22.2, "AT" to -20.4, "TA" to -21.3, "CA" to -22.7,
        "GT" to -22.4, "CT" to -21.0, "GA" to -22.2, "CG" to -27.2,
        "GC" to -24.4, "GG" to -19.9, "TT" to -22.2, "TG" to -22.7,
        "AG" to -21.0, "TC" to -22.2, "AC" to -22.4, "CC" to -19.9
    )

    fun hybridize(seq1: DNASequence, seq2: DNASequence,
                  saltConc: Double = 0.05,   // M NaCl
                  strandConc: Double = 250e-9 // 250 nM
    ): HybridizationResult {
        val s1 = seq1.sequence
        val s2 = seq2.reverseComplement().sequence
        val alignLen = minOf(s1.length, s2.length)

        var dH = 0.0; var dS = 0.0; var mismatches = 0

        for (i in 0 until alignLen - 1) {
            val dinuc = "${s1[i]}${s1[i + 1]}"
            dH += NN_DH[dinuc] ?: -8.0
            dS += NN_DS[dinuc] ?: -21.0
            if (i < s2.length && s1[i] != s2[i]) mismatches++
        }

        // Initiation parameters
        val initH  = if (s1.first() == 'G' || s1.first() == 'C') 0.1 else 2.3
        val initS  = if (s1.first() == 'G' || s1.first() == 'C') -2.8 else 4.1
        dH += initH * 2; dS += initS * 2

        // Salt correction (Owczarzy 2004)
        val saltCorr = 0.368 * (alignLen - 1) * ln(saltConc)
        dS += saltCorr

        val R = 1.987  // cal/mol·K
        val dG37 = dH - (310.15 * (dS / 1000.0))

        // Tm calculation
        val ct = if (mismatches == 0) strandConc else strandConc / 4
        val tm = (dH * 1000.0) / (dS + R * ln(ct)) - 273.15

        val gc = seq1.sequence.count { it == 'G' || it == 'C' }
        val at = seq1.sequence.count { it == 'A' || it == 'T' }

        return HybridizationResult(dH, dS, dG37, tm, gc, at, mismatches)
    }

    // ── Primer design ────────────────────────────────────────

    data class Primer(
        val sequence: DNASequence,
        val position: Int,
        val strand: Char,      // '+' or '-'
        val tm: Double,
        val gcContent: Double,
        val selfComplementarity: Double,
        val hairpinDeltaG: Double,
        val score: Double
    ) {
        override fun toString() = "Primer(${sequence.sequence}, pos=$position, Tm=${"%.1f".format(tm)}°C, GC=${"%.0f".format(gcContent)}%)"
    }

    fun designPrimers(
        template: DNASequence,
        targetStart: Int,
        targetEnd: Int,
        primerLen: IntRange = 18..25,
        tmRange: ClosedFloatingPointRange<Double> = 55.0..65.0,
        gcRange: ClosedFloatingPointRange<Double> = 40.0..60.0
    ): Pair<List<Primer>, List<Primer>> {
        val forward = mutableListOf<Primer>()
        val reverse = mutableListOf<Primer>()
        val seq = template.sequence

        // Forward primers (upstream of target)
        for (start in maxOf(0, targetStart - 50)..minOf(targetStart, seq.length - primerLen.first)) {
            for (len in primerLen) {
                if (start + len > seq.length) continue
                val pSeq = DNASequence(seq.substring(start, start + len))
                val primer = evaluatePrimer(pSeq, start, '+')
                if (primer.tm in tmRange && primer.gcContent in gcRange)
                    forward.add(primer)
            }
        }

        // Reverse primers (downstream of target)
        for (end in targetEnd..minOf(targetEnd + 50, seq.length)) {
            for (len in primerLen) {
                val start = end - len
                if (start < 0) continue
                val pSeq = DNASequence(seq.substring(start, end)).reverseComplement()
                val primer = evaluatePrimer(pSeq, start, '-')
                if (primer.tm in tmRange && primer.gcContent in gcRange)
                    reverse.add(primer)
            }
        }

        return forward.sortedByDescending { it.score } to reverse.sortedByDescending { it.score }
    }

    private fun evaluatePrimer(seq: DNASequence, pos: Int, strand: Char): Primer {
        val gc = seq.gcContent()
        val tm = seq.meltingTemperature()

        // Self-complementarity score (lower is better)
        val comp = seq.complement().sequence.reversed()
        var maxSelfCompInt = 0
        for (offset in 0 until seq.length) {
            val matches = (0 until minOf(seq.length - offset, comp.length)).count { i ->
                seq.sequence[i + offset] == comp[i]
            }
            if (matches > maxSelfCompInt) maxSelfCompInt = matches
        }
        val maxSelfComp = maxSelfCompInt.toDouble() / seq.length

        // Hairpin deltaG (simplified)
        val hairpinDG = -1.0 * seq.length * 0.1 + maxSelfComp * 2.0

        val score = 1.0 - abs(gc - 50.0) / 50.0 - abs(tm - 60.0) / 20.0 -
                    maxSelfComp * 0.5

        return Primer(seq, pos, strand, tm, gc, maxSelfComp, hairpinDG, score)
    }

    // ── DNA–DNA binding free energy ──────────────────────────

    fun bindingAffinity(
        seq1: DNASequence,
        seq2: DNASequence,
        topology: String = "linear"  // "linear" or "circular"
    ): Double {
        val hyb = hybridize(seq1, seq2)
        val topCorrection = if (topology == "circular") -3.5 else 0.0  // kcal/mol
        return hyb.deltaG37 + topCorrection
    }

    // ── FRET efficiency (for DNA-FRET experiments) ───────────

    fun fretEfficiency(r: Double, r0: Double): Double {
        // Förster equation: E = 1/(1 + (r/r0)^6)
        return 1.0 / (1.0 + (r / r0).pow(6))
    }
}

// ─── RNA–RNA Interaction (secondary structure aware) ─────────

object RNAInteraction {

    // ── RNA secondary structure (Zuker-style energy minimization) ─
    // Simplified Nussinov algorithm for base pair maximization

    data class RNAStructure(
        val sequence: RNASequence,
        val dotBracket: String,
        val numBasePairs: Int,
        val freeEnergy: Double
    ) {
        fun stems(): List<IntRange> {
            val stems = mutableListOf<IntRange>()
            var i = 0
            while (i < dotBracket.length) {
                if (dotBracket[i] == '(') {
                    val start = i
                    while (i < dotBracket.length && dotBracket[i] == '(') i++
                    stems.add(start until i)
                }
                i++
            }
            return stems
        }
    }

    fun fold(seq: RNASequence): RNAStructure {
        val s = seq.sequence
        val n = s.length

        // Nussinov DP for max base pairs
        val dp = Array(n) { IntArray(n) }

        fun canPair(i: Int, j: Int): Boolean {
            if (j - i < 3) return false
            val pairs = setOf("AU", "UA", "GC", "CG", "GU", "UG")
            return "${s[i]}${s[j]}" in pairs
        }

        for (len in 1 until n) {
            for (i in 0 until n - len) {
                val j = i + len
                dp[i][j] = dp[i][j - 1]
                if (canPair(i, j)) dp[i][j] = maxOf(dp[i][j], dp[i + 1][j - 1] + 1)
                for (k in i + 1 until j) {
                    dp[i][j] = maxOf(dp[i][j], dp[i][k] + dp[k + 1][j])
                }
            }
        }

        // Traceback
        val structure = CharArray(n) { '.' }
        fun traceback(i: Int, j: Int) {
            if (i >= j) return
            if (dp[i][j] == dp[i][j - 1]) { traceback(i, j - 1) }
            else if (canPair(i, j) && dp[i][j] == dp[i + 1][j - 1] + 1) {
                structure[i] = '('; structure[j] = ')'
                traceback(i + 1, j - 1)
            } else {
                for (k in i + 1 until j) {
                    if (dp[i][j] == dp[i][k] + dp[k + 1][j]) {
                        traceback(i, k); traceback(k + 1, j); return
                    }
                }
            }
        }
        traceback(0, n - 1)
        val dotBracket = String(structure)
        val numPairs = dotBracket.count { it == '(' }

        // Approximate free energy: -1.5 kcal per base pair (simplified)
        val freeEnergy = numPairs * -1.5 + (n - numPairs * 2) * 0.1

        return RNAStructure(seq, dotBracket, numPairs, freeEnergy)
    }

    // ── RNA–RNA kissing loop interaction ─────────────────────

    fun kissLoop(struct1: RNAStructure, struct2: RNAStructure): Double {
        val loop1 = extractLoops(struct1)
        val loop2 = extractLoops(struct2)

        var bestScore = 0.0
        for (l1 in loop1) for (l2 in loop2) {
            val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.LOCAL)
            val r1 = RNASequence(l1)
            val r2 = r1.complement()
            val aln = aligner.align(l1, r2.sequence)
            if (aln.score > bestScore) bestScore = aln.score
        }
        return bestScore
    }

    private fun extractLoops(struct: RNAStructure): List<String> {
        val loops = mutableListOf<String>()
        val s = struct.sequence.sequence
        val db = struct.dotBracket
        var i = 0
        while (i < db.length) {
            if (db[i] == '.') {
                val start = i
                while (i < db.length && db[i] == '.') i++
                if (i - start >= 3) loops.add(s.substring(start, i))
            } else i++
        }
        return loops
    }

    // ── RNA–RNA intermolecular base pairing ──────────────────

    fun hybridize(seq1: RNASequence, seq2: RNASequence): Double {
        val comp = seq2.complement()
        val aligner = PairwiseAligner(
            mode = PairwiseAligner.Mode.LOCAL,
            matchScore = 3.0, mismatchScore = -2.0,
            gapOpen = -5.0, gapExtend = -1.0
        )
        val aln = aligner.align(seq1.sequence, comp.sequence)
        // Convert score to approximate free energy
        return aln.score * -0.7
    }

    // ── miRNA target prediction ──────────────────────────────

    data class MiRNATarget(
        val miRNA: RNASequence,
        val targetMRNA: RNASequence,
        val seedMatchPosition: Int,
        val seedLength: Int,
        val hybridizationScore: Double,
        val freeEnergy: Double,
        val siteType: String  // "8mer", "7mer-m8", "7mer-A1", "6mer"
    )

    fun predictMiRNATargets(
        miRNA: RNASequence,
        mRNA: RNASequence,
        minSeedLength: Int = 6
    ): List<MiRNATarget> {
        val targets = mutableListOf<MiRNATarget>()
        val seed = miRNA.sequence.substring(1, 8)  // positions 2-8 (seed region)
        val seedComplement = RNASequence(seed).complement().sequence

        // Search for seed matches in 3' UTR proxy (full mRNA here)
        var pos = 0
        while (pos < mRNA.length - seed.length) {
            val window = mRNA.sequence.substring(pos, minOf(pos + seed.length + 2, mRNA.length))
            val matches = window.indices.count { i ->
                i < seed.length && window[i] == seedComplement.getOrElse(i) { ' ' }
            }
            if (matches >= minSeedLength) {
                val hybScore = matches.toDouble() / seed.length
                val siteType = when {
                    matches == 8 -> "8mer"
                    matches == 7 && window.getOrElse(0) { 'X' } == 'A' -> "7mer-A1"
                    matches == 7 -> "7mer-m8"
                    else -> "6mer"
                }
                targets.add(MiRNATarget(
                    miRNA, mRNA, pos, matches, hybScore,
                    freeEnergy = hybScore * -3.5, siteType = siteType
                ))
            }
            pos++
        }
        return targets.sortedByDescending { it.hybridizationScore }
    }
}

// ─── DNA–RNA Interaction ─────────────────────────────────────

object DNARNAInteraction {

    // R-loop formation propensity
    fun rLoopPropensity(dna: DNASequence, rna: RNASequence): Double {
        val dnaStr = dna.sequence
        val rnaStr = rna.backTranscribe().sequence  // convert to DNA for comparison

        val aligner = PairwiseAligner(mode = PairwiseAligner.Mode.LOCAL,
            matchScore = 2.0, mismatchScore = -1.0)
        val aln = aligner.align(dnaStr, rnaStr)

        val gcContent = dna.gcContent()
        // R-loops favor GC-rich sequences
        return aln.identity * (gcContent / 100.0) * aln.score / maxOf(dna.length, rna.length)
    }

    // Transcription bubble stability
    fun transcriptionBubbleStability(dna: DNASequence, windowSize: Int = 10): List<Pair<Int, Double>> {
        val seq = dna.sequence
        val results = mutableListOf<Pair<Int, Double>>()
        for (i in 0..seq.length - windowSize) {
            val window = seq.substring(i, i + windowSize)
            val at = window.count { it == 'A' || it == 'T' }.toDouble() / windowSize
            // AT-rich windows are less stable (easier to melt = easier transcription)
            results.add(i to at)
        }
        return results
    }

    // CRISPR guide RNA efficiency prediction
    fun crisprGuideEfficiency(guide: RNASequence, target: DNASequence): Double {
        if (guide.length < 20) return 0.0
        val guideSeq = guide.backTranscribe().sequence.takeLast(20)
        val targetSeq = target.sequence

        // Check PAM (NGG for SpCas9)
        val hasPam = targetSeq.endsWith("GG")
        if (!hasPam) return 0.0

        // GC content of guide (optimal 40-70%)
        val gc = DNASequence(guideSeq).gcContent()
        val gcScore = if (gc in 40.0..70.0) 1.0 else 0.5

        // Position-specific nucleotide preferences (simplified Doench 2016)
        val scoreMap = mapOf(0 to 0.1, 1 to 0.1, 19 to 0.2, 18 to 0.2, 17 to 0.15)
        var posScore = 1.0
        for ((pos, bonus) in scoreMap) {
            if (pos < guideSeq.length && guideSeq[pos] == 'G') posScore += bonus
        }

        // Seed region mismatches (positions 1-12 from PAM) matter most
        val seed = guideSeq.takeLast(12)
        val targetSeed = targetSeq.dropLast(3).takeLast(12)
        val mismatches = seed.zip(targetSeed).count { (a, b) -> a != b }
        val mismatchPenalty = mismatches * 0.2

        return maxOf(0.0, minOf(1.0, gcScore * posScore - mismatchPenalty))
    }
}
