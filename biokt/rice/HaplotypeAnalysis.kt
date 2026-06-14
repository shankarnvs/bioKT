package biokt.rice

import biokt.*
import kotlin.math.*

// ============================================================
// HAPLOTYPE ANALYSIS — haplotype blocks, diversity, selection
// Subpackage: biokt.rice
//
// Gabriel haplotype blocks, pi/theta_W, windowed Tajima's D,
// haplotype frequency/diversity, and selection sweep scoring.
// ============================================================

// ─────────────────────────────────────────────────────────────
// HAPLOTYPE REPRESENTATION
// ─────────────────────────────────────────────────────────────

data class Haplotype(
    val id: String,
    val alleles: IntArray,   // 0=ref, 1=alt, -1=missing; one per variant in block
    val sampleName: String
) {
    val length: Int get() = alleles.size

    /** Pairwise difference count with another haplotype */
    fun differencesWith(other: Haplotype): Int {
        require(alleles.size == other.alleles.size)
        return alleles.indices.count { i ->
            alleles[i] >= 0 && other.alleles[i] >= 0 && alleles[i] != other.alleles[i]
        }
    }

    /** Fraction of sites that are missing */
    fun missingFraction(): Double = alleles.count { it < 0 }.toDouble() / alleles.size

    override fun toString(): String =
        "$sampleName: ${alleles.joinToString("") { if (it < 0) "?" else it.toString() }}"
}

data class HaplotypeBlock(
    val chrom: String,
    val start: Int,
    val end: Int,
    val snpPositions: List<Int>,
    val haplotypes: List<Haplotype>,
    val r2Matrix: Array<DoubleArray>   // LD matrix within block
) {
    val length: Int get() = end - start + 1
    val numSnps: Int get() = snpPositions.size
    val numSamples: Int get() = haplotypes.size

    /** Unique haplotype sequences and their frequencies */
    fun haplotypeFrequencies(): Map<String, Double> {
        val seqs = haplotypes.groupBy { it.alleles.joinToString("") { if (it < 0) "?" else it.toString() } }
        val total = haplotypes.size.toDouble()
        return seqs.mapValues { it.value.size / total }
    }

    /** Haplotype diversity (1 - Σpi²) */
    fun diversity(): Double {
        val freqs = haplotypeFrequencies().values
        return 1.0 - freqs.map { it * it }.sum()
    }

    /** Number of distinct haplotypes */
    fun numDistinctHaplotypes(): Int = haplotypeFrequencies().size

    override fun toString() =
        "HaploBlock[$chrom:$start-$end, ${numSnps} SNPs, ${numDistinctHaplotypes()} haplotypes]"
}

// ─────────────────────────────────────────────────────────────
// HAPLOTYPE BLOCK DEFINITION
// Gabriel et al. 2002 method
// ─────────────────────────────────────────────────────────────

object HaplotypeBlocks {

    /**
     * Define haplotype blocks using the Gabriel et al. 2002 rule.
     * A block is defined where 95% of SNP pairs show strong LD (D' ≥ 0.98)
     * and at most 5% show evidence of recombination (D' ≤ 0.70).
     *
     * @param vcf         VcfFile (must be biallelic SNPs)
     * @param chrom       chromosome to analyse
     * @param strongLD    minimum D' to count as "strong LD" (default 0.98)
     * @param recombBound maximum D' to count as "recombination" (default 0.70)
     * @param maxSpan     maximum block span in bp (default 200 kb)
     */
    fun defineBlocks(
        vcf: VcfFile,
        chrom: String,
        strongLD: Double = 0.98,
        recombBound: Double = 0.70,
        maxSpan: Int = 200_000
    ): List<HaplotypeBlock> {
        val snps = vcf.snps().biallelic().onChrom(chrom).records
            .sortedBy { it.pos }
        if (snps.size < 2) return emptyList()

        val blocks = mutableListOf<HaplotypeBlock>()
        var blockStart = 0

        var i = 0
        while (i < snps.size) {
            var blockEnd = i
            // Try to extend the block
            var j = i + 1
            while (j < snps.size && snps[j].pos - snps[i].pos <= maxSpan) {
                val ld = vcf.ldBetween(snps[i], snps[j])
                val dprime = ld.second

                // Check all pairs in candidate block
                var strongCount = 0; var recombCount = 0; var totalPairs = 0
                for (a in i..blockEnd) {
                    for (b in a + 1..j) {
                        val ldPair = vcf.ldBetween(snps[a], snps[b])
                        totalPairs++
                        if (ldPair.second >= strongLD)    strongCount++
                        if (ldPair.second <= recombBound) recombCount++
                    }
                }

                if (totalPairs > 0) {
                    val strongFrac  = strongCount.toDouble() / totalPairs
                    val recombFrac  = recombCount.toDouble() / totalPairs
                    if (strongFrac >= 0.95 && recombFrac <= 0.05) {
                        blockEnd = j
                    } else {
                        break
                    }
                }
                j++
            }

            if (blockEnd > i) {
                val blockSnps = snps.subList(i, blockEnd + 1)
                blocks.add(buildBlock(vcf, chrom, blockSnps))
                i = blockEnd + 1
            } else {
                i++
            }
        }
        return blocks
    }

    private fun buildBlock(vcf: VcfFile, chrom: String, snps: List<VcfRecord>): HaplotypeBlock {
        val positions = snps.map { it.pos }
        val start = positions.first()
        val end   = positions.last()

        // Build haplotypes from genotypes
        val haplotypes = vcf.samples.flatMap { sample ->
            val alleles0 = snps.map { rec ->
                val gt = rec.genotypes[sample]
                when {
                    gt == null || gt.isMissing() -> intArrayOf(-1, -1)
                    else -> {
                        val als = gt.alleles()
                        intArrayOf(
                            als.getOrNull(0)?.toIntOrNull() ?: -1,
                            als.getOrNull(1)?.toIntOrNull() ?: -1
                        )
                    }
                }
            }
            // Create two phased haplotypes per diploid sample
            listOf(
                Haplotype("${sample}_0", IntArray(snps.size) { i -> alleles0[i][0] }, sample),
                Haplotype("${sample}_1", IntArray(snps.size) { i -> alleles0[i][1] }, sample)
            )
        }

        // LD matrix for this block
        val n = snps.size
        val r2 = Array(n) { DoubleArray(n) { 0.0 } }
        for (a in 0 until n) {
            r2[a][a] = 1.0
            for (b in a + 1 until n) {
                val ld = vcf.ldBetween(snps[a], snps[b])
                r2[a][b] = ld.first; r2[b][a] = ld.first
            }
        }

        return HaplotypeBlock(chrom, start, end, positions, haplotypes, r2)
    }

    /** Block statistics summary */
    fun summarise(blocks: List<HaplotypeBlock>): String = buildString {
        append("=== Haplotype Blocks ===\n")
        append("Total blocks     : ${blocks.size}\n")
        append("Mean block size  : ${"%.1f".format(blocks.map { it.length }.average())} bp\n")
        append("Mean SNPs/block  : ${"%.1f".format(blocks.map { it.numSnps }.average())}\n")
        append("Mean haplotypes  : ${"%.1f".format(blocks.map { it.numDistinctHaplotypes() }.average())}\n")
        append("Mean diversity   : ${"%.4f".format(blocks.map { it.diversity() }.average())}\n")
        val totalBp = blocks.map { it.length }.sum()
        append("Total bp covered : $totalBp\n")
    }
}

// ─────────────────────────────────────────────────────────────
// WINDOWED DIVERSITY STATISTICS
// ─────────────────────────────────────────────────────────────

data class DiversityWindow(
    val chrom: String,
    val start: Int,
    val end: Int,
    val numSnps: Int,
    val pi: Double,              // nucleotide diversity
    val thetaW: Double,          // Watterson's theta
    val tajimasD: Double,        // Tajima's D
    val numHaplotypes: Int,
    val haplotypeDiv: Double     // haplotype diversity
)

object WindowedDiversity {

    /**
     * Compute windowed diversity statistics from a VCF.
     * @param vcf        VcfFile
     * @param chrom      chromosome name
     * @param windowSize window size in bp
     * @param stepSize   step between windows
     * @param nSeq       number of sequences/haplotypes (diploid: 2 × nSamples)
     */
    fun compute(
        vcf: VcfFile,
        chrom: String,
        windowSize: Int = 100_000,
        stepSize: Int = 50_000,
        nSeq: Int? = null
    ): List<DiversityWindow> {
        val n = nSeq ?: (vcf.samples.size * 2)
        val snps = vcf.snps().biallelic().onChrom(chrom).records.sortedBy { it.pos }
        if (snps.isEmpty()) return emptyList()

        val maxPos  = snps.maxBy { it.pos }!!.pos
        val results = mutableListOf<DiversityWindow>()
        var winStart = 1

        while (winStart <= maxPos) {
            val winEnd = winStart + windowSize - 1
            val winSnps = snps.filter { it.pos in winStart..winEnd }
            val s = winSnps.size   // segregating sites

            if (s >= 1) {
                // π (pairwise nucleotide diversity)
                var piSum = 0.0
                winSnps.forEach { rec ->
                    val p = rec.alleleFrequencies()[rec.ref] ?: 0.5
                    val q = 1.0 - p
                    piSum += 2 * p * q
                }
                val pi = piSum / windowSize

                // θ_W (Watterson)
                val an = (1 until n).map { 1.0 / it }.sum()
                val thetaW = if (an > 0) s.toDouble() / (an * windowSize) else 0.0

                // Tajima's D
                val tajimasD = computeTajimasD(winSnps, n)

                // Haplotype diversity in this window
                val blocks = mutableListOf<HaplotypeBlock>()
                val subVcf = VcfFile(vcf.header, winSnps)
                if (winSnps.size >= 2) {
                    val haploBlock = HaplotypeBlocks.defineBlocks(subVcf, chrom, maxSpan = windowSize)
                    val hapDiv = if (haploBlock.isNotEmpty())
                        haploBlock.map { it.diversity() }.average() else 0.0
                    val numH  = if (haploBlock.isNotEmpty())
                        haploBlock.map { it.numDistinctHaplotypes() }.max()!! ?: 0 else 0

                    results.add(DiversityWindow(chrom, winStart, winEnd, s,
                        pi, thetaW, tajimasD, numH, hapDiv))
                } else {
                    results.add(DiversityWindow(chrom, winStart, winEnd, s,
                        pi, thetaW, tajimasD, 0, 0.0))
                }
            } else {
                results.add(DiversityWindow(chrom, winStart, winEnd, 0, 0.0, 0.0, 0.0, 0, 0.0))
            }
            winStart += stepSize
        }
        return results
    }

    private fun computeTajimasD(snps: List<VcfRecord>, n: Int): Double {
        val s = snps.size.toDouble()
        if (s < 2 || n < 4) return 0.0

        val an = (1 until n).map { 1.0 / it }.sum()
        val bn = (1 until n).map { 1.0 / (it * it) }.sum()

        var piSum = 0.0
        snps.forEach { rec ->
            val p = rec.alleleFrequencies()[rec.ref] ?: 0.5
            val q = 1.0 - p
            piSum += 2 * p * q * n.toDouble() / (n - 1)
        }

        val theta_pi = piSum
        val theta_w  = s / an
        val diff     = theta_pi - theta_w

        // Variance of (π - θ_W) under neutrality
        val a1 = an; val a2 = bn
        val b1 = (n + 1).toDouble() / (3 * (n - 1))
        val b2 = 2.0 * (n * n + n + 3) / (9 * n * (n - 1))
        val c1 = b1 - 1.0 / a1
        val c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
        val e1 = c1 / a1
        val e2 = c2 / (a1 * a1 + a2)
        val variance = e1 * s + e2 * s * (s - 1)

        return if (variance > 0) diff / sqrt(variance) else 0.0
    }

    /** Find windows with extreme Tajima's D (potential selection sweeps) */
    fun selectionWindows(
        windows: List<DiversityWindow>,
        minSnps: Int = 5,
        zScoreThreshold: Double = 2.0
    ): List<DiversityWindow> {
        val valid = windows.filter { it.numSnps >= minSnps }
        if (valid.isEmpty()) return emptyList()
        val mean = valid.map { it.tajimasD }.average()
        val sd   = sqrt(valid.map { (it.tajimasD - mean).pow(2) }.average())
        return valid.filter { abs(it.tajimasD - mean) / sd.coerceAtLeast(0.001) >= zScoreThreshold }
            .sortedBy { it.tajimasD }  // most negative first (positive selection)
    }
}

// ─────────────────────────────────────────────────────────────
// XP-CLR — CROSS-POPULATION COMPOSITE LIKELIHOOD RATIO
// Simplified proxy for selection sweep detection between
// two populations (e.g. Indica vs Japonica)
// ─────────────────────────────────────────────────────────────

data class SelectionScore(
    val chrom: String,
    val start: Int,
    val end: Int,
    val score: Double,           // XP-CLR-like score (higher = stronger signal)
    val fstMean: Double,         // mean Fst in window
    val freqDiff: Double,        // mean allele frequency difference
    val numSnps: Int
) {
    override fun toString() = "SelectionSweep[$chrom:$start-$end score=${score.format(3)}]"
    private fun Double.format(n: Int) = "%.${n}f".format(this)
}

object SelectionSweepDetector {

    /**
     * Detect selection sweeps between two VCF populations.
     * Uses a composite score: Fst + allele frequency differentiation.
     * Higher scores indicate stronger signals of selection.
     */
    fun detect(
        pop1: VcfFile,
        pop2: VcfFile,
        chrom: String,
        windowSize: Int = 200_000,
        stepSize: Int = 50_000,
        minSnps: Int = 10
    ): List<SelectionScore> {
        // Get SNPs present in both populations on this chromosome
        val pos1 = pop1.snps().biallelic().onChrom(chrom).records
            .associateBy { "${it.chrom}:${it.pos}:${it.ref}:${it.alt.joinToString(",")}" }
        val pos2 = pop2.snps().biallelic().onChrom(chrom).records
            .associateBy { "${it.chrom}:${it.pos}:${it.ref}:${it.alt.joinToString(",")}" }

        val sharedKeys = pos1.keys.filter { it in pos2.keys }
        if (sharedKeys.size < minSnps) return emptyList()

        val shared = sharedKeys.mapNotNull { key ->
            val r1 = pos1[key] ?: return@mapNotNull null
            val r2 = pos2[key] ?: return@mapNotNull null
            Triple(r1.pos, r1.alleleFrequencies()[r1.ref] ?: 0.5,
                           r2.alleleFrequencies()[r2.ref] ?: 0.5)
        }.sortedBy { it.first }

        val maxPos  = shared.maxBy { it.first }?.first ?: return emptyList()
        val results = mutableListOf<SelectionScore>()
        var ws = 1

        while (ws <= maxPos) {
            val we    = ws + windowSize - 1
            val winSNPs = shared.filter { it.first in ws..we }
            if (winSNPs.size >= minSnps) {
                // Mean allele freq difference
                val freqDiff = winSNPs.map { abs(it.second - it.third) }.average()
                // Fst proxy: using allele freq differentiation
                val fstValues = winSNPs.map { (_, p1, p2) ->
                    val pMean = (p1 + p2) / 2
                    val qMean = 1 - pMean
                    if (pMean == 0.0 || qMean == 0.0) 0.0
                    else {
                        val ht = 2 * pMean * qMean
                        val hs = (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2
                        if (ht > 0) (ht - hs) / ht else 0.0
                    }
                }
                val fstMean = fstValues.average()
                // Composite score: Fst × freq_diff × log(nSNPs)
                val score = fstMean * freqDiff * ln(winSNPs.size.toDouble() + 1)
                results.add(SelectionScore(chrom, ws, we, score, fstMean, freqDiff, winSNPs.size))
            }
            ws += stepSize
        }

        return results.sortedByDescending { it.score }
    }

    /** Report top selection sweep candidates */
    fun topSweeps(scores: List<SelectionScore>, topN: Int = 20): List<SelectionScore> =
        scores.sortedByDescending { it.score }.take(topN)
}

// ─────────────────────────────────────────────────────────────
// HAPLOTYPE NETWORK BUILDER (simplified — star topology)
// ─────────────────────────────────────────────────────────────

data class HaploNetEdge(val from: String, val to: String, val mutations: Int, val samples: Int)

object HaplotypeNetwork {

    /** Build a median-joining proxy network from haplotype block */
    fun build(block: HaplotypeBlock): List<HaploNetEdge> {
        val freqs = block.haplotypeFrequencies()
        val sorted = freqs.entries.sortedByDescending { it.value }
        if (sorted.size < 2) return emptyList()

        val dominant = sorted.first().key
        val edges = mutableListOf<HaploNetEdge>()

        sorted.drop(1).forEach { (hapSeq, freq) ->
            val mutations = dominant.zip(hapSeq).count { (a, b) -> a != b && a != '?' && b != '?' }
            val sampCount = (freq * block.numSamples).toInt()
            edges.add(HaploNetEdge(dominant, hapSeq, mutations, sampCount))
        }
        return edges.sortedBy { it.mutations }
    }
}
