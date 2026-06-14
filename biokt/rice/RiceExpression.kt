package biokt.rice

import biokt.*
import kotlin.math.*

// ============================================================
// RICE EXPRESSION — RNA-seq normalisation and expression analysis
// Subpackage: biokt.rice
//
// TPM / FPKM / RPKM normalisation, differential expression,
// co-expression, RiceXPro tissue atlas, rice stress gene sets.
// ============================================================

// ─────────────────────────────────────────────────────────────
// EXPRESSION MATRIX
// ─────────────────────────────────────────────────────────────

class ExpressionMatrix(
    val geneIds: List<String>,
    val sampleIds: List<String>,
    private val data: Array<DoubleArray>   // [geneIdx][sampleIdx]
) {
    val numGenes: Int   get() = geneIds.size
    val numSamples: Int get() = sampleIds.size

    operator fun get(geneIdx: Int, sampleIdx: Int): Double = data[geneIdx][sampleIdx]
    operator fun get(geneId: String, sampleId: String): Double? {
        val gi = geneIds.indexOf(geneId).takeIf { it >= 0 } ?: return null
        val si = sampleIds.indexOf(sampleId).takeIf { it >= 0 } ?: return null
        return data[gi][si]
    }

    fun row(geneIdx: Int): DoubleArray = data[geneIdx].clone()
    fun row(geneId: String): DoubleArray? {
        val i = geneIds.indexOf(geneId).takeIf { it >= 0 } ?: return null
        return data[i].clone()
    }

    fun column(sampleIdx: Int): DoubleArray = DoubleArray(numGenes) { data[it][sampleIdx] }
    fun column(sampleId: String): DoubleArray? {
        val i = sampleIds.indexOf(sampleId).takeIf { it >= 0 } ?: return null
        return DoubleArray(numGenes) { data[it][i] }
    }

    fun submatrix(genes: List<String>): ExpressionMatrix {
        val indices = genes.mapNotNull { g -> geneIds.indexOf(g).takeIf { it >= 0 } }
        val newData = Array(indices.size) { data[indices[it]] }
        return ExpressionMatrix(indices.map { geneIds[it] }, sampleIds, newData)
    }

    fun log2transform(pseudocount: Double = 1.0): ExpressionMatrix {
        val newData = Array(numGenes) { i ->
            DoubleArray(numSamples) { j -> log2(data[i][j] + pseudocount) }
        }
        return ExpressionMatrix(geneIds, sampleIds, newData)
    }

    fun toCsv(): String = buildString {
        append("GeneID,${sampleIds.joinToString(",")}\n")
        geneIds.forEachIndexed { i, gId ->
            append("$gId,${data[i].joinToString(",") { "%.4f".format(it) }}\n")
        }
    }

    companion object {
        fun fromCsv(csv: String): ExpressionMatrix {
            val lines  = csv.trim().lines()
            val header = lines[0].split(",")
            val samples = header.drop(1)
            val genes  = mutableListOf<String>()
            val rows   = mutableListOf<DoubleArray>()
            lines.drop(1).forEach { line ->
                val cols = line.split(",")
                if (cols.size > 1) {
                    genes.add(cols[0])
                    rows.add(DoubleArray(samples.size) { i ->
                        cols.getOrElse(i + 1) { "0" }.toDoubleOrNull() ?: 0.0
                    })
                }
            }
            return ExpressionMatrix(genes, samples, rows.toTypedArray())
        }
    }
}

// ─────────────────────────────────────────────────────────────
// NORMALISATION
// ─────────────────────────────────────────────────────────────

object ExpressionNormaliser {

    /**
     * Normalise raw read counts to TPM (Transcripts Per Million).
     * @param counts   raw count matrix [gene × sample]
     * @param lengths  gene lengths in base pairs (same order as genes)
     */
    fun toTPM(counts: Array<DoubleArray>, lengths: DoubleArray): Array<DoubleArray> {
        require(counts.size == lengths.size) { "counts and lengths must have same number of genes" }
        val numSamples = counts[0].size
        val rpk = Array(counts.size) { i ->
            DoubleArray(numSamples) { j -> counts[i][j] / (lengths[i] / 1000.0) }
        }
        val scalingFactors = DoubleArray(numSamples) { j ->
            rpk.map { row -> row[j] }.sum() / 1_000_000.0
        }
        return Array(counts.size) { i ->
            DoubleArray(numSamples) { j ->
                if (scalingFactors[j] > 0) rpk[i][j] / scalingFactors[j] else 0.0
            }
        }
    }

    /**
     * Normalise to FPKM (Fragments Per Kilobase per Million mapped reads).
     * @param counts     raw count matrix
     * @param lengths    gene lengths in bp
     * @param totalReads total mapped reads per sample (if null, computed from counts)
     */
    fun toFPKM(
        counts: Array<DoubleArray>,
        lengths: DoubleArray,
        totalReads: DoubleArray? = null
    ): Array<DoubleArray> {
        val numSamples = counts[0].size
        val totals = totalReads ?: DoubleArray(numSamples) { j -> counts.map { it[j] }.sum() }
        return Array(counts.size) { i ->
            DoubleArray(numSamples) { j ->
                if (totals[j] > 0 && lengths[i] > 0)
                    (counts[i][j] * 1e9) / (lengths[i] * totals[j])
                else 0.0
            }
        }
    }

    /**
     * Normalise to RPKM (Reads Per Kilobase per Million).
     * Functionally identical to FPKM for single-end reads.
     */
    fun toRPKM(
        counts: Array<DoubleArray>,
        lengths: DoubleArray,
        totalReads: DoubleArray? = null
    ): Array<DoubleArray> = toFPKM(counts, lengths, totalReads)

    /**
     * DESeq2-style median-of-ratios normalisation (size factor normalisation).
     * Robust to outliers. Recommended for DE analysis.
     */
    fun sizeFactorNormalise(counts: Array<DoubleArray>): Pair<Array<DoubleArray>, DoubleArray> {
        val numGenes   = counts.size
        val numSamples = counts[0].size

        // Geometric mean per gene (skip zeros)
        val geomMeans = DoubleArray(numGenes) { i ->
            val nonZero = counts[i].filter { it > 0 }
            if (nonZero.isEmpty()) 0.0
            else exp(nonZero.map { ln(it) }.sum() / nonZero.size)
        }

        // Ratios for each sample
        val sizeFactors = DoubleArray(numSamples) { j ->
            val ratios = (0 until numGenes)
                .filter { i -> geomMeans[i] > 0 && counts[i][j] > 0 }
                .map { i -> counts[i][j] / geomMeans[i] }
            if (ratios.isEmpty()) 1.0 else {
                val sorted = ratios.sorted()
                sorted[sorted.size / 2]  // median
            }
        }

        val normalised = Array(numGenes) { i ->
            DoubleArray(numSamples) { j ->
                if (sizeFactors[j] > 0) counts[i][j] / sizeFactors[j] else 0.0
            }
        }
        return normalised to sizeFactors
    }

    /** Quantile normalise — make all samples have the same distribution */
    fun quantileNormalise(matrix: Array<DoubleArray>): Array<DoubleArray> {
        val numGenes   = matrix.size
        val numSamples = matrix[0].size
        // Sort each column, compute mean ranks, then substitute back
        val sorted = Array(numSamples) { j ->
            matrix.map { it[j] }.sorted()
        }
        val rankMeans = DoubleArray(numGenes) { rank ->
            sorted.map { it[rank] }.average()
        }
        return Array(numGenes) { i ->
            DoubleArray(numSamples) { j ->
                val rank = matrix.map { it[j] }.sortedDescending().indexOf(matrix[i][j])
                rankMeans[maxOf(0, numGenes - 1 - rank)]
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
// DIFFERENTIAL EXPRESSION
// ─────────────────────────────────────────────────────────────

data class DEResult(
    val geneId: String,
    val log2FoldChange: Double,
    val meanExprCase: Double,
    val meanExprControl: Double,
    val pValue: Double,
    val adjustedPValue: Double,
    val zScore: Double,
    val isSignificant: Boolean
) {
    val isUpregulated: Boolean   get() = isSignificant && log2FoldChange > 0
    val isDownregulated: Boolean get() = isSignificant && log2FoldChange < 0
    override fun toString() = "$geneId  log2FC=${"%.3f".format(log2FoldChange)}  padj=${"%.4f".format(adjustedPValue)}"
}

object DifferentialExpression {

    /**
     * Compute differential expression between two groups of samples.
     * Uses Welch's t-test + Benjamini-Hochberg FDR correction.
     *
     * @param matrix   ExpressionMatrix (log2-transformed recommended)
     * @param caseIds  column IDs of case/treatment samples
     * @param ctrlIds  column IDs of control samples
     * @param minExpr  minimum mean expression to consider a gene
     */
    fun analyse(
        matrix: ExpressionMatrix,
        caseIds: List<String>,
        ctrlIds: List<String>,
        minExpr: Double = 1.0,
        log2fcThreshold: Double = 1.0,
        padjThreshold: Double = 0.05
    ): List<DEResult> {
        val caseIdx = caseIds.mapNotNull { id -> matrix.sampleIds.indexOf(id).takeIf { it >= 0 } }
        val ctrlIdx = ctrlIds.mapNotNull { id -> matrix.sampleIds.indexOf(id).takeIf { it >= 0 } }
        require(caseIdx.isNotEmpty() && ctrlIdx.isNotEmpty()) { "No valid sample indices found" }

        val results = mutableListOf<DEResult>()
        matrix.geneIds.forEachIndexed { i, geneId ->
            val row    = matrix.row(i)
            val case   = caseIdx.map { row[it] }
            val ctrl   = ctrlIdx.map { row[it] }
            val meanC  = case.average()
            val meanK  = ctrl.average()
            if (meanC < minExpr && meanK < minExpr) return@forEachIndexed

            val log2FC = meanC - meanK  // assumes log2-transformed input
            val (tStat, pVal) = welchTTest(case, ctrl)
            val zScore = if (pVal > 0) -sign(log2FC) * log10(pVal).coerceAtLeast(-300.0) else 0.0
            results.add(DEResult(geneId, log2FC, meanC, meanK, pVal, pVal, zScore, false))
        }

        // Benjamini-Hochberg FDR adjustment
        val sorted   = results.sortedBy { it.pValue }
        val n        = sorted.size
        val adjusted = sorted.mapIndexed { rank, r ->
            val adjP = (r.pValue * n / (rank + 1)).coerceAtMost(1.0)
            r.copy(adjustedPValue = adjP,
                isSignificant = adjP <= padjThreshold && abs(r.log2FoldChange) >= log2fcThreshold)
        }
        return adjusted.sortedBy { it.adjustedPValue }
    }

    /** Welch's t-test for two independent samples */
    fun welchTTest(a: List<Double>, b: List<Double>): Pair<Double, Double> {
        if (a.size < 2 || b.size < 2) return 0.0 to 1.0
        val ma = a.average(); val mb = b.average()
        val va = a.map { (it - ma).pow(2) }.sum() / (a.size - 1)
        val vb = b.map { (it - mb).pow(2) }.sum() / (b.size - 1)
        val se = sqrt(va / a.size + vb / b.size)
        if (se == 0.0) return 0.0 to 1.0
        val t = (ma - mb) / se
        // Welch-Satterthwaite degrees of freedom
        val df = (va / a.size + vb / b.size).pow(2) /
                 ((va / a.size).pow(2) / (a.size - 1) + (vb / b.size).pow(2) / (b.size - 1))
        val p  = tDistPValue(t, df)
        return t to p
    }

    /** Approximate p-value from t-distribution (two-tailed) */
    private fun tDistPValue(t: Double, df: Double): Double {
        val x = df / (df + t * t)
        val ibeta = regularisedIncompleteBeta(df / 2.0, 0.5, x)
        return ibeta.coerceIn(1e-300, 1.0)
    }

    /** Regularised incomplete beta function approximation */
    private fun regularisedIncompleteBeta(a: Double, b: Double, x: Double): Double {
        if (x <= 0.0) return 0.0; if (x >= 1.0) return 1.0
        val lbeta = lgamma(a + b) - lgamma(a) - lgamma(b)
        val logX  = ln(x); val log1X = ln(1 - x)
        var sum = 0.0; var term = 1.0
        for (k in 0..100) {
            if (k > 0) term *= (a + k - 1) * x / k
            sum += term / (a + k)
            if (abs(term / (a + k)) < 1e-10) break
        }
        return exp(a * logX + b * log1X - lbeta + ln(sum))
    }

    private fun lgamma(x: Double): Double {
        // Lanczos approximation
        val g = 7.0
        val c = doubleArrayOf(0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
        if (x < 0.5) return ln(PI / sin(PI * x)) - lgamma(1 - x)
        var xv = x - 1; var t = xv + g + 0.5; var s = c[0]
        for (i in 1..8) s += c[i] / (xv + i)
        return 0.5 * ln(2 * PI) + (xv + 0.5) * ln(t) - t + ln(s)
    }
}

// ─────────────────────────────────────────────────────────────
// CO-EXPRESSION
// ─────────────────────────────────────────────────────────────

data class CoexprEdge(val geneA: String, val geneB: String, val correlation: Double) {
    override fun toString() = "$geneA — $geneB: ${"%.4f".format(correlation)}"
}

object CoexpressionAnalysis {

    /** Compute Pearson correlation matrix for all gene pairs */
    fun correlationMatrix(matrix: ExpressionMatrix): Array<DoubleArray> {
        val n = matrix.numGenes
        val corr = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            corr[i][i] = 1.0
            for (j in i + 1 until n) {
                val r = pearsonR(matrix.row(i), matrix.row(j))
                corr[i][j] = r; corr[j][i] = r
            }
        }
        return corr
    }

    /** Get top co-expressed gene pairs above threshold */
    fun topEdges(
        matrix: ExpressionMatrix,
        threshold: Double = 0.8,
        topN: Int = 100
    ): List<CoexprEdge> {
        val edges = mutableListOf<CoexprEdge>()
        for (i in 0 until matrix.numGenes) {
            for (j in i + 1 until matrix.numGenes) {
                val r = pearsonR(matrix.row(i), matrix.row(j))
                if (abs(r) >= threshold)
                    edges.add(CoexprEdge(matrix.geneIds[i], matrix.geneIds[j], r))
            }
        }
        return edges.sortedByDescending { abs(it.correlation) }.take(topN)
    }

    /** Find genes correlated with a query gene above threshold */
    fun coexpressedWith(
        matrix: ExpressionMatrix,
        geneId: String,
        threshold: Double = 0.7,
        topN: Int = 20
    ): List<Pair<String, Double>> {
        val row = matrix.row(geneId) ?: return emptyList()
        return matrix.geneIds.mapIndexedNotNull { i, id ->
            if (id == geneId) null
            else {
                val r = pearsonR(row, matrix.row(i))
                if (abs(r) >= threshold) id to r else null
            }
        }.sortedByDescending { abs(it.second) }.take(topN)
    }

    fun pearsonR(x: DoubleArray, y: DoubleArray): Double {
        if (x.size != y.size || x.size < 3) return 0.0
        val mx = x.average(); val my = y.average()
        var num = 0.0; var dx = 0.0; var dy = 0.0
        for (i in x.indices) {
            val ax = x[i] - mx; val ay = y[i] - my
            num += ax * ay; dx += ax * ax; dy += ay * ay
        }
        return if (dx == 0.0 || dy == 0.0) 0.0 else num / sqrt(dx * dy)
    }
}

// ─────────────────────────────────────────────────────────────
// RICE STRESS GENE SETS
// Curated from published literature and RiceXPro
// ─────────────────────────────────────────────────────────────

object RiceStressGeneSets {

    data class GeneSet(val name: String, val description: String, val genes: List<String>)

    // RAP-DB IDs for key stress-responsive genes
    val DROUGHT = GeneSet("Drought response",
        "Genes up-regulated under drought/water deficit stress",
        listOf(
            "Os01g0843200",  // SNAC1/OsNAC9
            "Os06g0127800",  // OsDREB1A
            "Os06g0194800",  // OsDREB2A
            "Os11g0639100",  // OsDREB1B
            "Os09g0522200",  // OsRD29A
            "Os05g0329700",  // OsAP37
            "Os01g0841800",  // OsNAC45
            "Os07g0154100",  // OsWRKY45
            "Os04g0448400",  // OsER1
            "Os02g0818000"   // OsNAC063
        )
    )

    val SALINITY = GeneSet("Salinity tolerance",
        "Genes associated with salt stress response",
        listOf(
            "Os01g0307700",  // OsSOS1
            "Os12g0476200",  // OsHKT1;5
            "Os04g0448400",  // OsNHX1
            "Os03g0401200",  // OsMPK5
            "Os07g0154100",  // OsWRKY45
            "Os05g0479700",  // OsGSK1
            "Os06g0682500",  // OsMYB2
            "Os02g0818000",  // OsNAC063
            "Os11g0574700",  // OsCPS2
            "Os11g0547800"   // OsNAC45
        )
    )

    val COLD = GeneSet("Cold stress",
        "Genes responsive to low temperature / cold acclimation",
        listOf(
            "Os06g0127800",  // OsDREB1A
            "Os06g0194800",  // OsDREB1B — cold-inducible
            "Os03g0236100",  // OsCBF1/DREB1B
            "Os09g0315700",  // OsCBF3
            "Os02g0585500",  // OsICE1
            "Os12g0422900",  // OsMYBS3
            "Os04g0448900",  // OsEIN3
            "Os03g0765300",  // OsSOC1
            "Os01g0884100",  // OsHSP90
            "Os09g0447800"   // OsERF1
        )
    )

    val BLAST = GeneSet("Blast resistance",
        "Genes involved in response to Magnaporthe oryzae",
        listOf(
            "Os11g0508800",  // Pi-ta
            "Os11g0695900",  // Pi-d2
            "Os12g0476200",  // Pi9
            "Os06g0140700",  // Pi2/Pi-z
            "Os01g0159500",  // Pit
            "Os11g0616200",  // Pi37
            "Os11g0459800",  // Pik-m
            "Os11g0484200",  // Pikp
            "Os11g0603300",  // Pib
            "Os07g0154100",  // OsWRKY45
            "Os04g0576800",  // OsNPR1
            "Os12g0441800"   // OsPR1a
        )
    )

    val SUBMERGENCE = GeneSet("Submergence tolerance",
        "Genes responsive to flooding / submergence",
        listOf(
            "Os09g0449500",  // Sub1A
            "Os09g0449800",  // Sub1B
            "Os09g0450200",  // Sub1C
            "Os09g0432700",  // OsERF2
            "Os06g0568700",  // OsADH1
            "Os01g0884300",  // OsSUS3
            "Os09g0432700",  // OsERF3
            "Os04g0448400",  // OsERF1
            "Os03g0829100"   // OsPDC2
        )
    )

    val HEAT = GeneSet("Heat stress",
        "Genes up-regulated under high temperature",
        listOf(
            "Os01g0884100",  // OsHSP90
            "Os03g0288500",  // OsHSP101
            "Os05g0402800",  // OsHsfA2a
            "Os03g0322600",  // OsHsfB2b
            "Os02g0558300",  // OsHSP70
            "Os05g0569100",  // OsHSP17.9
            "Os01g0633100",  // OsHSP18.0
            "Os07g0668100"   // OsGrpE
        )
    )

    val ALL_SETS = listOf(DROUGHT, SALINITY, COLD, BLAST, SUBMERGENCE, HEAT)

    /** Find which stress gene sets contain a given RAP-DB ID */
    fun setsContaining(rapId: String): List<GeneSet> =
        ALL_SETS.filter { rapId in it.genes }

    /** Overlap a list of DE gene IDs with stress gene sets */
    fun enrichment(
        deGenes: List<String>,
        background: Int = 35_000   // approximate rice gene count
    ): List<EnrichmentResult> {
        return ALL_SETS.map { gs ->
            val overlap = deGenes.count { it in gs.genes }
            val setSize = gs.genes.size
            val listSize = deGenes.size
            // Fisher's exact test approximation (hypergeometric)
            val expected = listSize.toDouble() * setSize / background
            val enrichmentRatio = if (expected > 0) overlap / expected else 0.0
            val pVal = hypergeometricPValue(background, setSize, listSize, overlap)
            EnrichmentResult(gs.name, overlap, setSize, listSize, enrichmentRatio, pVal)
        }.sortedBy { it.pValue }
    }

    data class EnrichmentResult(
        val setName: String,
        val overlap: Int,
        val setSize: Int,
        val querySize: Int,
        val enrichmentRatio: Double,
        val pValue: Double
    ) {
        override fun toString() =
            "$setName: $overlap/$setSize genes  enrichment=${"%.2f".format(enrichmentRatio)}x  p=${"%.4f".format(pValue)}"
    }

    private fun hypergeometricPValue(N: Int, K: Int, n: Int, k: Int): Double {
        // P(X >= k) under hypergeometric(N, K, n)
        var pVal = 0.0
        for (i in k..minOf(K, n)) {
            pVal += hypergeometricPMF(N, K, n, i)
        }
        return pVal.coerceIn(0.0, 1.0)
    }

    private fun hypergeometricPMF(N: Int, K: Int, n: Int, k: Int): Double {
        return exp(
            logComb(K, k) + logComb(N - K, n - k) - logComb(N, n)
        )
    }

    private fun logComb(n: Int, k: Int): Double {
        if (k < 0 || k > n) return Double.NEGATIVE_INFINITY
        if (k == 0 || k == n) return 0.0
        return lgamma((n + 1).toDouble()) - lgamma((k + 1).toDouble()) - lgamma((n - k + 1).toDouble())
    }

    private fun lgamma(x: Double): Double {
        val g = 7.0
        val c = doubleArrayOf(0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
        if (x < 0.5) return ln(PI / sin(PI * x)) - lgamma(1 - x)
        var xv = x - 1; var t = xv + g + 0.5; var s = c[0]
        for (i in 1..8) s += c[i] / (xv + i)
        return 0.5 * ln(2 * PI) + (xv + 0.5) * ln(t) - t + ln(s)
    }
}

// ─────────────────────────────────────────────────────────────
// TISSUE SPECIFICITY
// ─────────────────────────────────────────────────────────────

object TissueSpecificity {

    /** τ (tau) index — tissue specificity score (0=housekeeping, 1=tissue-specific) */
    fun tau(expressions: DoubleArray): Double {
        val max = expressions.max()!! ?: return 0.0
        if (max == 0.0) return 0.0
        val n = expressions.size
        val sum = expressions.map { 1.0 - it / max }.sum()
        return sum / (n - 1).toDouble()
    }

    /** RiceXPro tissue labels (simplified — 15 major categories) */
    val RICEXPRO_TISSUES = listOf(
        "Root", "Shoot", "Leaf_blade", "Leaf_sheath", "Stem",
        "Flag_leaf", "SAM", "Panicle", "Spikelet", "Endosperm",
        "Embryo", "Pollen", "Seed", "Callus", "Seedling"
    )

    /** Classify expression pattern based on tau value */
    fun classifyPattern(tauValue: Double): String = when {
        tauValue >= 0.85 -> "Highly tissue-specific"
        tauValue >= 0.65 -> "Tissue-enhanced"
        tauValue >= 0.35 -> "Mixed/tissue-preferred"
        else             -> "Housekeeping / broadly expressed"
    }
}
