package biokt.rice

import biokt.*
import kotlin.math.*

// ============================================================
// RICE VARIANTS — VCF / BCF parsing and variant analysis
// Subpackage: biokt.rice
//
// Covers the 3K-RGP, RiceVarMap, IRRI germplasm panel formats.
// All functionality works with standard VCF 4.1/4.2 files.
// ============================================================

// ─────────────────────────────────────────────────────────────
// VCF RECORD
// ─────────────────────────────────────────────────────────────

/** A single VCF record */
data class VcfRecord(
    val chrom: String,
    val pos: Int,                    // 1-based
    val id: String,
    val ref: String,
    val alt: List<String>,
    val qual: Double?,
    val filter: List<String>,
    val info: Map<String, String>,
    val format: List<String>,
    val genotypes: Map<String, Genotype>  // sample name → genotype
) {
    val isSNP: Boolean get() = ref.length == 1 && alt.all { it.length == 1 && it != "." }
    val isIndel: Boolean get() = !isSNP && alt.any { it != "." }
    val isMultiallelic: Boolean get() = alt.size > 1
    val isTransition: Boolean get() {
        if (!isSNP || alt.isEmpty()) return false
        val transitions = setOf("AG", "GA", "CT", "TC")
        return transitions.contains(ref + alt[0])
    }
    val isTransversion: Boolean get() = isSNP && !isTransition

    /** Minor allele frequency across all samples */
    fun maf(): Double {
        if (genotypes.isEmpty()) return 0.0
        val alleles = genotypes.values.flatMap { it.alleles() }
            .filter { it != "." && it != "*" }
        if (alleles.isEmpty()) return 0.0
        val counts = alleles.groupingBy { it }.eachCount()
        val sorted = counts.values.sortedDescending()
        if (sorted.size < 2) return 0.0
        return sorted[1].toDouble() / alleles.size
    }

    /** Observed heterozygosity */
    fun hetFrequency(): Double {
        val called = genotypes.values.filter { !it.isMissing() }
        if (called.isEmpty()) return 0.0
        return called.count { it.isHet() }.toDouble() / called.size
    }

    /** Allele frequency for each allele (ref first, then alts) */
    fun alleleFrequencies(): Map<String, Double> {
        val allAlleles = listOf(ref) + alt.filter { it != "." }
        val alleles = genotypes.values.flatMap { it.alleles() }.filter { it != "." && it != "*" }
        if (alleles.isEmpty()) return allAlleles.associateWith { 0.0 }
        return allAlleles.associateWith { a ->
            alleles.count { it == a }.toDouble() / alleles.size
        }
    }

    /** Chromosome index for sorting (handles chr1, Chr01, 1, Os01) */
    fun chromIndex(): Int {
        val num = chrom.replace(Regex("[^0-9]"), "")
        return num.toIntOrNull() ?: chrom.hashCode()
    }
}

/** Genotype for one sample at one locus */
data class Genotype(
    val gt: String,           // e.g. "0/1", "1/1", "./.", "0|1"
    val fields: Map<String, String> = emptyMap()
) {
    val isPhased: Boolean get() = gt.contains('|')
    fun isMissing(): Boolean = gt == "./." || gt == "." || gt == ".|."
    fun isHom(): Boolean = !isMissing() && alleles().toSet().size == 1
    fun isHet(): Boolean = !isMissing() && alleles().toSet().size > 1
    fun isHomRef(): Boolean = isHom() && alleles().all { it == "0" }
    fun isHomAlt(): Boolean = isHom() && alleles().all { it != "0" }

    fun alleles(): List<String> {
        if (isMissing()) return emptyList()
        val sep = if (isPhased) '|' else '/'
        return gt.split(sep)
    }

    fun dosage(): Int? {
        if (isMissing()) return null
        return alleles().count { it != "0" }
    }

    val dp: Int? get() = fields["DP"]?.toIntOrNull()
    val gq: Int? get() = fields["GQ"]?.toIntOrNull()
    val ad: List<Int> get() = fields["AD"]?.split(",")?.mapNotNull { it.toIntOrNull() } ?: emptyList()
}

// ─────────────────────────────────────────────────────────────
// VCF HEADER
// ─────────────────────────────────────────────────────────────

data class VcfHeader(
    val fileFormat: String = "VCFv4.2",
    val samples: List<String> = emptyList(),
    val contigs: Map<String, Int> = emptyMap(),       // chrom → length
    val infos: Map<String, String> = emptyMap(),       // ID → description
    val filters: Map<String, String> = emptyMap(),
    val formats: Map<String, String> = emptyMap(),
    val metadata: Map<String, String> = emptyMap()
)

// ─────────────────────────────────────────────────────────────
// VCF PARSER
// ─────────────────────────────────────────────────────────────

object VcfParser {

    /** Parse a VCF file from disk */
    fun parse(
        path: String,
        regionChrom: String? = null,
        regionStart: Int? = null,
        regionEnd: Int? = null,
        passOnly: Boolean = false,
        snpsOnly: Boolean = false,
        minMaf: Double = 0.0,
        maxMissing: Double = 1.0
    ): VcfFile {
        val lines = java.io.File(path).readLines()
        return parseLines(lines, regionChrom, regionStart, regionEnd,
            passOnly, snpsOnly, minMaf, maxMissing)
    }

    /** Parse VCF from a String (for testing / streaming) */
    fun parseString(
        text: String,
        regionChrom: String? = null,
        regionStart: Int? = null,
        regionEnd: Int? = null,
        passOnly: Boolean = false,
        snpsOnly: Boolean = false,
        minMaf: Double = 0.0,
        maxMissing: Double = 1.0
    ): VcfFile {
        val lines = text.lines()
        return parseLines(lines, regionChrom, regionStart, regionEnd,
            passOnly, snpsOnly, minMaf, maxMissing)
    }

    private fun parseLines(
        lines: List<String>,
        regionChrom: String?,
        regionStart: Int?,
        regionEnd: Int?,
        passOnly: Boolean,
        snpsOnly: Boolean,
        minMaf: Double,
        maxMissing: Double
    ): VcfFile {
        val headerLines = mutableListOf<String>()
        val records     = mutableListOf<VcfRecord>()
        var samples     = emptyList<String>()
        var headerParsed = false

        // Meta maps
        val contigs   = mutableMapOf<String, Int>()
        val infos     = mutableMapOf<String, String>()
        val filters   = mutableMapOf<String, String>()
        val formats   = mutableMapOf<String, String>()
        val metadata  = mutableMapOf<String, String>()
        var fileFormat = "VCFv4.2"

        for (line in lines) {
            when {
                line.startsWith("##") -> {
                    headerLines.add(line)
                    parseMetaLine(line, contigs, infos, filters, formats, metadata)
                    if (line.startsWith("##fileformat=")) fileFormat = line.substringAfter("=")
                }
                line.startsWith("#CHROM") -> {
                    headerLines.add(line)
                    val cols = line.trimStart('#').split("\t")
                    samples = if (cols.size > 9) cols.drop(9) else emptyList()
                    headerParsed = true
                }
                else -> { /* non-header line, process below */ }
            }
            // Process data lines (outside when block - continue not allowed inside when)
            if (headerParsed && !line.startsWith("#") && line.isNotBlank()) {
                val rec = parseRecord(line, samples)
                if (rec != null) {
                    val pass = (regionChrom == null || rec.chrom == regionChrom)
                        && (regionStart == null || rec.pos >= regionStart)
                        && (regionEnd   == null || rec.pos <= regionEnd)
                        && (!passOnly || rec.filter.isEmpty() || rec.filter.contains("PASS"))
                        && (!snpsOnly || rec.isSNP)
                    if (pass) {
                        val miss = rec.genotypes.values.count { it.isMissing() }.toDouble() /
                                   rec.genotypes.size.toDouble().coerceAtLeast(1.0)
                        if (miss <= maxMissing && (minMaf <= 0 || rec.maf() >= minMaf))
                            records.add(rec)
                    }
                }
            }
        }

        val header = VcfHeader(fileFormat, samples, contigs, infos, filters, formats, metadata)
        return VcfFile(header, records)
    }

    private fun parseMetaLine(
        line: String,
        contigs: MutableMap<String, Int>,
        infos: MutableMap<String, String>,
        filters: MutableMap<String, String>,
        formats: MutableMap<String, String>,
        metadata: MutableMap<String, String>
    ) {
        when {
            line.startsWith("##contig=") -> {
                val id  = Regex("ID=([^,>]+)").find(line)?.groupValues?.get(1) ?: return
                val len = Regex("length=(\\d+)").find(line)?.groupValues?.get(1)?.toIntOrNull()
                contigs[id] = len ?: 0
            }
            line.startsWith("##INFO=") -> {
                val id   = Regex("ID=([^,]+)").find(line)?.groupValues?.get(1) ?: return
                val desc = Regex("Description=\"([^\"]+)\"").find(line)?.groupValues?.get(1) ?: ""
                infos[id] = desc
            }
            line.startsWith("##FILTER=") -> {
                val id   = Regex("ID=([^,]+)").find(line)?.groupValues?.get(1) ?: return
                val desc = Regex("Description=\"([^\"]+)\"").find(line)?.groupValues?.get(1) ?: ""
                filters[id] = desc
            }
            line.startsWith("##FORMAT=") -> {
                val id   = Regex("ID=([^,]+)").find(line)?.groupValues?.get(1) ?: return
                val desc = Regex("Description=\"([^\"]+)\"").find(line)?.groupValues?.get(1) ?: ""
                formats[id] = desc
            }
            else -> {
                val kv = line.removePrefix("##").split("=", limit = 2)
                if (kv.size == 2) metadata[kv[0]] = kv[1]
            }
        }
    }

    private fun parseRecord(line: String, samples: List<String>): VcfRecord? {
        val cols = line.split("\t")
        if (cols.size < 8) return null
        return try {
            val chrom  = cols[0]
            val pos    = cols[1].toInt()
            val id     = cols[2]
            val ref    = cols[3].toUpperCase()
            val alt    = cols[4].split(",").map { it.toUpperCase() }
            val qual   = cols[5].takeIf { it != "." }?.toDoubleOrNull()
            val filter = if (cols[6] == "." || cols[6] == "PASS") emptyList() else cols[6].split(";")
            val info   = parseInfo(cols[7])
            val format = if (cols.size > 8) cols[8].split(":") else emptyList()
            val gts    = if (cols.size > 9 && samples.isNotEmpty()) {
                samples.zip(cols.drop(9)).associate { (s, g) ->
                    s to parseGenotype(g, format)
                }
            } else emptyMap()
            VcfRecord(chrom, pos, id, ref, alt, qual, filter, info, format, gts)
        } catch (e: Exception) { null }
    }

    private fun parseInfo(info: String): Map<String, String> {
        if (info == ".") return emptyMap()
        return info.split(";").associate { entry ->
            val kv = entry.split("=", limit = 2)
            if (kv.size == 2) kv[0] to kv[1] else kv[0] to "true"
        }
    }

    private fun parseGenotype(raw: String, format: List<String>): Genotype {
        val vals = raw.split(":")
        val gt   = vals.firstOrNull() ?: "./."
        val flds = format.drop(1).zip(vals.drop(1)).associate { (k, v) -> k to v }
        return Genotype(gt, flds)
    }
}

// ─────────────────────────────────────────────────────────────
// VCF FILE — container + analysis methods
// ─────────────────────────────────────────────────────────────

class VcfFile(
    val header: VcfHeader,
    val records: List<VcfRecord>
) {
    val size: Int get() = records.size
    val samples: List<String> get() = header.samples

    // ── Filtering ────────────────────────────────────────────

    fun snps() = VcfFile(header, records.filter { it.isSNP })
    fun indels() = VcfFile(header, records.filter { it.isIndel })
    fun biallelic() = VcfFile(header, records.filter { !it.isMultiallelic && it.alt.size == 1 })
    fun onChrom(chrom: String) = VcfFile(header, records.filter { it.chrom == chrom })
    fun inRegion(chrom: String, start: Int, end: Int) =
        VcfFile(header, records.filter { it.chrom == chrom && it.pos in start..end })
    fun minMaf(threshold: Double) = VcfFile(header, records.filter { it.maf() >= threshold })
    fun maxMissing(fraction: Double) = VcfFile(header, records.filter { r ->
        r.genotypes.values.count { it.isMissing() }.toDouble() /
        r.genotypes.size.toDouble().coerceAtLeast(1.0) <= fraction
    })

    // ── Summary statistics ────────────────────────────────────

    fun summarise(): VcfSummary {
        val snps       = records.count { it.isSNP }
        val indels     = records.count { it.isIndel }
        val multiall   = records.count { it.isMultiallelic }
        val ts         = records.count { it.isTransition }
        val tv         = records.count { it.isTransversion }
        val tsTV       = if (tv > 0) ts.toDouble() / tv else 0.0
        val mafs       = records.map { it.maf() }
        val avgMaf     = if (mafs.isEmpty()) 0.0 else mafs.average()
        val singletons = mafs.count { it > 0 && it <= 1.0 / (2 * samples.size).coerceAtLeast(1) }
        val missingRates = records.map { r ->
            r.genotypes.values.count { it.isMissing() }.toDouble() /
            r.genotypes.size.toDouble().coerceAtLeast(1.0)
        }
        val avgMissing = if (missingRates.isEmpty()) 0.0 else missingRates.average()
        return VcfSummary(records.size, snps, indels, multiall, ts, tv, tsTV,
            avgMaf, singletons, avgMissing, samples.size)
    }

    // ── Genotype matrix ───────────────────────────────────────
    // Returns dosage matrix: rows=variants, cols=samples, values=0/1/2 (-1=missing)

    fun dosageMatrix(): Array<IntArray> {
        return Array(records.size) { i ->
            val rec = records[i]
            IntArray(samples.size) { j ->
                rec.genotypes[samples[j]]?.dosage() ?: -1
            }
        }
    }

    // ── LD (Linkage Disequilibrium) ───────────────────────────

    fun ldMatrix(maxPairs: Int = 500): LDMatrix {
        val n = minOf(records.size, maxPairs)
        val subset = records.take(n)
        val r2 = Array(n) { DoubleArray(n) { 0.0 } }
        val dprime = Array(n) { DoubleArray(n) { 0.0 } }
        for (i in 0 until n) {
            r2[i][i] = 1.0; dprime[i][i] = 1.0
            for (j in i + 1 until n) {
                val ld = computeLD(subset[i], subset[j])
                r2[i][j] = ld.first; r2[j][i] = ld.first
                dprime[i][j] = ld.second; dprime[j][i] = ld.second
            }
        }
        return LDMatrix(subset.map { "${it.chrom}:${it.pos}" }, r2, dprime)
    }

    fun ldBetween(variantA: VcfRecord, variantB: VcfRecord): Pair<Double, Double> =
        computeLD(variantA, variantB)

    private fun computeLD(a: VcfRecord, b: VcfRecord): Pair<Double, Double> {
        val sharedSamples = samples.filter {
            !a.genotypes[it]?.isMissing()!! && !b.genotypes[it]?.isMissing()!!
        }
        if (sharedSamples.size < 4) return 0.0 to 0.0

        // Haplotype counts (assuming Hardy-Weinberg, using expected haplotype frequencies)
        var n11 = 0.0; var n10 = 0.0; var n01 = 0.0; var n00 = 0.0

        sharedSamples.forEach { s ->
            val gtA = a.genotypes[s]!!.dosage() ?: return@forEach
            val gtB = b.genotypes[s]!!.dosage() ?: return@forEach
            // Dosage: 0=ref/ref, 1=ref/alt, 2=alt/alt
            val pA = gtA / 2.0  // estimated alt allele freq in this individual
            val pB = gtB / 2.0
            n11 += pA * pB
            n10 += pA * (1 - pB)
            n01 += (1 - pA) * pB
            n00 += (1 - pA) * (1 - pB)
        }

        val n = sharedSamples.size.toDouble()
        val p1 = (n11 + n10) / n  // alt freq at A
        val p2 = (n11 + n01) / n  // alt freq at B
        val q1 = 1 - p1; val q2 = 1 - p2

        val D = n11 / n - p1 * p2
        val denom = p1 * q1 * p2 * q2
        if (denom <= 0) return 0.0 to 0.0

        val r2 = (D * D) / denom
        val dMax = if (D > 0) minOf(p1 * q2, q1 * p2) else minOf(p1 * p2, q1 * q2)
        val dprime = if (dMax == 0.0) 0.0 else (D / dMax).absoluteValue

        return r2.coerceIn(0.0, 1.0) to dprime.coerceIn(0.0, 1.0)
    }

    // ── Population statistics ─────────────────────────────────

    fun snpDensity(windowSize: Int = 100_000): List<WindowStat> {
        val grouped = records.filter { it.isSNP }.groupBy { it.chrom }
        val result  = mutableListOf<WindowStat>()
        grouped.forEach { (chrom, snps) ->
            val maxPos  = snps.maxBy { it.pos }?.pos ?: return@forEach
            var winStart = 1
            while (winStart <= maxPos) {
                val winEnd = winStart + windowSize - 1
                val count  = snps.count { it.pos in winStart..winEnd }
                result.add(WindowStat(chrom, winStart, winEnd, count.toDouble()))
                winStart = winEnd + 1
            }
        }
        return result.sortedWith(compareBy({ it.chrom }, { it.start }))
    }

    fun nucleotideDiversity(windowSize: Int = 100_000): List<WindowStat> {
        if (samples.isEmpty()) return emptyList()
        val n    = samples.size
        val snps = records.filter { it.isSNP }
        val grouped = snps.groupBy { it.chrom }
        val result  = mutableListOf<WindowStat>()
        grouped.forEach { (chrom, chromSnps) ->
            val maxPos = chromSnps.maxBy { it.pos }?.pos ?: return@forEach
            var ws = 1
            while (ws <= maxPos) {
                val we   = ws + windowSize - 1
                val win  = chromSnps.filter { it.pos in ws..we }
                var piSum = 0.0
                win.forEach { snp ->
                    val p = snp.alleleFrequencies()[snp.ref] ?: 0.0
                    val q = 1 - p
                    piSum += 2 * p * q
                }
                val pi = if (win.isNotEmpty()) piSum / windowSize else 0.0
                result.add(WindowStat(chrom, ws, we, pi))
                ws = we + 1
            }
        }
        return result
    }

    // ── Sample subsetting ─────────────────────────────────────

    fun subset(sampleNames: List<String>): VcfFile {
        val filteredRecords = records.map { rec ->
            rec.copy(genotypes = rec.genotypes.filterKeys { it in sampleNames })
        }
        return VcfFile(header.copy(samples = sampleNames), filteredRecords)
    }

    // ── VCF output ────────────────────────────────────────────

    fun writeVcf(path: String) {
        java.io.File(path).bufferedWriter().use { w ->
            w.write("##fileformat=${header.fileFormat}\n")
            header.metadata.forEach { (k, v) -> w.write("##$k=$v\n") }
            header.contigs.forEach { (id, len) ->
                w.write("##contig=<ID=$id${if (len > 0) ",length=$len" else ""}>\n")
            }
            val colHeader = listOf("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") + samples
            w.write("#${colHeader.joinToString("\t")}\n")
            records.forEach { rec ->
                val filterStr = if (rec.filter.isEmpty()) "PASS" else rec.filter.joinToString(";")
                val infoStr   = if (rec.info.isEmpty()) "." else rec.info.entries.joinToString(";") { "${it.key}=${it.value}" }
                val fmtStr    = if (rec.format.isEmpty()) "GT" else rec.format.joinToString(":")
                val gtFields  = samples.map { s ->
                    rec.genotypes[s]?.let { gt ->
                        if (rec.format.size <= 1) gt.gt
                        else {
                            val vals = mutableListOf(gt.gt)
                            rec.format.drop(1).forEach { f -> vals.add(gt.fields[f] ?: ".") }
                            vals.joinToString(":")
                        }
                    } ?: "./."
                }
                val row = listOf(rec.chrom, rec.pos.toString(), rec.id, rec.ref,
                    rec.alt.joinToString(","), rec.qual?.toString() ?: ".",
                    filterStr, infoStr, fmtStr) + gtFields
                w.write(row.joinToString("\t") + "\n")
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
// RESULT DATA CLASSES
// ─────────────────────────────────────────────────────────────

data class VcfSummary(
    val totalVariants: Int,
    val snps: Int,
    val indels: Int,
    val multiallelic: Int,
    val transitions: Int,
    val transversions: Int,
    val tsTvRatio: Double,
    val avgMaf: Double,
    val singletons: Int,
    val avgMissingRate: Double,
    val numSamples: Int
) {
    fun print() = buildString {
        append("=== VCF Summary ===\n")
        append("Samples        : $numSamples\n")
        append("Total variants : $totalVariants\n")
        append("SNPs           : $snps\n")
        append("Indels         : $indels\n")
        append("Multiallelic   : $multiallelic\n")
        append("Ts/Tv ratio    : ${"%.3f".format(tsTvRatio)}\n")
        append("Avg MAF        : ${"%.4f".format(avgMaf)}\n")
        append("Singletons     : $singletons\n")
        append("Avg missing    : ${"%.3f".format(avgMissingRate)}\n")
    }
}

data class LDMatrix(
    val loci: List<String>,
    val r2: Array<DoubleArray>,
    val dprime: Array<DoubleArray>
) {
    fun topPairs(n: Int = 10, metric: String = "r2"): List<Triple<String, String, Double>> {
        val mat = if (metric == "dprime") dprime else r2
        val pairs = mutableListOf<Triple<String, String, Double>>()
        for (i in loci.indices) for (j in i + 1 until loci.size)
            pairs.add(Triple(loci[i], loci[j], mat[i][j]))
        return pairs.sortedByDescending { it.third }.take(n)
    }
}

data class WindowStat(
    val chrom: String,
    val start: Int,
    val end: Int,
    val value: Double
)

// ─────────────────────────────────────────────────────────────
// RICE-SPECIFIC VARIANT UTILITIES
// ─────────────────────────────────────────────────────────────

object RiceVariantUtils {

    // Known rice chromosome names and lengths (Nipponbare IRGSP-1.0)
    val RICE_CHROMOSOMES = mapOf(
        "chr01" to 43_270_923, "chr02" to 35_937_250, "chr03" to 36_413_819,
        "chr04" to 35_502_694, "chr05" to 29_958_434, "chr06" to 31_248_787,
        "chr07" to 29_697_621, "chr08" to 28_443_022, "chr09" to 23_012_720,
        "chr10" to 23_207_287, "chr11" to 29_021_106, "chr12" to 27_531_856,
        // Also accept Os prefix and number-only
        "Os01" to 43_270_923, "Os02" to 35_937_250, "Os03" to 36_413_819,
        "1" to 43_270_923, "2" to 35_937_250, "3" to 36_413_819
    )

    /** Normalise chromosome names to chr01..chr12 */
    fun normaliseChrom(raw: String): String {
        val n = raw.replace(Regex("[^0-9]"), "").trimStart('0')
        val num = n.toIntOrNull() ?: return raw
        return "chr${num.toString().padStart(2, '0')}"
    }

    /** Known rice domestication SNPs (selected examples from published studies) */
    data class DomLocus(val chrom: String, val pos: Int, val gene: String, val effect: String)
    val DOMESTICATION_SNPS = listOf(
        // Format: chrom, pos, gene, effect
        DomLocus("chr01", 38_286_560, "Rc", "Red pericarp loss (f deletion, japonica)"),
        DomLocus("chr04", 3_545_070, "Wx", "Waxy starch — G→T splice site, glutinous"),
        DomLocus("chr07", 13_476_000, "sh4", "Shattering — K→N substitution, domesticated"),
        DomLocus("chr03", 9_027_843, "Prog1", "Prostrate growth — domestication sweep"),
        DomLocus("chr01", 6_691_563, "GS3", "Grain size control"),
        DomLocus("chr06", 1_621_227, "GW5", "Grain width — major QTL"),
        DomLocus("chr08", 25_563_000, "OsC1", "Anthocyanin accumulation"),
        DomLocus("chr11", 4_267_500, "Hd1", "Heading date — photoperiod response")
    )

    /** Check if a VCF record overlaps a known rice domestication locus */
    fun isDomesticationSNP(rec: VcfRecord, windowBp: Int = 5000): String? {
        val normChrom = normaliseChrom(rec.chrom)
        return DOMESTICATION_SNPS.firstOrNull { d ->
            d.chrom == normChrom && kotlin.math.abs(d.pos - rec.pos) <= windowBp
        }?.let { "${it.gene}: ${it.effect}" }
    }

    /** Classify variants from a 3K-RGP style panel into subpopulations
     *  Uses simple MAF and heterozygosity thresholds as proxies */
    fun inferSubpopGroup(
        vcf: VcfFile,
        sampleName: String,
        minVariants: Int = 100
    ): String {
        val gts = vcf.records.take(minVariants).mapNotNull { rec ->
            rec.genotypes[sampleName]
        }
        if (gts.size < minVariants / 2) return "Unknown"
        val hetRate  = gts.count { it.isHet() }.toDouble() / gts.size
        val homAlt   = gts.count { it.isHomAlt() }.toDouble() / gts.size
        return if (hetRate > 0.15) "Admixed"
               else if (homAlt > 0.45) "Indica"
               else if (homAlt < 0.20) "Japonica"
               else "Intermediate"
    }
}
