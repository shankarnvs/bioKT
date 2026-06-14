package biokt.rice

import biokt.*
import kotlin.math.*

// ============================================================
// RICE GENOMICS - rice-specific biological knowledge base
// Subpackage: biokt.rice
//
// SSR detection, cis-element scanning, transposable element
// classification, R-gene domain patterns, heading date gene
// catalogue, chromosome coordinates, and rice ontology.
// ============================================================

// ─────────────────────────────────────────────────────────────
// SSR - Simple Sequence Repeat (microsatellite) detection
// ─────────────────────────────────────────────────────────────

data class SSRMarker(
    val id: String,
    val start: Int,            // 1-based
    val end: Int,
    val motif: String,         // e.g. "AT", "AGC"
    val motifLength: Int,      // 1–6
    val repeats: Int,
    val sequence: String,      // full SSR sequence
    val type: String           // "mono", "di", "tri", "tetra", "penta", "hexa"
) {
    val length: Int get() = end - start + 1
    override fun toString() = "SSR[$id] ${motif}×$repeats at ${start}-$end"
}

object SSRDetector {

    private val TYPE_NAMES = mapOf(1 to "mono", 2 to "di", 3 to "tri",
        4 to "tetra", 5 to "penta", 6 to "hexa")

    // Minimum repeat counts per motif length (MISA defaults)
    private val MIN_REPEATS = mapOf(1 to 10, 2 to 6, 3 to 5, 4 to 5, 5 to 4, 6 to 4)

    /**
     * Detect all SSRs in a DNA sequence.
     * motifLengths: which repeat unit lengths to scan (default 1–6)
     * minRepeats: minimum number of repeats (uses MISA defaults if null)
     */
    fun detect(
        seq: DNASequence,
        motifLengths: IntRange = 1..6,
        minRepeats: Map<Int, Int>? = null
    ): List<SSRMarker> {
        val s = seq.sequence
        val n = s.length
        val ssrs = mutableListOf<SSRMarker>()
        val minR = minRepeats ?: MIN_REPEATS
        var ssrCounter = 0

        for (motifLen in motifLengths) {
            val minRep = minR[motifLen] ?: 4
            var i = 0
            while (i <= n - motifLen * minRep) {
                val motif = s.substring(i, i + motifLen)
                if (motif.length != motifLen) { i++; continue }

                // Count how many times this motif repeats
                var repeats = 1
                var j = i + motifLen
                while (j + motifLen <= n && s.substring(j, j + motifLen) == motif) {
                    repeats++
                    j += motifLen
                }

                if (repeats >= minRep) {
                    ssrCounter++
                    val end = i + repeats * motifLen
                    ssrs.add(SSRMarker(
                        id         = "${seq.id.ifEmpty{"seq"}}_ssr$ssrCounter",
                        start      = i + 1,
                        end        = end,
                        motif      = motif,
                        motifLength = motifLen,
                        repeats    = repeats,
                        sequence   = s.substring(i, end),
                        type       = TYPE_NAMES[motifLen] ?: "hexa"
                    ))
                    i = end  // skip past this SSR
                } else {
                    i++
                }
            }
        }
        return ssrs.sortedBy { it.start }
    }

    /** Summary statistics for detected SSRs */
    fun summarise(ssrs: List<SSRMarker>): Map<String, Any> {
        val byType = ssrs.groupBy { it.type }
        return mapOf(
            "total"         to ssrs.size,
            "by_type"       to byType.mapValues { it.value.size },
            "most_common_motif" to (ssrs.groupBy { it.motif }.maxBy { it.value.size }?.key ?: "none"),
            "avg_repeats"   to if (ssrs.isEmpty()) 0.0 else ssrs.map { it.repeats }.average(),
            "total_length"  to ssrs.map { it.length }.sum()
        )
    }
}

// ─────────────────────────────────────────────────────────────
// CIS-ELEMENT SCANNER
// Plant promoter cis-regulatory element detection
// Based on PLACE, JASPAR, and rice-specific literature
// ─────────────────────────────────────────────────────────────

data class CisElement(
    val name: String,
    val motif: String,
    val position: Int,
    val strand: Char,
    val function: String,
    val family: String
)

object CisElementScanner {

    /** Rice promoter cis-regulatory elements database */
    val ELEMENTS = listOf(
        // ABA-responsive elements
        CisElementDef("ABRE",       "ACGTGG|ACGTGT|CACGTG|TACGTC", "ABA Response Element - drought/stress", "bZIP"),
        CisElementDef("ABRE-like",  "ACGTG",                        "ABA Response Element variant", "bZIP"),
        CisElementDef("DRE/CRT",    "TACCGACAT|RCCGAC",            "Drought Response Element / C-repeat", "DREB/CBF"),
        CisElementDef("MYCBS",      "CATGTG",                       "MYC binding site - ABA signalling", "MYC"),
        CisElementDef("MYBBS",      "TAACTG",                       "MYB binding site - ABA/drought", "MYB"),
        // Light-responsive
        CisElementDef("G-box",      "CACGTG",                       "G-box - light response, circadian", "bHLH/bZIP"),
        CisElementDef("I-box",      "GATAAG|GATAAN",               "I-box - light-responsive promoter", "GATA"),
        CisElementDef("GT-1",       "GRWAAW",                       "GT-1 element - light response", "GT factor"),
        CisElementDef("GATA-motif", "TATCTA",                       "GATA motif - light regulation", "GATA"),
        // Defence/pathogen
        CisElementDef("W-box",      "TTGACY|TTGACT|TTGACC",        "W-box - WRKY binding, defence/JA signalling", "WRKY"),
        CisElementDef("GCC-box",    "AGCCGCC",                      "GCC box - ethylene response / defence", "AP2/ERF"),
        CisElementDef("as-1/TCA",   "TGACG",                        "as-1 element - SA response, TCA motif", "TGA/bZIP"),
        CisElementDef("CGTCA",      "CGTCA",                        "CGTCA motif - jasmonate response", "MYC/MYB"),
        CisElementDef("TGACG",      "TGACG",                        "TGACG motif - SA/JA signalling", "TGA"),
        // Hormones
        CisElementDef("AuxRE",      "TGASTC|TGACGC",               "Auxin Response Element", "ARF"),
        CisElementDef("GARE",       "TAACAAA|TAACAAR",             "Gibberellin Response Element", "GAMYB"),
        CisElementDef("CARE",       "TGACG",                        "Cytokinin Activation Response Element", "ARR"),
        CisElementDef("P-box",      "CCTTTTG",                      "P-box / GARE2 - gibberellin response", "GAMYB"),
        // Development
        CisElementDef("RY-element", "CATGCAT|CATGCAC",             "RY/Sph element - seed-specific expression", "VP1/B3"),
        CisElementDef("GCN4",       "TGACTCAT|TGAGTCAT",          "GCN4 motif - endosperm expression", "bZIP"),
        CisElementDef("AACA",       "AACAAAC",                      "AACA motif - nitrogen starvation / leaf senescence", "MYB"),
        CisElementDef("CCAAT-box",  "CCAAT",                        "CCAAT box - constitutive / tissue-specific", "NF-Y"),
        CisElementDef("TATA-box",   "TATAAATA|TATAAAAT|TATAAAT",  "TATA box - core promoter element", "TBP"),
        CisElementDef("Inr",        "TCANTCT",                      "Initiator element - transcription start", "TBP"),
        // Circadian
        CisElementDef("EE",         "AAATATCT",                     "Evening Element - circadian clock output", "RVE/CCA1"),
        CisElementDef("Morning-E",  "AACCACGAAAC",                  "Morning Element - circadian clock", "PRR"),
        // Stress-specific (rice)
        CisElementDef("TC-rich",    "ATTTTCTTCA",                   "TC-rich repeat - defence/stress", "RAV"),
        CisElementDef("HSE",        "AGAANTTCT|AGAANNTTCNAGAA",   "Heat Shock Element", "HSF"),
        CisElementDef("ARE",        "AAACCA",                       "Anaerobic Response Element - submergence", "RAV/ERF"),
        CisElementDef("SURE",       "TGATAATG",                     "SUlphur Responsive Element", "bZIP"),
        CisElementDef("PIRF",       "CACGCAACGC",                   "Phosphate starvation response element", "PHR1/MYB")
    )

    data class CisElementDef(val name: String, val patterns: String, val function: String, val family: String) {
        val motifList: List<String> get() = patterns.split("|")
    }

    private val IUPAC_MAP = mapOf(
        'R' to "[AG]", 'Y' to "[CT]", 'S' to "[GC]", 'W' to "[AT]",
        'K' to "[GT]", 'M' to "[AC]", 'B' to "[CGT]", 'D' to "[AGT]",
        'H' to "[ACT]", 'V' to "[ACG]", 'N' to "[ACGT]"
    )

    private fun motifToRegex(motif: String): Regex {
        val pattern = motif.toUpperCase().map { c ->
            IUPAC_MAP[c] ?: c.toString()
        }.joinToString("")
        return Regex(pattern)
    }

    /**
     * Scan a sequence for cis-regulatory elements.
     * Returns all elements found on both strands.
     */
    fun scan(
        seq: DNASequence,
        elements: List<CisElementDef> = ELEMENTS,
        bothStrands: Boolean = true
    ): List<CisElement> {
        val results = mutableListOf<CisElement>()
        val s = seq.sequence

        elements.forEach { el ->
            el.motifList.forEach { motif ->
                val re = motifToRegex(motif)
                // Forward strand
                re.findAll(s).forEach { m ->
                    results.add(CisElement(el.name, motif, m.range.first + 1, '+', el.function, el.family))
                }
                // Reverse complement strand
                if (bothStrands) {
                    val rc = DNASequence(s).reverseComplement().sequence
                    re.findAll(rc).forEach { m ->
                        val pos = s.length - m.range.last  // convert to forward coords
                        results.add(CisElement(el.name, motif, pos, '-', el.function, el.family))
                    }
                }
            }
        }
        return results.sortedBy { it.position }
    }

    /** Summarise scan results by element family */
    fun summarise(elements: List<CisElement>): Map<String, Int> =
        elements.groupBy { it.name }.mapValues { it.value.size }.entries
            .sortedByDescending { it.value }
            .associate { it.key to it.value }
}

// ─────────────────────────────────────────────────────────────
// TRANSPOSABLE ELEMENT CLASSIFIER
// Rice genome is ~45% TEs - classification is critical for
// repeat masking, annotation quality, and genome evolution
// ─────────────────────────────────────────────────────────────

enum class TEClass { LTR_RETROTRANSPOSON, NON_LTR_RETROTRANSPOSON, DNA_TRANSPOSON, UNKNOWN }
enum class TEFamily { GYPSY, COPIA, LINE, SINE, MUTATOR, CACTA, HARBINGER, PIF_HARBINGER,
    TOURIST_MINIATURE, STOWAWAY, UNKNOWN }

data class TEAnnotation(
    val id: String,
    val start: Int,
    val end: Int,
    val strand: Char,
    val teClass: TEClass,
    val family: TEFamily,
    val score: Double,
    val evidence: List<String>
) {
    val length: Int get() = end - start + 1
    override fun toString() = "TE[$id] ${teClass.name}/${family.name} ${start}-${end}"
}

object TEClassifier {

    // Protein domain signatures for TE classification
    private val DOMAIN_SIGNATURES = mapOf(
        // LTR retrotransposons
        "RVT" to Pair(TEClass.LTR_RETROTRANSPOSON, TEFamily.UNKNOWN),
        "GYPSY" to Pair(TEClass.LTR_RETROTRANSPOSON, TEFamily.GYPSY),
        "COPIA" to Pair(TEClass.LTR_RETROTRANSPOSON, TEFamily.COPIA),
        "GAG" to Pair(TEClass.LTR_RETROTRANSPOSON, TEFamily.UNKNOWN),
        "INTEGRASE" to Pair(TEClass.LTR_RETROTRANSPOSON, TEFamily.UNKNOWN),
        // Non-LTR
        "LINE" to Pair(TEClass.NON_LTR_RETROTRANSPOSON, TEFamily.LINE),
        "SINE" to Pair(TEClass.NON_LTR_RETROTRANSPOSON, TEFamily.SINE),
        "ENDONUCLEASE" to Pair(TEClass.NON_LTR_RETROTRANSPOSON, TEFamily.LINE),
        // DNA transposons
        "TRANSPOSASE" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.UNKNOWN),
        "MUTATOR" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.MUTATOR),
        "CACTA" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.CACTA),
        "HARBINGER" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.HARBINGER),
        "TOURIST" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.TOURIST_MINIATURE),
        "STOWAWAY" to Pair(TEClass.DNA_TRANSPOSON, TEFamily.STOWAWAY)
    )

    // TIR (Terminal Inverted Repeat) motifs characteristic of DNA TEs
    private val TIR_MOTIFS = mapOf(
        TEFamily.MUTATOR    to listOf("CAGTGGATG"),
        TEFamily.CACTA      to listOf("CACT", "GTGA"),
        TEFamily.HARBINGER  to listOf("TTAC", "GTAA"),
        TEFamily.TOURIST_MINIATURE to listOf("TGGGGG", "CCCCA")
    )

        /**
     * Classify a sequence as a TE based on:
     * 1. Known domain keywords in sequence ID/name
     * 2. TIR motif detection at sequence ends
     * 3. Sequence composition (high proportion of tandem repeats)
     */
    fun classify(seq: DNASequence, name: String = seq.id): TEAnnotation? {
        val evidence = mutableListOf<String>()

        // Check name-based signatures
        val upperName = name.toUpperCase()
        var bestClass = TEClass.UNKNOWN
        var bestFamily = TEFamily.UNKNOWN
        var bestScore = 0.0

        DOMAIN_SIGNATURES.forEach { (keyword, classFam) ->
            if (upperName.contains(keyword)) {
                evidence.add("Name contains '$keyword'")
                bestClass = classFam.first
                bestFamily = classFam.second
                bestScore = 0.8
            }
        }

        // Check TIR motifs
        val s = seq.sequence
        TIR_MOTIFS.forEach { (family, motifs) ->
            val hasTIR = motifs.any { m ->
                s.startsWith(m) || s.endsWith(m) ||
                s.startsWith(DNASequence(m).reverseComplement().sequence) ||
                s.endsWith(DNASequence(m).reverseComplement().sequence)
            }
            if (hasTIR) {
                evidence.add("TIR motif for $family")
                if (bestScore < 0.7) {
                    bestClass = TEClass.DNA_TRANSPOSON
                    bestFamily = family
                    bestScore = 0.7
                }
            }
        }

        // Compositional check - TEs often have low complexity
        val complexity = SeqStats.linguisticComplexity(seq)
        if (complexity < 0.35) {
            evidence.add("Low sequence complexity (${"%.2f".format(complexity)}) - repeat-rich")
            if (bestScore < 0.5) bestScore = 0.4
        }

        if (bestScore == 0.0) return null

        return TEAnnotation(name, 1, seq.length, '+', bestClass, bestFamily, bestScore, evidence)
    }

    /**
     * Estimate TE content in a genome sequence (large window scan).
     * Returns fraction of bases that are in TE-like regions.
     */
    fun estimateTEContent(seq: DNASequence, windowSize: Int = 500): Double {
        val s = seq.sequence
        var teWindows = 0
        val totalWindows = s.length / windowSize
        if (totalWindows == 0) return 0.0

        for (i in 0 until totalWindows) {
            val win = s.substring(i * windowSize, minOf((i + 1) * windowSize, s.length))
            val winSeq = DNASequence(win)
            val complexity = SeqStats.linguisticComplexity(winSeq)
            val gcContent  = winSeq.gcContent()
            // Heuristic: low complexity OR extreme GC → likely repetitive/TE
            if (complexity < 0.40 || gcContent > 75 || gcContent < 25) teWindows++
        }
        return teWindows.toDouble() / totalWindows
    }
}

// ─────────────────────────────────────────────────────────────
// R-GENE (DISEASE RESISTANCE GENE) PATTERNS
// ─────────────────────────────────────────────────────────────

data class RGenePattern(
    val type: String,
    val domains: List<String>,
    val description: String,
    val examples: List<String>
)

object RGeneClassifier {

    val R_GENE_TYPES = listOf(
        RGenePattern("NBS-LRR (CC-NBS-LRR)", listOf("CC","NBS","LRR"),
            "Coiled-coil + NBS + LRR - largest class in rice (~480 genes)",
            listOf("Pita", "Pi9", "Pi2", "Pik-h")),
        RGenePattern("NBS-LRR (TIR-NBS-LRR)", listOf("TIR","NBS","LRR"),
            "TIR domain + NBS + LRR - rare in monocots",
            listOf()),
        RGenePattern("RLK (LRR-RLK)", listOf("SP","LRR","TM","KD"),
            "Leucine-rich repeat receptor-like kinase",
            listOf("Pi-d2","Xa21","Xa26","Xa3/Xa26")),
        RGenePattern("RLP (LRR-RLP)", listOf("SP","LRR","TM"),
            "LRR receptor-like protein (no kinase domain)",
            listOf()),
        RGenePattern("Lectin-RLK", listOf("LECTIN","TM","KD"),
            "Lectin receptor-like kinase",
            listOf("Xa26","OsWAK")),
        RGenePattern("Executor (OsHOX)", listOf("Executor"),
            "Transcription factor - executor R gene",
            listOf("Xa10","Xa23","Xa27")),
        RGenePattern("Unknown", listOf(),
            "Non-canonical resistance gene",
            listOf("Pi21","Pid3"))
    )

    /** Classify a protein sequence as a potential R-gene */
    fun classify(prot: ProteinSequence): List<String> {
        val evidence = mutableListOf<String>()
        val s = prot.sequence

        // NBS domain hallmarks (P-loop: GXXXXGKS/T, kinase-2b: RNBS-D: FLHIACFP)
        val pLoopRe   = Regex("G.{4}GK[ST]")
        val kinase2b  = Regex("RNBS")
        val hhamplet  = Regex("MHDV|MHDD|CFLH")
        if (pLoopRe.containsMatchIn(s))  evidence.add("P-loop (NBS domain)")
        if (kinase2b.containsMatchIn(s)) evidence.add("RNBS motif (NBS domain)")
        if (hhamplet.containsMatchIn(s)) evidence.add("HHAMPLET/MHD motif (NBS domain)")

        // LRR hallmarks (LxxLxL pattern repeated)
        val lrr = Regex("L.{1,2}L.{1,2}L.{1,2}L")
        if (lrr.containsMatchIn(s) && lrr.findAll(s).count() >= 2)
            evidence.add("LRR repeats (≥2 found)")

        // CC domain proxy - leucine zipper
        val lz = Regex("L.{6}L.{6}L.{6}L")
        if (lz.containsMatchIn(s)) evidence.add("Putative coiled-coil / CC domain")

        // TIR domain - TIR signature in N-terminal region
        val tir = Regex("[FW]LR[ND][LI]")
        if (tir.containsMatchIn(s.take(300))) evidence.add("TIR domain signature (N-terminal)")

        // Kinase domain
        val kd = Regex("DFG|HRDL|VAIK")
        if (kd.containsMatchIn(s)) evidence.add("Kinase domain motif (RLK candidate)")

        // Signal peptide proxy - hydrophobic N-terminal stretch
        if (s.length > 25) {
            val nTerm = s.take(25)
            val hydro = nTerm.count { it in "AVILMFWP" }.toDouble() / 25
            if (hydro > 0.5) evidence.add("Hydrophobic N-terminus (signal peptide / TM?)")
        }

        return evidence
    }

    /** Is this protein likely an R-gene based on domain analysis? */
    fun isLikelyRGene(prot: ProteinSequence): Pair<Boolean, List<String>> {
        val evidence = classify(prot)
        val hasNBS   = evidence.any { it.contains("NBS") || it.contains("P-loop") }
        val hasLRR   = evidence.any { it.contains("LRR") }
        val hasKD    = evidence.any { it.contains("Kinase") }
        return (hasNBS || (hasLRR && hasKD)) to evidence
    }
}

// ─────────────────────────────────────────────────────────────
// HEADING DATE & FLOWERING GENE CATALOGUE
// ─────────────────────────────────────────────────────────────

data class FloweringGene(
    val symbol: String,
    val rapId: String,
    val chromosome: String,
    val pathway: String,
    val function: String,
    val photoperiodSensitive: Boolean,
    val effect: String            // "promotes" or "represses" flowering
)

object FloweringGenes {

    val CATALOGUE = listOf(
        FloweringGene("Hd1",  "Os06g0275000", "chr06",
            "Photoperiod pathway",
            "CONSTANS ortholog - represses Hd3a under LD, promotes under SD",
            true, "dual"),
        FloweringGene("Hd3a", "Os06g0163700", "chr06",
            "Florigen pathway",
            "FT ortholog - mobile flowering signal, promotes heading under SD",
            true, "promotes"),
        FloweringGene("RFT1", "Os06g0242200", "chr06",
            "Florigen pathway",
            "Florigen 1 (FT-L) - promotes flowering under LD in temperate cultivars",
            true, "promotes"),
        FloweringGene("Ghd7", "Os07g0261200", "chr07",
            "Circadian / LD repressor",
            "CO-like protein - represses Hd3a under LD, major yield component",
            true, "represses"),
        FloweringGene("Ghd8/DTH8","Os08g0424300","chr08",
            "LD repressor",
            "HAP3 subunit - represses Hd3a under LD, increases grain number",
            true, "represses"),
        FloweringGene("OsPRR1/Toc1","Os02g0724000","chr02",
            "Circadian clock",
            "Pseudo-response regulator 1 - circadian clock component",
            true, "represses"),
        FloweringGene("OsPRR37/Hd2","Os07g0695100","chr07",
            "Circadian / photoperiod",
            "Pseudo-response regulator 37 - represses flowering under LD",
            true, "represses"),
        FloweringGene("OsGI",  "Os01g0830800", "chr01",
            "Photoperiod pathway",
            "GIGANTEA ortholog - regulates Hd1 and Hd3a",
            true, "promotes"),
        FloweringGene("SE5",   "Os06g0227700", "chr06",
            "Light signalling",
            "Heme oxygenase - phytochrome chromophore biosynthesis, early flowering",
            false, "promotes"),
        FloweringGene("Ehd1",  "Os10g0392400", "chr10",
            "SD promotion pathway",
            "B-type response regulator - promotes Hd3a/RFT1 under SD",
            true, "promotes"),
        FloweringGene("Ehd2/RID1","Os1g0182800","chr01",
            "Ehd1 regulation",
            "C2H2 zinc finger - represses Ehd1",
            true, "represses"),
        FloweringGene("Ehd3",  "Os01g0747900", "chr01",
            "Ehd1 regulation",
            "PHD finger protein - activates Ehd1 under SD",
            true, "promotes"),
        FloweringGene("Ehd4",  "Os01g0182800", "chr01",
            "Ehd1 regulation",
            "CCCH zinc finger - represses Ehd1 under LD",
            true, "represses"),
        FloweringGene("OsFKF1","Os05g0536000", "chr05",
            "Light-F-box pathway",
            "FKF1 ortholog - promotes Hd3a via Hd1 degradation under LD",
            true, "promotes"),
        FloweringGene("OsCOL4","Os03g0819800", "chr03",
            "CO-like repressor",
            "COL4 - represses heading under SD, promotes under LD",
            true, "represses"),
        FloweringGene("IDS1/OsAP2-39","Os04g0448900","chr04",
            "Autonomous pathway",
            "AP2-domain - induces heading, reduces spikelet number",
            false, "promotes"),
        FloweringGene("RCN1/OsAGL6","Os06g0193700","chr06",
            "MADS-box pathway",
            "AGAMOUS-LIKE 6 - controls heading and panicle development",
            false, "promotes"),
        FloweringGene("OsMADS51","Os05g0450480","chr05",
            "Autonomous/SD pathway",
            "MADS-box gene - activates Ehd1 under SD",
            true, "promotes")
    )

    val bySymbol: Map<String, FloweringGene>  by lazy { CATALOGUE.associateBy { it.symbol } }
    val byRapId: Map<String, FloweringGene>   by lazy { CATALOGUE.associateBy { it.rapId } }
    val byChrom: Map<String, List<FloweringGene>> by lazy { CATALOGUE.groupBy { it.chromosome } }

    fun getBySymbol(symbol: String) = bySymbol[symbol]
    fun getByRapId(rapId: String)   = byRapId[rapId]
    fun getByChromosome(chrom: String) = byChrom[RiceVariantUtils.normaliseChrom(chrom)] ?: emptyList()

    fun promoters(): List<FloweringGene> = CATALOGUE.filter { it.effect == "promotes" }
    fun repressors(): List<FloweringGene> = CATALOGUE.filter { it.effect == "represses" }
    fun photoperiodGenes(): List<FloweringGene> = CATALOGUE.filter { it.photoperiodSensitive }

    fun pathwaySummary(): String = buildString {
        append("=== Rice Flowering Time Gene Catalogue ===\n")
        val byPath = CATALOGUE.groupBy { it.pathway }
        byPath.forEach { (path, genes) ->
            append("\n$path (${genes.size} genes):\n")
            genes.forEach { g ->
                append("  ${g.symbol.padEnd(15)} ${g.rapId}  ${g.effect}  ")
                append(if (g.photoperiodSensitive) "[photoperiod]" else "[autonomous]")
                append("\n")
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
// RICE CHROMOSOME INFORMATION
// ─────────────────────────────────────────────────────────────

object RiceChromosomes {

    data class ChromInfo(
        val name: String,
        val length: Int,
        val centromereStart: Int,
        val centromereEnd: Int,
        val telomereSize: Int = 100_000
    ) {
        val armRatio: Double get() {
            val p = centromereStart.toDouble()
            val q = length - centromereEnd.toDouble()
            return if (q > 0) p / q else 0.0
        }
    }

    // IRGSP-1.0 Nipponbare reference coordinates
    val CHROMOSOMES = listOf(
        ChromInfo("chr01", 43_270_923, 16_589_552, 17_151_912),
        ChromInfo("chr02", 35_937_250, 10_000_000, 10_500_000),
        ChromInfo("chr03", 36_413_819, 11_000_000, 11_500_000),
        ChromInfo("chr04", 35_502_694,  2_000_000,  2_800_000),
        ChromInfo("chr05", 29_958_434,  9_000_000,  9_500_000),
        ChromInfo("chr06", 31_248_787, 10_000_000, 10_600_000),
        ChromInfo("chr07", 29_697_621,  7_000_000,  7_500_000),
        ChromInfo("chr08", 28_443_022,  8_000_000,  8_500_000),
        ChromInfo("chr09", 23_012_720,  6_000_000,  6_500_000),
        ChromInfo("chr10", 23_207_287,  8_000_000,  8_500_000),
        ChromInfo("chr11", 29_021_106,  9_000_000,  9_500_000),
        ChromInfo("chr12", 27_531_856,  8_000_000,  8_500_000)
    )

    val byName: Map<String, ChromInfo> by lazy {
        CHROMOSOMES.flatMap { c ->
            listOf(c.name, c.name.replace("chr0","chr").replace("chr",""),
                "Os${c.name.replace("chr","")}")
                .map { it to c }
        }.toMap()
    }

    fun get(chrom: String): ChromInfo? =
        byName[RiceVariantUtils.normaliseChrom(chrom)] ?: byName[chrom]

    fun inCentromere(chrom: String, pos: Int): Boolean {
        val info = get(chrom) ?: return false
        return pos in info.centromereStart..info.centromereEnd
    }

    fun inTelomere(chrom: String, pos: Int): Boolean {
        val info = get(chrom) ?: return false
        return pos <= info.telomereSize || pos >= info.length - info.telomereSize
    }

    fun position(chrom: String, pos: Int): String {
        return when {
            inTelomere(chrom, pos) -> "telomeric"
            inCentromere(chrom, pos) -> "centromeric"
            else -> {
                val info = get(chrom)
                if (info != null && pos < info.centromereStart) "short arm (p)"
                else "long arm (q)"
            }
        }
    }

    fun totalGenomeSize(): Long = CHROMOSOMES.map { it.length.toLong() }.sum()
}

// ─────────────────────────────────────────────────────────────
// RICE SUBPOPULATION DEFINITIONS
// ─────────────────────────────────────────────────────────────

object RiceSubpopulations {

    data class Subpopulation(
        val name: String,
        val code: String,
        val ecotype: String,
        val origin: String,
        val traits: List<String>
    )

    val SUBPOPS = listOf(
        Subpopulation("Aus", "aus", "tropical",
            "South Asia (Bangladesh, India)",
            listOf("Early heading", "Drought tolerance", "High protein content")),
        Subpopulation("Indica", "ind", "tropical",
            "Tropical Asia (India, SE Asia, China south)",
            listOf("Long grain", "High amylose", "Photoperiod sensitive")),
        Subpopulation("Aromatic", "aro", "tropical",
            "South Asia (Basmati, Jasmine types)",
            listOf("Fragrance (2-AP)", "Long slender grain", "Low amylose")),
        Subpopulation("Temperate Japonica", "trop_jap",
            "temperate",
            "East Asia (China north, Japan, Korea)",
            listOf("Short grain", "Low amylose", "Cold tolerance")),
        Subpopulation("Tropical Japonica", "trj", "tropical",
            "SE Asia, Latin America, Africa",
            listOf("Wide adaptation", "Photoperiod insensitive", "Blast resistance")),
        Subpopulation("Admixed", "adm", "mixed",
            "Various",
            listOf("Intermediate phenotypes", "Breeding materials"))
    )

    val byCode: Map<String, Subpopulation> by lazy { SUBPOPS.associateBy { it.code } }
    val byName: Map<String, Subpopulation> by lazy { SUBPOPS.associateBy { it.name } }
}
