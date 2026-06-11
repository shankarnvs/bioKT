package biokt.rice

import biokt.*
import kotlin.math.*

// ============================================================
// GENOME ANNOTATION — GFF3 / GTF parsing
// Subpackage: biokt.rice
//
// Handles IRGSP-1.0 / MSU7 / RAP-DB annotation formats.
// Supports both RAP-DB (Os01g0100100) and MSU (LOC_Os01g01010)
// gene ID schemes, and cross-maps between them where possible.
// ============================================================

// ─────────────────────────────────────────────────────────────
// ANNOTATION FEATURE
// ─────────────────────────────────────────────────────────────

data class GffFeature(
    val seqname: String,
    val source: String,
    val type: String,
    val start: Int,          // 1-based, inclusive
    val end: Int,            // 1-based, inclusive
    val score: Double?,
    val strand: Char,        // '+', '-', '.'
    val phase: Int?,         // 0, 1, 2, or null
    val attributes: Map<String, List<String>>
) {
    val length: Int get() = end - start + 1

    fun attr(key: String): String? = attributes[key]?.firstOrNull()
    fun attrList(key: String): List<String> = attributes[key] ?: emptyList()

    val id: String? get() = attr("ID")
    val name: String? get() = attr("Name") ?: attr("gene_name") ?: attr("gene_id")
    val parentId: String? get() = attr("Parent")

    fun isGene(): Boolean = type.toLowerCase() in setOf("gene", "pseudogene")
    fun isMRNA(): Boolean = type.toLowerCase() in setOf("mrna", "transcript", "processed_transcript")
    fun isCDS(): Boolean = type.toLowerCase() == "cds"
    fun isExon(): Boolean = type.toLowerCase() == "exon"
    fun isUTR(): Boolean = type.toLowerCase() in setOf("utr", "five_prime_utr", "three_prime_utr",
        "5utr", "3utr", "5'utr", "3'utr")

    /** Return true if this feature overlaps [queryStart, queryEnd] on the same sequence */
    fun overlaps(seq: String, queryStart: Int, queryEnd: Int): Boolean =
        seqname == seq && start <= queryEnd && end >= queryStart

    /** GFF3 formatted output line */
    fun toGff3(): String {
        val scoreStr = score?.toString() ?: "."
        val phaseStr = phase?.toString() ?: "."
        val attrStr  = attributes.entries.joinToString(";") { (k, vs) ->
            vs.joinToString(",") { "$k=${gff3Escape(it)}" }
        }
        return listOf(seqname, source, type, start, end, scoreStr, strand, phaseStr, attrStr)
            .joinToString("\t")
    }

    private fun gff3Escape(s: String) = s.replace(";", "%3B").replace(",", "%2C")
        .replace("=", "%3D").replace("&", "%26").replace(" ", "%20")
}

// ─────────────────────────────────────────────────────────────
// GENE MODEL — hierarchical gene structure
// ─────────────────────────────────────────────────────────────

data class GeneModel(
    val gene: GffFeature,
    val transcripts: List<TranscriptModel>
) {
    val id: String get() = gene.id ?: gene.name ?: ""
    val chrom: String get() = gene.seqname
    val start: Int get() = gene.start
    val end: Int get() = gene.end
    val strand: Char get() = gene.strand
    val length: Int get() = gene.length

    val rapId: String? get() = extractRapId(id) ?: transcripts.mapNotNull { extractRapId(it.mrna.id ?: "") }.firstOrNull()
    val msuId: String? get() = extractMsuId(id) ?: transcripts.mapNotNull { extractMsuId(it.mrna.id ?: "") }.firstOrNull()

    val primaryTranscript: TranscriptModel? get() = transcripts.maxBy { it.cdsLength }

    /** Upstream promoter region coordinates */
    fun promoterRegion(upstreamBp: Int = 2000, downstreamBp: Int = 200): Triple<String, Int, Int> {
        return when (strand) {
            '+' -> Triple(chrom, maxOf(1, start - upstreamBp), start + downstreamBp)
            '-' -> Triple(chrom, maxOf(1, end - downstreamBp), end + upstreamBp)
            else -> Triple(chrom, maxOf(1, start - upstreamBp), start + downstreamBp)
        }
    }

    fun summary(): String = buildString {
        append("Gene: $id  ${chrom}:${start}-${end} ($strand)\n")
        append("  RAP-DB: ${rapId ?: "—"}  MSU: ${msuId ?: "—"}\n")
        append("  Length: ${length} bp  Transcripts: ${transcripts.size}\n")
        primaryTranscript?.let { t ->
            append("  Primary transcript: ${t.id}\n")
            append("  Exons: ${t.exons.size}  CDS length: ${t.cdsLength} bp\n")
        }
    }
}

data class TranscriptModel(
    val mrna: GffFeature,
    val exons: List<GffFeature>,
    val cdss: List<GffFeature>,
    val utrs: List<GffFeature>
) {
    val id: String get() = mrna.id ?: ""
    val strand: Char get() = mrna.strand

    val cdsLength: Int get() = cdss.map { it.length }.sum()
    val exonLength: Int get() = exons.map { it.length }.sum()
    val introns: List<Pair<Int, Int>> get() {
        if (exons.size < 2) return emptyList()
        val sorted = exons.sortedBy { it.start }
        return sorted.zipWithNext { a, b -> Pair(a.end + 1, b.start - 1) }
    }
    val spliceSites: List<Pair<Int, Char>> get() =
        introns.flatMap { (s, e) -> listOf(s to strand, e to strand) }
}

// ─────────────────────────────────────────────────────────────
// GFF3 / GTF PARSER
// ─────────────────────────────────────────────────────────────

object AnnotationParser {

    /** Parse a GFF3 file */
    fun parseGff3(path: String): AnnotationDb =
        parseGff3Lines(java.io.File(path).readLines())

    /** Parse GFF3 from String */
    fun parseGff3String(text: String): AnnotationDb =
        parseGff3Lines(text.lines())

    /** Parse a GTF file (Ensembl/GENCODE format) */
    fun parseGtf(path: String): AnnotationDb =
        parseGtfLines(java.io.File(path).readLines())

    /** Parse GTF from String */
    fun parseGtfString(text: String): AnnotationDb =
        parseGtfLines(text.lines())

    private fun parseGff3Lines(lines: List<String>): AnnotationDb {
        val features = mutableListOf<GffFeature>()
        for (line in lines) {
            if (line.startsWith("#") || line.isBlank()) continue
            val f = parseGff3Line(line) ?: continue
            features.add(f)
        }
        return buildAnnotationDb(features)
    }

    private fun parseGtfLines(lines: List<String>): AnnotationDb {
        val features = mutableListOf<GffFeature>()
        for (line in lines) {
            if (line.startsWith("#") || line.isBlank()) continue
            val f = parseGtfLine(line) ?: continue
            features.add(f)
        }
        return buildAnnotationDb(features)
    }

    private fun parseGff3Line(line: String): GffFeature? {
        val cols = line.split("\t")
        if (cols.size < 9) return null
        return try {
            val attrs = parseGff3Attrs(cols[8])
            GffFeature(
                seqname    = cols[0],
                source     = cols[1],
                type       = cols[2],
                start      = cols[3].toInt(),
                end        = cols[4].toInt(),
                score      = cols[5].takeIf { it != "." }?.toDoubleOrNull(),
                strand     = cols[6].firstOrNull() ?: '.',
                phase      = cols[7].toIntOrNull(),
                attributes = attrs
            )
        } catch (e: Exception) { null }
    }

    private fun parseGff3Attrs(raw: String): Map<String, List<String>> {
        val attrs = mutableMapOf<String, MutableList<String>>()
        raw.split(";").forEach { pair ->
            val kv = pair.trim().split("=", limit = 2)
            if (kv.size == 2) {
                val key  = kv[0].trim()
                val vals = kv[1].split(",").map { gff3Unescape(it.trim()) }
                attrs.getOrPut(key) { mutableListOf() }.addAll(vals)
            }
        }
        return attrs
    }

    private fun gff3Unescape(s: String) = s
        .replace("%3B", ";").replace("%2C", ",").replace("%3D", "=")
        .replace("%26", "&").replace("%20", " ").replace("%09", "\t")

    private fun parseGtfLine(line: String): GffFeature? {
        val cols = line.split("\t")
        if (cols.size < 9) return null
        return try {
            val attrs = parseGtfAttrs(cols[8])
            GffFeature(
                seqname    = cols[0],
                source     = cols[1],
                type       = cols[2],
                start      = cols[3].toInt(),
                end        = cols[4].toInt(),
                score      = cols[5].takeIf { it != "." }?.toDoubleOrNull(),
                strand     = cols[6].firstOrNull() ?: '.',
                phase      = cols[7].toIntOrNull(),
                attributes = attrs
            )
        } catch (e: Exception) { null }
    }

    private fun parseGtfAttrs(raw: String): Map<String, List<String>> {
        val attrs = mutableMapOf<String, MutableList<String>>()
        val pattern = Regex("""(\w+)\s+"([^"]+)"""")
        pattern.findAll(raw).forEach { m ->
            attrs.getOrPut(m.groupValues[1]) { mutableListOf() }.add(m.groupValues[2])
        }
        // GTF uses gene_id/transcript_id → map to ID/Parent for GFF3 compatibility
        attrs["gene_id"]?.firstOrNull()?.let { attrs.getOrPut("ID") { mutableListOf() }.add(it) }
        attrs["transcript_id"]?.firstOrNull()?.let {
            attrs.getOrPut("ID") { mutableListOf() }.clear()
            attrs.getOrPut("ID") { mutableListOf() }.add(it)
        }
        return attrs
    }

    private fun buildAnnotationDb(features: List<GffFeature>): AnnotationDb {
        // Index by ID
        val byId = features.filter { it.id != null }.associateBy { it.id!! }

        // Build hierarchy: gene → mRNA → exon/CDS/UTR
        val genes = mutableListOf<GeneModel>()
        val geneFeatures = features.filter { it.isGene() }

        geneFeatures.forEach { gene ->
            val geneId = gene.id ?: return@forEach
            // Find transcripts whose Parent = geneId
            val txFeatures = features.filter { it.isMRNA() && it.parentId == geneId }
            val transcripts = txFeatures.map { tx ->
                val txId = tx.id ?: ""
                val exons  = features.filter { it.isExon() && it.parentId == txId }
                    .sortedBy { if (tx.strand == '-') -it.start else it.start }
                val cdss   = features.filter { it.isCDS() && it.parentId == txId }
                    .sortedBy { if (tx.strand == '-') -it.start else it.start }
                val utrs   = features.filter { it.isUTR() && it.parentId == txId }
                TranscriptModel(tx, exons, cdss, utrs)
            }
            genes.add(GeneModel(gene, transcripts))
        }

        // Handle flat GFF3 (no gene features, just mRNA)
        if (genes.isEmpty()) {
            features.filter { it.isMRNA() }.forEach { tx ->
                val txId = tx.id ?: ""
                val exons  = features.filter { it.isExon() && it.parentId == txId }
                val cdss   = features.filter { it.isCDS() && it.parentId == txId }
                val utrs   = features.filter { it.isUTR() && it.parentId == txId }
                val syntheticGene = tx.copy(type = "gene", attributes =
                    tx.attributes + mapOf("ID" to listOf("gene_$txId")))
                genes.add(GeneModel(syntheticGene, listOf(TranscriptModel(tx, exons, cdss, utrs))))
            }
        }

        return AnnotationDb(genes, features)
    }
}

// ─────────────────────────────────────────────────────────────
// ANNOTATION DATABASE — query interface
// ─────────────────────────────────────────────────────────────

class AnnotationDb(
    val genes: List<GeneModel>,
    val allFeatures: List<GffFeature>
) {
    val size: Int get() = genes.size

    private val byId: Map<String, GeneModel>   by lazy { genes.associateBy { it.id } }
    private val byRap: Map<String, GeneModel>  by lazy {
        genes.mapNotNull { g -> g.rapId?.let { it to g } }.toMap()
    }
    private val byMsu: Map<String, GeneModel>  by lazy {
        genes.mapNotNull { g -> g.msuId?.let { it to g } }.toMap()
    }
    private val byChrom: Map<String, List<GeneModel>> by lazy {
        genes.groupBy { it.chrom }
    }

    // ── Lookup by ID ─────────────────────────────────────────

    fun getById(id: String): GeneModel? = byId[id] ?: byRap[id] ?: byMsu[id]

    fun getByRapId(rapId: String): GeneModel? = byRap[rapId]
    fun getByMsuId(msuId: String): GeneModel? = byMsu[msuId]

    /** Look up a gene by partial ID match */
    fun search(query: String): List<GeneModel> =
        genes.filter { g ->
            g.id.contains(query, ignoreCase = true) ||
            g.rapId?.contains(query, ignoreCase = true) == true ||
            g.msuId?.contains(query, ignoreCase = true) == true
        }

    // ── Region queries ────────────────────────────────────────

    fun genesInRegion(chrom: String, start: Int, end: Int): List<GeneModel> =
        byChrom[chrom]?.filter { it.start <= end && it.end >= start } ?: emptyList()

    fun genesOnChrom(chrom: String): List<GeneModel> =
        byChrom[chrom]?.sortedBy { it.start } ?: emptyList()

    fun nearestGene(chrom: String, pos: Int): GeneModel? =
        byChrom[chrom]?.minBy { minOf(abs(it.start - pos), abs(it.end - pos)) }

    fun featuresInRegion(chrom: String, start: Int, end: Int, type: String? = null): List<GffFeature> =
        allFeatures.filter { f ->
            f.seqname == chrom && f.start <= end && f.end >= start &&
            (type == null || f.type.toLowerCase() == type.toLowerCase())
        }

    // ── Statistics ────────────────────────────────────────────

    fun summary(): String = buildString {
        append("=== Annotation Database ===\n")
        append("Genes            : ${genes.size}\n")
        append("Chromosomes      : ${byChrom.keys.sorted().joinToString(", ")}\n")
        val withRap = genes.count { it.rapId != null }
        val withMsu = genes.count { it.msuId != null }
        append("With RAP-DB IDs  : $withRap\n")
        append("With MSU IDs     : $withMsu\n")
        val totalTx = genes.map { it.transcripts.size }.sum()
        append("Total transcripts: $totalTx\n")
        val avgTx = if (genes.isNotEmpty()) totalTx.toDouble() / genes.size else 0.0
        append("Avg tx/gene      : ${"%.2f".format(avgTx)}\n")
        val intronCounts = genes.flatMap { it.transcripts }.flatMap { it.introns }
        append("Total introns    : ${intronCounts.size}\n")
    }

    // ── Promoter extraction ───────────────────────────────────

    /** Extract promoter region coordinates for all genes */
    fun promoterRegions(upstreamBp: Int = 2000, downstreamBp: Int = 200): List<Triple<String, Int, Int>> =
        genes.map { it.promoterRegion(upstreamBp, downstreamBp) }

    /** Get splice site positions for all genes */
    fun spliceSites(): List<Triple<String, Int, Char>> =
        genes.flatMap { gene ->
            gene.transcripts.flatMap { tx ->
                tx.spliceSites.map { (pos, strand) -> Triple(gene.chrom, pos, strand) }
            }
        }.distinct()

    // ── GFF3 output ───────────────────────────────────────────

    fun writeGff3(path: String) {
        java.io.File(path).bufferedWriter().use { w ->
            w.write("##gff-version 3\n")
            allFeatures.sortedWith(compareBy({ it.seqname }, { it.start }))
                .forEach { w.write(it.toGff3() + "\n") }
        }
    }
}

// ─────────────────────────────────────────────────────────────
// RICE ID UTILITIES
// ─────────────────────────────────────────────────────────────

private val RAP_PATTERN = Regex("""Os(\d{2})g(\d{7})""", RegexOption.IGNORE_CASE)
private val MSU_PATTERN = Regex("""LOC_Os(\d{2})g(\d{5})""", RegexOption.IGNORE_CASE)

fun extractRapId(s: String): String? = RAP_PATTERN.find(s)?.value?.let {
    "Os${it.substring(2, 4).trimStart('0').padStart(2, '0')}g${it.substring(5)}"
}

fun extractMsuId(s: String): String? = MSU_PATTERN.find(s)?.value

object RiceGeneIds {

    /** Parse a RAP-DB ID like Os01g0100100 → chrom=1, locus=100100 */
    fun parseRapId(id: String): Pair<Int, Int>? {
        val m = RAP_PATTERN.find(id) ?: return null
        val chrom = m.groupValues[1].toIntOrNull() ?: return null
        val locus = m.groupValues[2].toIntOrNull() ?: return null
        return chrom to locus
    }

    /** Parse an MSU ID like LOC_Os01g01010 → chrom=1, locus=1010 */
    fun parseMsuId(id: String): Pair<Int, Int>? {
        val m = MSU_PATTERN.find(id) ?: return null
        val chrom = m.groupValues[1].toIntOrNull() ?: return null
        val locus = m.groupValues[2].toIntOrNull() ?: return null
        return chrom to locus
    }

    /** Sort gene IDs by chromosomal position (RAP-DB format) */
    fun sortByPosition(ids: List<String>): List<String> =
        ids.sortedWith(Comparator { a, b ->
            val pA = parseRapId(a) ?: parseMsuId(a) ?: return@Comparator a.compareTo(b)
            val pB = parseRapId(b) ?: parseMsuId(b) ?: return@Comparator a.compareTo(b)
            if (pA.first != pB.first) pA.first - pB.first else pA.second - pB.second
        })

    /** Guess chromosome from a RAP or MSU gene ID */
    fun chromFromId(id: String): String? {
        val chrom = parseRapId(id)?.first ?: parseMsuId(id)?.first ?: return null
        return "chr${chrom.toString().padStart(2, '0')}"
    }

    /** Built-in catalogue of key rice agronomic genes with RAP IDs */
    val KNOWN_GENES = mapOf(
        // Grain quality
        "Os06g0133000" to GeneInfo("Wx",   "Waxy",   "Starch synthesis — glutinous phenotype"),
        "Os07g0695100" to GeneInfo("OsSSI","Starch synthase I","Starch quality"),
        // Heading date / flowering
        "Os06g0275000" to GeneInfo("Hd1",  "Heading date 1", "CONSTANS homolog — photoperiod flowering"),
        "Os06g0163700" to GeneInfo("Hd3a", "Heading date 3a","Florigen — day-length insensitive"),
        "Os06g0242200" to GeneInfo("RFT1", "RFT1",   "Florigen — promotes flowering under LD"),
        "Os07g0261200" to GeneInfo("Ghd7", "Ghd7",   "Grain number, plant height, heading date"),
        "Os08g0424300" to GeneInfo("OsDTH8","DTH8",  "Days-to-heading under short days"),
        // Stress resistance
        "Os11g0508800" to GeneInfo("Pita",  "Pi-ta",  "Blast resistance gene (NBS-LRR)"),
        "Os11g0695900" to GeneInfo("Pid2",  "Pi-d2",  "Blast resistance — B-lectin kinase"),
        "Os04g0448900" to GeneInfo("OsSULT","OsSULT", "Sulphate transporter"),
        "Os04g0448400" to GeneInfo("Sub1A", "Submergence 1A","Submergence tolerance — ERF"),
        "Os11g0615700" to GeneInfo("OsNramp5","OsNramp5","Cadmium uptake — transporter"),
        // Domestication
        "Os04g0670800" to GeneInfo("sh4",  "Shattering 4","Grain shattering control"),
        "Os07g0153600" to GeneInfo("Prog1","PROG1",   "Prostrate growth 1 — domestication"),
        "Os09g0310900" to GeneInfo("Rc",   "Rc",     "Red coleoptile / pericarp colour"),
        // Grain size / yield
        "Os03g0407400" to GeneInfo("GS3",  "Grain size 3","Major grain length QTL"),
        "Os05g0187500" to GeneInfo("GW5",  "Grain width 5","Grain width and weight"),
        "Os02g0814300" to GeneInfo("GW2",  "Grain width 2","RING-type E3 ubiquitin ligase"),
        "Os07g0603300" to GeneInfo("GIF1", "GIF1",   "Cell wall invertase — grain filling"),
        "Os01g0885100" to GeneInfo("OsLOGL7","LOG7", "Cytokinin activation — tiller number")
    )

    data class GeneInfo(val symbol: String, val fullName: String, val function: String)

    fun lookupGene(rapId: String): GeneInfo? = KNOWN_GENES[rapId]

    fun lookupBySymbol(symbol: String): Pair<String, GeneInfo>? =
        KNOWN_GENES.entries.firstOrNull { it.value.symbol.equals(symbol, ignoreCase = true) }
            ?.let { it.key to it.value }
}

// ─────────────────────────────────────────────────────────────
// ANNOTATION × VCF OVERLAP
// ─────────────────────────────────────────────────────────────

object AnnotationOverlap {

    data class VariantAnnotation(
        val record: VcfRecord,
        val geneId: String,
        val geneSymbol: String?,
        val featureType: String,   // CDS, exon, UTR, intron, upstream, downstream
        val consequence: String    // missense, synonymous, stop_gain, splice, intergenic
    )

    /** Annotate VCF records against a gene annotation database */
    fun annotateVariants(vcf: VcfFile, db: AnnotationDb): List<VariantAnnotation> {
        return vcf.records.flatMap { rec ->
            val overlapping = db.genesInRegion(rec.chrom, rec.pos, rec.pos)
            if (overlapping.isEmpty()) {
                // Check 2 kb upstream/downstream
                val nearest = db.nearestGene(rec.chrom, rec.pos)
                if (nearest != null) {
                    val dist = minOf(abs(nearest.start - rec.pos), abs(nearest.end - rec.pos))
                    if (dist <= 2000) {
                        val ftype = if ((nearest.strand == '+' && rec.pos < nearest.start) ||
                                       (nearest.strand == '-' && rec.pos > nearest.end)) "upstream"
                                    else "downstream"
                        listOf(VariantAnnotation(rec, nearest.id, nearest.gene.attr("Name"),
                            ftype, "regulatory_region"))
                    } else emptyList()
                } else emptyList()
            } else {
                overlapping.map { gene ->
                    val ftype = determineFeatureType(rec, gene)
                    val csq   = estimateConsequence(rec, ftype)
                    VariantAnnotation(rec, gene.id, gene.gene.attr("Name"), ftype, csq)
                }
            }
        }
    }

    private fun determineFeatureType(rec: VcfRecord, gene: GeneModel): String {
        val tx = gene.primaryTranscript ?: return "intron"
        return when {
            tx.cdss.any { rec.pos in it.start..it.end } -> "CDS"
            tx.exons.any { rec.pos in it.start..it.end } -> "exon"
            tx.utrs.any  { rec.pos in it.start..it.end } -> "UTR"
            else -> "intron"
        }
    }

    private fun estimateConsequence(rec: VcfRecord, featureType: String): String {
        if (!rec.isSNP) return when (featureType) {
            "CDS"   -> "frameshift"
            "exon"  -> "exon_indel"
            else    -> "indel"
        }
        return when (featureType) {
            "CDS"        -> "missense_or_synonymous"   // would need codon context for exact class
            "exon"       -> "exon_SNP"
            "UTR"        -> "UTR_SNP"
            "intron"     -> if (isNearSpliceSite(rec)) "splice_region" else "intronic"
            "upstream"   -> "upstream_variant"
            "downstream" -> "downstream_variant"
            else         -> "unknown"
        }
    }

    private fun isNearSpliceSite(rec: VcfRecord): Boolean {
        // Canonical splice site positions ±2 bp from exon boundary
        // Without exon coords here we use a proxy: intronic variants close to round numbers
        return false  // placeholder — full implementation requires exon coords
    }
}
