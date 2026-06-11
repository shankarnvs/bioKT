package biokt.rice

// ============================================================
// biokt.rice — Rice Genomics Subpackage
// ============================================================
//
// A specialised extension of BioKt for rice (Oryza sativa)
// genomics research. Builds on top of the core biokt package.
//
// Modules:
//   RiceVariants.kt     — VCF/BCF parsing, SNP/INDEL handling,
//                         LD calculation, allele frequencies,
//                         genotype matrices, population stats
//
//   GenomeAnnotation.kt — GFF3/GTF parsing, gene model hierarchy,
//                         RAP-DB + MSU ID handling,
//                         promoter extraction, variant annotation
//
//   RiceGenomics.kt     — SSR/microsatellite detection,
//                         cis-element scanner (W-box, ABRE, etc.),
//                         TE classification (Gypsy, Copia, MuDR…),
//                         R-gene domain patterns,
//                         heading date gene catalogue,
//                         chromosome coordinates (IRGSP-1.0),
//                         subpopulation definitions
//
//   RiceExpression.kt   — TPM/FPKM/RPKM normalisation,
//                         DESeq2-style size factor normalisation,
//                         differential expression (Welch t-test + BH),
//                         co-expression (Pearson r),
//                         stress gene sets (drought/cold/blast/…),
//                         gene set enrichment, tissue specificity (τ)
//
//   HaplotypeAnalysis.kt — Gabriel haplotype blocks,
//                          windowed π / θ_W / Tajima's D,
//                          selection sweep detection (XP-CLR proxy),
//                          haplotype network
//
// Usage:
//   import biokt.rice.*   — imports all rice subpackage
//   import biokt.*        — imports core BioKt classes
//
// All rice subpackage classes use core BioKt types (DNASequence,
// ProteinSequence, etc.) from the biokt package.
// ============================================================

// ─────────────────────────────────────────────────────────────
// PACKAGE-LEVEL CONVENIENCE FUNCTIONS
// ─────────────────────────────────────────────────────────────

import biokt.*

/** Quick summary of a rice genome annotation database */
fun AnnotationDb.riceStats(): String = buildString {
    append(summary())
    val knownGenes = genes.count { g ->
        RiceGeneIds.KNOWN_GENES.containsKey(g.rapId)
    }
    append("Known agronomic genes: $knownGenes\n")
    val floweringGenes = genes.count { g ->
        FloweringGenes.byRapId.containsKey(g.rapId)
    }
    append("Flowering time genes: $floweringGenes\n")
}

/** Scan a promoter sequence for rice cis-elements and return summary */
fun DNASequence.scanRicePromoter(): Map<String, Int> {
    val results = CisElementScanner.scan(this)
    return CisElementScanner.summarise(results)
}

/** Detect SSRs in a rice sequence */
fun DNASequence.detectSSRs(motifLengths: IntRange = 1..6): List<SSRMarker> =
    SSRDetector.detect(this, motifLengths)

/** Check if this protein has R-gene characteristics */
fun biokt.ProteinSequence.isRGene(): Pair<Boolean, List<String>> =
    RGeneClassifier.isLikelyRGene(this)

/** Classify TE content of a rice sequence */
fun DNASequence.estimateTEContent(): Double =
    TEClassifier.estimateTEContent(this)

/** Look up this gene in the rice agronomic gene database */
fun String.riceGeneInfo(): RiceGeneIds.GeneInfo? = RiceGeneIds.lookupGene(this)

/** Parse chromosome name to normalised rice format (chr01..chr12) */
fun String.normaliseRiceChrom(): String = RiceVariantUtils.normaliseChrom(this)

/** Check which stress gene sets contain this RAP-DB gene */
fun String.stressGeneSets(): List<RiceStressGeneSets.GeneSet> =
    RiceStressGeneSets.setsContaining(this)

// ─────────────────────────────────────────────────────────────
// RICE GENOME CONSTANTS
// ─────────────────────────────────────────────────────────────

object Rice {
    /** Reference genome version */
    const val REFERENCE = "IRGSP-1.0 (Oryza sativa ssp. japonica cv. Nipponbare)"

    /** Number of chromosomes */
    const val NUM_CHROMOSOMES = 12

    /** Approximate genome size */
    const val GENOME_SIZE_BP = 389_840_340L

    /** Approximate gene count */
    const val APPROX_GENE_COUNT = 35_679

    /** TE fraction */
    const val TE_FRACTION = 0.45

    /** Average gene density (genes per Mb) */
    const val GENE_DENSITY = 91.5

    /** Species */
    const val SPECIES = "Oryza sativa L."
    const val SUBSPECIES_JAPONICA = "Oryza sativa ssp. japonica"
    const val SUBSPECIES_INDICA   = "Oryza sativa ssp. indica"

    /** Key databases */
    const val RAPDB_URL   = "https://rapdb.dna.affrc.go.jp"
    const val MSU_URL     = "https://rice.uga.edu"
    const val RICEXPRO_URL = "https://ricexpro.dna.affrc.go.jp"
    const val IRRI_URL    = "https://www.irri.org"
    const val RICEVARMAP_URL = "http://ricevarmap.ncpgr.cn"
}
