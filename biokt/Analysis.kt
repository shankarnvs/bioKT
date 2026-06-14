package biokt

import kotlin.math.*

// ─── Restriction Enzymes ────────────────────────────────────

data class RestrictionSite(
    val enzyme: String,
    val position: Int,
    val recognition: String,
    val cutPosition: Int   // cut position within recognition site
)

object RestrictionEnzymes {

    // Recognition site → (enzyme name, cut offset from start)
    val ENZYMES: Map<String, Pair<String, Int>> = mapOf(
        "GAATTC"  to ("EcoRI"   to 1),
        "AAGCTT"  to ("HindIII" to 1),
        "GGATCC"  to ("BamHI"   to 1),
        "CTCGAG"  to ("XhoI"    to 1),
        "GCGGCCGC" to ("NotI"   to 2),
        "GTCGAC"  to ("SalI"    to 1),
        "CCCGGG"  to ("SmaI"    to 3),
        "ACGCGT"  to ("MluI"    to 1),
        "CCATGG"  to ("NcoI"    to 1),
        "CATATG"  to ("NdeI"    to 2),
        "TCTAGA"  to ("XbaI"    to 1),
        "AGATCT"  to ("BglII"   to 1),
        "GGTACC"  to ("KpnI"    to 5),
        "GAGCTC"  to ("SacI"    to 5),
        "CAATTG"  to ("MfeI"    to 1),
        "TTAATTAA" to ("PacI"   to 3),
        "GCATGC"  to ("SphI"    to 5),
        "ACTAGT"  to ("SpeI"    to 1),
        "GGATCN"  to ("BamHI-like" to 1),
        "AGGCCT"  to ("StuI"    to 3),
        "ATGCAT"  to ("NsiI"    to 5),
        "GGCGCGCC" to ("AscI"  to 2),
        "TTTTAAAA" to ("SwaI"  to 4),
        "TAATTAA" to ("MssI"   to 3),
        "RGATCY"  to ("BclI-like" to 1)
    )

    fun digest(dna: DNASequence): List<RestrictionSite> {
        val sites = mutableListOf<RestrictionSite>()
        val seq = dna.sequence
        for ((recognition, info) in ENZYMES) {
            val (name, cutOffset) = info
            if (recognition.contains('N') || recognition.contains('R') || recognition.contains('Y')) continue
            var idx = 0
            while (true) {
                val found = seq.indexOf(recognition, idx)
                if (found == -1) break
                sites.add(RestrictionSite(name, found, recognition, found + cutOffset))
                idx = found + 1
            }
        }
        return sites.sortedBy { it.position }
    }

    fun virtualDigest(dna: DNASequence, enzymes: List<String>): List<DNASequence> {
        val allSites = digest(dna).filter { it.enzyme in enzymes }
        if (allSites.isEmpty()) return listOf(dna)

        val cutPositions = allSites.map { it.cutPosition }.sorted()
        val fragments = mutableListOf<DNASequence>()
        var start = 0
        for (cut in cutPositions) {
            if (cut > start && cut <= dna.length)
                fragments.add(dna.subSequence(start, cut))
            start = cut
        }
        if (start < dna.length) fragments.add(dna.subSequence(start, dna.length))
        return fragments
    }
}

// ─── Sequence Statistics ────────────────────────────────────

object SeqStats {

    data class NStats(
        val n50: Int, val n90: Int, val l50: Int, val l90: Int,
        val totalLength: Long, val numSeqs: Int,
        val maxLen: Int, val minLen: Int, val meanLen: Double
    )

    fun nStats(sequences: List<BioSequence>): NStats {
        val lengths = sequences.map { it.length }.sortedDescending()
        val total = lengths.map { it.toLong() }.sum()
        var cumLen = 0L; var n50 = 0; var n90 = 0; var l50 = 0; var l90 = 0
        var foundN50 = false; var foundN90 = false
        for ((i, len) in lengths.withIndex()) {
            cumLen += len
            if (!foundN50 && cumLen >= total * 0.5) { n50 = len; l50 = i + 1; foundN50 = true }
            if (!foundN90 && cumLen >= total * 0.9) { n90 = len; l90 = i + 1; foundN90 = true }
        }
        return NStats(n50, n90, l50, l90, total, lengths.size,
            lengths.max() ?: 0, lengths.min() ?: 0,
            if (lengths.isEmpty()) 0.0 else total.toDouble() / lengths.size)
    }

    fun codonUsageTable(dna: DNASequence): Map<String, Int> {
        val counts = mutableMapOf<String, Int>()
        val seq = dna.sequence
        var i = 0
        while (i + 3 <= seq.length) {
            val codon = seq.substring(i, i + 3)
            counts[codon] = (counts[codon] ?: 0) + 1
            i += 3
        }
        return counts.toSortedMap()
    }

    fun relativeCodonAdaptation(codonCounts: Map<String, Int>, table: CodonTable = CodonTables.STANDARD): Map<String, Double> {
        // Group synonymous codons
        val synonymGroups = mutableMapOf<Char, MutableList<String>>()
        val allCodons = listOf("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG",
                               "ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG",
                               "TCT","TCC","TCA","TCG","AGT","AGC","CCT","CCC","CCA","CCG",
                               "ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG",
                               "TAT","TAC","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG",
                               "GAT","GAC","GAA","GAG","TGT","TGC","TGG","CGT","CGC","CGA","CGG",
                               "AGA","AGG","GGT","GGC","GGA","GGG")
        for (codon in allCodons) {
            val aa = try { table.translate(codon) } catch (e: Exception) { continue }
            if (aa != '*') synonymGroups.getOrPut(aa) { mutableListOf() }.add(codon)
        }

        val rcsu = mutableMapOf<String, Double>()
        for ((_, codons) in synonymGroups) {
            val total = codons.map { codonCounts[it] ?: 0 }.sum().toDouble()
            val expected = total / codons.size
            for (codon in codons) {
                val observed = (codonCounts[codon] ?: 0).toDouble()
                rcsu[codon] = if (expected == 0.0) 0.0 else observed / expected
            }
        }
        return rcsu
    }

    fun nucleotideFrequency(seq: BioSequence): Map<Char, Double> {
        val counts = seq.countAll()
        val total = counts.values.sum().toDouble()
        return counts.mapValues { (_, v) -> v / total }
    }

    fun dinucleotideFrequency(seq: BioSequence): Map<String, Double> {
        val counts = mutableMapOf<String, Int>()
        for (i in 0 until seq.length - 1) {
            val di = "${seq[i]}${seq[i+1]}"
            counts[di] = (counts[di] ?: 0) + 1
        }
        val total = counts.values.sum().toDouble()
        return counts.mapValues { (_, v) -> v / total }
    }

    // Shannon entropy
    fun entropy(sequence: BioSequence): Double {
        val freqs = nucleotideFrequency(sequence)
        return -freqs.values.filter { it > 0 }.map { it * log2(it) }.sum()
    }

    // Linguistic complexity
    fun linguisticComplexity(sequence: BioSequence): Double {
        val seq = sequence.sequence
        var totalObserved = 0; var totalPossible = 0
        for (k in 1..minOf(seq.length, 8)) {
            val observed = (0..seq.length - k).map { seq.substring(it, it + k) }.toSet().size
            val possible = minOf(4.0.pow(k).toInt(), seq.length - k + 1)
            totalObserved += observed; totalPossible += possible
        }
        return if (totalPossible == 0) 0.0 else totalObserved.toDouble() / totalPossible
    }
}

// ─── Population Genetics ────────────────────────────────────

object PopGen {

    // Tajima's D
    fun tajimasD(sequences: List<DNASequence>): Double {
        val n = sequences.size
        if (n < 4) throw IllegalArgumentException("Need at least 4 sequences for Tajima's D")
        val len = sequences[0].length

        // Count segregating sites
        var S = 0
        for (col in 0 until len) {
            val bases = sequences.map { it.sequence[col] }.toSet()
            if (bases.size > 1) S++
        }
        if (S == 0) return 0.0

        // Nucleotide diversity (pi)
        var pi = 0.0
        for (i in 0 until n) for (j in i + 1 until n) {
            var diff = 0
            for (col in 0 until len) if (sequences[i].sequence[col] != sequences[j].sequence[col]) diff++
            pi += diff.toDouble() / len
        }
        pi = 2.0 * pi / (n * (n - 1))

        // Watterson's theta
        val a1 = (1 until n).map { 1.0 / it }.sum()
        val theta = S / a1 / len

        // Variance (simplified)
        val a2 = (1 until n).map { 1.0 / (it * it) }.sum()
        val b1 = (n + 1).toDouble() / (3 * (n - 1))
        val b2 = 2.0 * (n * n + n + 3) / (9 * n * (n - 1))
        val c1 = b1 - 1.0 / a1
        val c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
        val e1 = c1 / a1; val e2 = c2 / (a1 * a1 + a2)
        val varD = e1 * S + e2 * S * (S - 1)

        return (pi - theta) / sqrt(varD)
    }

    // Fst between two populations
    fun fst(pop1: List<DNASequence>, pop2: List<DNASequence>): Double {
        val allPops = pop1 + pop2
        val n1 = pop1.size; val n2 = pop2.size
        val len = allPops[0].length

        var totalHt = 0.0; var totalHs = 0.0
        for (col in 0 until len) {
            val basesPop1 = pop1.map { it.sequence[col] }
            val basesPop2 = pop2.map { it.sequence[col] }
            val allBases = allPops.map { it.sequence[col] }

            val ht = 1.0 - allBases.groupingBy { it }.eachCount().values.map { (it.toDouble() / allBases.size).pow(2) }.sum()
            val hs1 = 1.0 - basesPop1.groupingBy { it }.eachCount().values.map { (it.toDouble() / n1).pow(2) }.sum()
            val hs2 = 1.0 - basesPop2.groupingBy { it }.eachCount().values.map { (it.toDouble() / n2).pow(2) }.sum()
            val hs = (hs1 + hs2) / 2.0

            totalHt += ht; totalHs += hs
        }
        val ht = totalHt / len; val hs = totalHs / len
        return if (ht == 0.0) 0.0 else (ht - hs) / ht
    }

    // dN/dS (Ka/Ks)
    fun dNdS(seq1: DNASequence, seq2: DNASequence, table: CodonTable = CodonTables.STANDARD): Pair<Double, Double> {
        val s1 = seq1.sequence; val s2 = seq2.sequence
        var synonymous = 0; var nonSynonymous = 0
        var totalS = 0.0; var totalN = 0.0

        var i = 0
        while (i + 3 <= minOf(s1.length, s2.length)) {
            val c1 = s1.substring(i, i + 3)
            val c2 = s2.substring(i, i + 3)
            val aa1 = try { table.translate(c1) } catch (e: Exception) { '?'; i += 3; continue }
            val aa2 = try { table.translate(c2) } catch (e: Exception) { '?'; i += 3; continue }

            val diffSites = c1.zip(c2).count { (a, b) -> a != b }
            totalS += 1.0 / 3; totalN += 2.0 / 3   // rough estimate

            if (diffSites > 0) {
                if (aa1 == aa2) synonymous++ else nonSynonymous++
            }
            i += 3
        }

        val ks = if (totalS == 0.0) 0.0 else {
            val p = synonymous / totalS
            if (p >= 0.75) Double.MAX_VALUE else -0.75 * ln(1 - 4 * p / 3)
        }
        val ka = if (totalN == 0.0) 0.0 else {
            val p = nonSynonymous / totalN
            if (p >= 0.75) Double.MAX_VALUE else -0.75 * ln(1 - 4 * p / 3)
        }
        return Pair(ka, ks)
    }
}

// ─── IUPAC Utilities ────────────────────────────────────────

object IUPAC {
    private val AMBIGUITY_MAP = mapOf(
        'R' to setOf('A', 'G'), 'Y' to setOf('C', 'T'), 'S' to setOf('G', 'C'),
        'W' to setOf('A', 'T'), 'K' to setOf('G', 'T'), 'M' to setOf('A', 'C'),
        'B' to setOf('C', 'G', 'T'), 'D' to setOf('A', 'G', 'T'),
        'H' to setOf('A', 'C', 'T'), 'V' to setOf('A', 'C', 'G'),
        'N' to setOf('A', 'T', 'G', 'C')
    )

    fun expandAmbiguous(sequence: String): List<String> {
        var results = listOf(StringBuilder())
        for (base in sequence.toUpperCase()) {
            val expansions = AMBIGUITY_MAP[base]?.toList() ?: listOf(base)
            results = results.flatMap { sb -> expansions.map { StringBuilder(sb).append(it) } }
        }
        return results.map { it.toString() }
    }

    fun ambiguityCode(bases: Set<Char>): Char {
        val upper = bases.map { it.toUpperCase() }.toSet()
        return when {
            upper == setOf('A', 'G') -> 'R'; upper == setOf('C', 'T') -> 'Y'
            upper == setOf('G', 'C') -> 'S'; upper == setOf('A', 'T') -> 'W'
            upper == setOf('G', 'T') -> 'K'; upper == setOf('A', 'C') -> 'M'
            upper == setOf('C', 'G', 'T') -> 'B'; upper == setOf('A', 'G', 'T') -> 'D'
            upper == setOf('A', 'C', 'T') -> 'H'; upper == setOf('A', 'C', 'G') -> 'V'
            upper.size == 4             -> 'N'
            upper.size == 1             -> upper.first()
            else -> 'N'
        }
    }
}
