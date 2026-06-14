package biokt

// ============================================================
// BIOKT — A Bioinformatics Library for Kotlin (inspired by BioPython)
// ============================================================

// ─── Alphabets ───────────────────────────────────────────────

enum class Alphabet { DNA, RNA, PROTEIN, AMBIGUOUS_DNA, AMBIGUOUS_RNA }

object Alphabets {
    val DNA_BASES      = setOf('A', 'T', 'G', 'C')
    val RNA_BASES      = setOf('A', 'U', 'G', 'C')
    val AMINO_ACIDS    = setOf('A','C','D','E','F','G','H','I','K','L',
                                'M','N','P','Q','R','S','T','V','W','Y','*')
    val AMBIGUOUS_DNA  = DNA_BASES + setOf('N','R','Y','S','W','K','M','B','D','H','V')
    val AMBIGUOUS_RNA  = RNA_BASES + setOf('N','R','Y','S','W','K','M','B','D','H','V')

    fun validate(seq: String, alphabet: Alphabet): Boolean {
        val allowed = when (alphabet) {
            Alphabet.DNA            -> DNA_BASES
            Alphabet.RNA            -> RNA_BASES
            Alphabet.PROTEIN        -> AMINO_ACIDS
            Alphabet.AMBIGUOUS_DNA  -> AMBIGUOUS_DNA
            Alphabet.AMBIGUOUS_RNA  -> AMBIGUOUS_RNA
        }
        return seq.toUpperCase().all { it in allowed }
    }
}

// ─── Base Sequence ──────────────────────────────────────────

sealed class BioSequence(
    open val sequence: String,
    open val id: String = "",
    open val description: String = ""
) {
    val length: Int get() = sequence.length

    operator fun get(index: Int): Char = sequence[index]
    operator fun get(range: IntRange): String = sequence.substring(range)

    fun count(base: Char): Int = sequence.toUpperCase().count { it == base.toUpperCase() }

    fun countAll(): Map<Char, Int> =
        sequence.toUpperCase().groupingBy { it }.eachCount()

    override fun toString(): String = sequence

    // Motif search — returns all start indices of pattern
    fun findMotif(pattern: String): List<Int> {
        val results = mutableListOf<Int>()
        val seq = sequence.toUpperCase()
        val pat = pattern.toUpperCase()
        var idx = 0
        while (true) {
            val found = seq.indexOf(pat, idx)
            if (found == -1) break
            results += found
            idx = found + 1
        }
        return results
    }

    // GC content (DNA/RNA)
    fun gcContent(): Double {
        if (length == 0) return 0.0
        val gc = sequence.toUpperCase().count { it == 'G' || it == 'C' }
        return gc.toDouble() / length * 100.0
    }

    // GC skew in a window
    fun gcSkew(windowSize: Int = 100, stepSize: Int = 10): List<Pair<Int, Double>> {
        val results = mutableListOf<Pair<Int, Double>>()
        val seq = sequence.toUpperCase()
        var pos = 0
        while (pos + windowSize <= seq.length) {
            val window = seq.substring(pos, pos + windowSize)
            val g = window.count { it == 'G' }.toDouble()
            val c = window.count { it == 'C' }.toDouble()
            val skew = if (g + c == 0.0) 0.0 else (g - c) / (g + c)
            results += Pair(pos, skew)
            pos += stepSize
        }
        return results
    }

    abstract fun complement(): BioSequence
    abstract fun reverseComplement(): BioSequence
}

// ─── DNA Sequence ───────────────────────────────────────────

class DNASequence(
    sequence: String,
    id: String = "",
    description: String = ""
) : BioSequence(sequence.toUpperCase(), id, description) {

    init {
        require(Alphabets.validate(sequence, Alphabet.AMBIGUOUS_DNA)) {
            "Invalid DNA sequence: contains non-DNA characters in '$sequence'"
        }
    }

    companion object {
        private val COMPLEMENT_MAP = mapOf(
            'A' to 'T', 'T' to 'A', 'G' to 'C', 'C' to 'G',
            'N' to 'N', 'R' to 'Y', 'Y' to 'R', 'S' to 'S',
            'W' to 'W', 'K' to 'M', 'M' to 'K', 'B' to 'V',
            'D' to 'H', 'H' to 'D', 'V' to 'B'
        )
    }

    override fun complement(): DNASequence {
        val comp = sequence.map { COMPLEMENT_MAP[it] ?: it }.joinToString("")
        return DNASequence(comp, id, description)
    }

    override fun reverseComplement(): DNASequence =
        DNASequence(complement().sequence.reversed(), id, description)

    fun transcribe(): RNASequence =
        RNASequence(sequence.replace('T', 'U'), id, description)

    fun translate(
        table: CodonTable = CodonTables.STANDARD,
        startCodon: Boolean = false
    ): ProteinSequence {
        val seq = if (startCodon) {
            val start = sequence.indexOf("ATG")
            if (start == -1) throw IllegalArgumentException("No start codon (ATG) found")
            sequence.substring(start)
        } else sequence
        val protein = StringBuilder()
        var i = 0
        while (i + 3 <= seq.length) {
            val codon = seq.substring(i, i + 3)
            val aa = table.translate(codon)
            if (aa == '*') break
            protein.append(aa)
            i += 3
        }
        return ProteinSequence(protein.toString(), id, description)
    }

    // Six-frame translation
    fun sixFrameTranslation(table: CodonTable = CodonTables.STANDARD): Map<String, ProteinSequence> {
        val rc = reverseComplement()
        return mapOf(
            "+1" to DNASequence(sequence).translate(table),
            "+2" to DNASequence(sequence.substring(1)).translate(table),
            "+3" to DNASequence(sequence.substring(2)).translate(table),
            "-1" to DNASequence(rc.sequence).translate(table),
            "-2" to DNASequence(rc.sequence.substring(1)).translate(table),
            "-3" to DNASequence(rc.sequence.substring(2)).translate(table)
        )
    }

    // Find all ORFs (Open Reading Frames)
    fun findOrfs(minLength: Int = 100, table: CodonTable = CodonTables.STANDARD): List<Orf> {
        val orfs = mutableListOf<Orf>()
        val stopCodons = setOf("TAA", "TAG", "TGA")
        for (frame in 0..2) {
            var i = frame
            while (i + 3 <= sequence.length) {
                if (sequence.substring(i, i + 3) == "ATG") {
                    val startPos = i
                    var j = i + 3
                    while (j + 3 <= sequence.length) {
                        val codon = sequence.substring(j, j + 3)
                        if (codon in stopCodons) {
                            val orfSeq = sequence.substring(startPos, j + 3)
                            if (orfSeq.length >= minLength) {
                                orfs += Orf(
                                    start = startPos,
                                    end = j + 3,
                                    frame = frame + 1,
                                    strand = '+',
                                    sequence = DNASequence(orfSeq),
                                    protein = DNASequence(orfSeq.dropLast(3)).translate(table)
                                )
                            }
                            break
                        }
                        j += 3
                    }
                }
                i += 3
            }
        }
        return orfs.sortedBy { it.start }
    }

    fun meltingTemperature(): Double {
        val a = count('A'); val t = count('T')
        val g = count('G'); val c = count('C')
        return if (length < 14)
            (a + t) * 2.0 + (g + c) * 4.0
        else
            64.9 + 41.0 * (g + c - 16.4) / length
    }

    fun subSequence(start: Int, end: Int): DNASequence =
        DNASequence(sequence.substring(start, end), id, description)
}

// ─── RNA Sequence ───────────────────────────────────────────

class RNASequence(
    sequence: String,
    id: String = "",
    description: String = ""
) : BioSequence(sequence.toUpperCase(), id, description) {

    init {
        require(Alphabets.validate(sequence, Alphabet.AMBIGUOUS_RNA)) {
            "Invalid RNA sequence: contains non-RNA characters in '$sequence'"
        }
    }

    companion object {
        private val COMPLEMENT_MAP = mapOf(
            'A' to 'U', 'U' to 'A', 'G' to 'C', 'C' to 'G'
        )
    }

    override fun complement(): RNASequence {
        val comp = sequence.map { COMPLEMENT_MAP[it] ?: it }.joinToString("")
        return RNASequence(comp, id, description)
    }

    override fun reverseComplement(): RNASequence =
        RNASequence(complement().sequence.reversed(), id, description)

    fun backTranscribe(): DNASequence =
        DNASequence(sequence.replace('U', 'T'), id, description)

    fun translate(table: CodonTable = CodonTables.STANDARD): ProteinSequence {
        val dna = backTranscribe()
        return dna.translate(table)
    }
}

// ─── Protein Sequence ───────────────────────────────────────

class ProteinSequence(
    sequence: String,
    id: String = "",
    description: String = ""
) : BioSequence(sequence.toUpperCase(), id, description) {

    override fun complement(): ProteinSequence =
        throw UnsupportedOperationException("Complement not defined for protein sequences")

    override fun reverseComplement(): ProteinSequence =
        throw UnsupportedOperationException("Reverse complement not defined for protein sequences")

    fun molecularWeight(): Double {
        val weights = mapOf(
            'A' to 89.09, 'R' to 174.20, 'N' to 132.12, 'D' to 133.10,
            'C' to 121.16, 'E' to 147.13, 'Q' to 146.15, 'G' to 75.03,
            'H' to 155.16, 'I' to 131.17, 'L' to 131.17, 'K' to 146.19,
            'M' to 149.21, 'F' to 165.19, 'P' to 115.13, 'S' to 105.09,
            'T' to 119.12, 'W' to 204.23, 'Y' to 181.19, 'V' to 117.15
        )
        val water = 18.02
        return sequence.map { weights[it] ?: 0.0 }.sum() - (length - 1) * water
    }

    fun isoelectricPoint(): Double {
        // Henderson-Hasselbalch approximation
        val pKa = mapOf(
            'D' to 3.9, 'E' to 4.1, 'H' to 6.0, 'C' to 8.3,
            'Y' to 10.1, 'K' to 10.5, 'R' to 12.5
        )
        val nTerm = 8.0; val cTerm = 3.1
        var pH = 7.0
        val step = 1.0
        repeat(1000) {
            // Binary search for charge = 0
        }
        // Simplified: return 7.0 as placeholder; production code would iterate
        var lo = 0.0; var hi = 14.0
        repeat(100) {
            val mid = (lo + hi) / 2.0
            val charge = chargeAtPH(mid, pKa, nTerm, cTerm)
            if (charge > 0) lo = mid else hi = mid
        }
        return (lo + hi) / 2.0
    }

    private fun chargeAtPH(pH: Double, pKa: Map<Char, Double>, nTerm: Double, cTerm: Double): Double {
        var charge = 1.0 / (1.0 + Math.pow(10.0, pH - nTerm))  // N-terminus
        charge -= 1.0 / (1.0 + Math.pow(10.0, cTerm - pH))     // C-terminus
        for (aa in sequence) {
            val pk = pKa[aa] ?: continue
            charge += when (aa) {
                'D', 'E', 'C', 'Y' -> -1.0 / (1.0 + Math.pow(10.0, pk - pH))
                'H', 'K', 'R'      ->  1.0 / (1.0 + Math.pow(10.0, pH - pk))
                else -> 0.0
            }
        }
        return charge
    }

    fun aminoAcidComposition(): Map<Char, Double> {
        val counts = countAll().filter { it.key != '*' }
        val total = counts.values.sum().toDouble()
        return counts.mapValues { (_, v) -> v / total * 100.0 }
    }

    fun aromaticity(): Double {
        val aromatic = sequence.count { it in setOf('F', 'H', 'W', 'Y') }
        return aromatic.toDouble() / length
    }

    fun instabilityIndex(): Double {
        // DIWV method (Guruprasad et al. 1990)
        val diwv = mapOf(
            "WW" to 1.0, "WC" to 1.0, "WE" to 1.0, "WT" to -0.7,
            "AA" to 1.0, "AC" to 44.94, "AE" to 0.0, "AD" to -7.49,
            "AF" to -14.03, "AG" to 1.0, "AH" to -7.49, "AI" to -7.49,
            "AK" to -7.49, "AL" to -7.49, "AM" to 1.0, "AN" to 23.44,
            "AP" to 20.26, "AQ" to 0.0, "AR" to 0.0, "AS" to 1.0
        )
        var score = 0.0
        for (i in 0 until length - 1) {
            val dipeptide = "${sequence[i]}${sequence[i+1]}"
            score += diwv[dipeptide] ?: 1.0
        }
        return 10.0 / length * score
    }

    fun subSequence(start: Int, end: Int): ProteinSequence =
        ProteinSequence(sequence.substring(start, end), id, description)
}

// ─── ORF Data Class ─────────────────────────────────────────

data class Orf(
    val start: Int,
    val end: Int,
    val frame: Int,
    val strand: Char,
    val sequence: DNASequence,
    val protein: ProteinSequence
) {
    val length: Int get() = end - start
    override fun toString() =
        "ORF[frame=$strand$frame, pos=$start..$end, len=$length, protein=${protein.sequence.take(20)}...]"
}
