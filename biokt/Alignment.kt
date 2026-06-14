package biokt

// ─── Alignment Results ──────────────────────────────────────

data class PairwiseAlignment(
    val seqA: String,
    val seqB: String,
    val alignedA: String,
    val alignedB: String,
    val score: Double,
    val identity: Double,      // fraction identical positions
    val similarity: Double,    // fraction similar positions
    val gaps: Int
) {
    fun prettyPrint(lineWidth: Int = 60): String {
        val sb = StringBuilder()
        sb.append("Score: %.2f  Identity: %.1f%%  Similarity: %.1f%%  Gaps: $gaps".format(score, identity * 100, similarity * 100)).append("\n")
        var i = 0
        while (i < alignedA.length) {
            val end = minOf(i + lineWidth, alignedA.length)
            val a = alignedA.substring(i, end)
            val b = alignedB.substring(i, end)
            val mid = a.zip(b).joinToString("") { (x, y) ->
                when {
                    x == y && x != '-'    -> "|"
                    x == '-' || y == '-'  -> " "
                    isSimilar(x, y)       -> "+"
                    else                  -> " "
                }
            }
            sb.append("Query  ${i+1}\t$a  ${i + a.replace("-","").length}").append("\n")
            sb.append("       \t$mid").append("\n")
            sb.append("Sbjct  ${i+1}\t$b  ${i + b.replace("-","").length}").append("\n")
            sb.append("\n")
            i += lineWidth
        }
        return sb.toString()
    }

    private fun isSimilar(a: Char, b: Char): Boolean {
        val groups = listOf(
            setOf('I','L','V','M'), setOf('F','Y','W'), setOf('K','R','H'),
            setOf('D','E'), setOf('S','T'), setOf('N','Q')
        )
        return groups.any { a in it && b in it }
    }
}

data class MultipleAlignment(
    val sequences: List<String>,
    val ids: List<String>,
    val alignedSequences: List<String>
) {
    val length: Int get() = alignedSequences.firstOrNull()?.length ?: 0
    val numSeqs: Int get() = alignedSequences.size

    fun consensusSequence(threshold: Double = 0.5): String {
        return (0 until length).map { col ->
            val chars = alignedSequences.map { it[col] }.filter { it != '-' }
            if (chars.isEmpty()) '-'
            else {
                val freq = chars.groupingBy { it }.eachCount()
                val best = freq.maxBy { it.value }!!
                if (best.value.toDouble() / chars.size >= threshold) best.key else 'X'
            }
        }.joinToString("")
    }

    fun conservationScore(): List<Double> {
        return (0 until length).map { col ->
            val chars = alignedSequences.map { it[col] }.filter { it != '-' }
            if (chars.isEmpty()) 0.0
            else {
                val freq = chars.groupingBy { it }.eachCount()
                val max = freq.values.max()!!.toDouble()
                max / chars.size
            }
        }
    }

    fun percentIdentity(): Double {
        if (numSeqs < 2) return 100.0
        var totalPairs = 0; var identicalPairs = 0
        for (i in 0 until numSeqs) {
            for (j in i + 1 until numSeqs) {
                for (col in 0 until length) {
                    val a = alignedSequences[i][col]; val b = alignedSequences[j][col]
                    if (a != '-' && b != '-') {
                        totalPairs++
                        if (a == b) identicalPairs++
                    }
                }
            }
        }
        return if (totalPairs == 0) 0.0 else identicalPairs.toDouble() / totalPairs * 100.0
    }

    fun prettyPrint(lineWidth: Int = 60): String {
        val maxIdLen = ids.map { it.length }.max()!!
        val sb = StringBuilder()
        var i = 0
        while (i < length) {
            val end = minOf(i + lineWidth, length)
            for (k in 0 until numSeqs) {
                sb.append(ids[k].padEnd(maxIdLen + 2) + alignedSequences[k].substring(i, end)).append("\n")
            }
            sb.append("\n")
            i += lineWidth
        }
        return sb.toString()
    }
}

// ─── Scoring Matrices ───────────────────────────────────────

object ScoringMatrices {

    // BLOSUM62
    val BLOSUM62: Map<Pair<Char, Char>, Int> by lazy {
        val aas = "ARNDCQEGHILKMFPSTWYV"
        val raw = arrayOf(
            intArrayOf( 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0),
            intArrayOf(-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3),
            intArrayOf(-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3),
            intArrayOf(-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3),
            intArrayOf( 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1),
            intArrayOf(-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2),
            intArrayOf(-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2),
            intArrayOf( 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3),
            intArrayOf(-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3),
            intArrayOf(-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3),
            intArrayOf(-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1),
            intArrayOf(-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2),
            intArrayOf(-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1),
            intArrayOf(-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1),
            intArrayOf(-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2),
            intArrayOf( 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2),
            intArrayOf( 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0),
            intArrayOf(-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3),
            intArrayOf(-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1),
            intArrayOf( 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4)
        )
        val m = mutableMapOf<Pair<Char, Char>, Int>()
        for (i in aas.indices) for (j in aas.indices) m[aas[i] to aas[j]] = raw[i][j]
        m
    }

    // PAM250
    val PAM250: Map<Pair<Char, Char>, Int> by lazy {
        // Simplified — diagonal only for brevity; full matrix in production
        val aas = "ARNDCQEGHILKMFPSTWYV"
        val m = mutableMapOf<Pair<Char, Char>, Int>()
        for (a in aas) for (b in aas) m[a to b] = if (a == b) 2 else -1
        m
    }

    // Simple DNA match/mismatch
    fun dnaScoringMatrix(match: Int = 2, mismatch: Int = -1): Map<Pair<Char, Char>, Int> {
        val bases = "ATGCN"
        val m = mutableMapOf<Pair<Char, Char>, Int>()
        for (a in bases) for (b in bases)
            m[a to b] = if (a == b) match else mismatch
        return m
    }
}

// ─── Pairwise Aligner ───────────────────────────────────────

class PairwiseAligner(
    val mode: Mode = Mode.GLOBAL,
    val matchScore: Double = 2.0,
    val mismatchScore: Double = -1.0,
    val gapOpen: Double = -2.0,
    val gapExtend: Double = -0.5,
    val substitutionMatrix: Map<Pair<Char, Char>, Int>? = null
) {
    enum class Mode { GLOBAL, LOCAL, SEMIGLOBAL }

    fun score(a: Char, b: Char): Double {
        return substitutionMatrix?.get(a to b)?.toDouble()
            ?: if (a == b) matchScore else mismatchScore
    }

    fun align(seqA: String, seqB: String): PairwiseAlignment {
        return when (mode) {
            Mode.GLOBAL     -> needlemanWunsch(seqA, seqB)
            Mode.LOCAL      -> smithWaterman(seqA, seqB)
            Mode.SEMIGLOBAL -> semiglobalAlign(seqA, seqB)
        }
    }

    // ── Needleman-Wunsch (global) ────────────────────────────
    private fun needlemanWunsch(seqA: String, seqB: String): PairwiseAlignment {
        val m = seqA.length; val n = seqB.length
        val dp = Array(m + 1) { DoubleArray(n + 1) }

        // Init
        for (i in 0..m) dp[i][0] = gapOpen + (i - 1) * gapExtend
        for (j in 0..n) dp[0][j] = gapOpen + (j - 1) * gapExtend
        dp[0][0] = 0.0

        for (i in 1..m) {
            for (j in 1..n) {
                val diag = dp[i-1][j-1] + score(seqA[i-1], seqB[j-1])
                val up   = dp[i-1][j]   + (if (i > 1 && dp[i-1][j] == dp[i-2][j] + gapExtend) gapExtend else gapOpen)
                val left = dp[i][j-1]   + (if (j > 1 && dp[i][j-1] == dp[i][j-2] + gapExtend) gapExtend else gapOpen)
                dp[i][j] = maxOf(diag, up, left)
            }
        }

        return traceback(seqA, seqB, dp, m, n, dp[m][n])
    }

    // ── Smith-Waterman (local) ───────────────────────────────
    private fun smithWaterman(seqA: String, seqB: String): PairwiseAlignment {
        val m = seqA.length; val n = seqB.length
        val dp = Array(m + 1) { DoubleArray(n + 1) { 0.0 } }
        var maxScore = 0.0; var maxI = 0; var maxJ = 0

        for (i in 1..m) {
            for (j in 1..n) {
                val diag = dp[i-1][j-1] + score(seqA[i-1], seqB[j-1])
                val up   = dp[i-1][j]   + gapOpen
                val left = dp[i][j-1]   + gapOpen
                dp[i][j] = listOf(0.0, diag, up, left).max()!!
                if (dp[i][j] > maxScore) { maxScore = dp[i][j]; maxI = i; maxJ = j }
            }
        }
        return traceback(seqA, seqB, dp, maxI, maxJ, maxScore, local = true)
    }

    // ── Semiglobal ───────────────────────────────────────────
    private fun semiglobalAlign(seqA: String, seqB: String): PairwiseAlignment {
        val m = seqA.length; val n = seqB.length
        val dp = Array(m + 1) { DoubleArray(n + 1) }
        for (i in 0..m) dp[i][0] = 0.0
        for (j in 0..n) dp[0][j] = 0.0

        for (i in 1..m) for (j in 1..n) {
            dp[i][j] = maxOf(
                dp[i-1][j-1] + score(seqA[i-1], seqB[j-1]),
                dp[i-1][j]   + gapOpen,
                dp[i][j-1]   + gapOpen
            )
        }
        val bestJ = (0..n).maxBy { dp[m][it] }!!
        return traceback(seqA, seqB, dp, m, bestJ, dp[m][bestJ])
    }

    private fun traceback(
        seqA: String, seqB: String,
        dp: Array<DoubleArray>, endI: Int, endJ: Int,
        totalScore: Double, local: Boolean = false
    ): PairwiseAlignment {
        val alignA = StringBuilder(); val alignB = StringBuilder()
        var i = endI; var j = endJ

        while (i > 0 && j > 0) {
            if (local && dp[i][j] == 0.0) break
            val s = score(seqA[i-1], seqB[j-1])
            when {
                dp[i][j] == dp[i-1][j-1] + s -> { alignA.append(seqA[i-1]); alignB.append(seqB[j-1]); i--; j-- }
                dp[i][j] == dp[i-1][j] + gapOpen -> { alignA.append(seqA[i-1]); alignB.append('-'); i-- }
                else -> { alignA.append('-'); alignB.append(seqB[j-1]); j-- }
            }
        }
        while (i > 0 && !local) { alignA.append(seqA[i-1]); alignB.append('-'); i-- }
        while (j > 0 && !local) { alignA.append('-'); alignB.append(seqB[j-1]); j-- }

        val a = alignA.reverse().toString()
        val b = alignB.reverse().toString()

        val identical = a.zip(b).count { (x, y) -> x == y && x != '-' }
        val gapCount = a.count { it == '-' } + b.count { it == '-' }
        val aligned = a.zip(b).count { (x, y) -> x != '-' && y != '-' }

        return PairwiseAlignment(
            seqA = seqA, seqB = seqB,
            alignedA = a, alignedB = b,
            score = totalScore,
            identity = if (aligned == 0) 0.0 else identical.toDouble() / aligned,
            similarity = if (aligned == 0) 0.0 else identical.toDouble() / aligned,
            gaps = gapCount
        )
    }
}

// ─── Multiple Sequence Alignment (progressive) ──────────────

class MultipleSequenceAligner {
    private val pairwiseAligner = PairwiseAligner(mode = PairwiseAligner.Mode.GLOBAL)

    fun align(sequences: List<BioSequence>): MultipleAlignment {
        if (sequences.isEmpty()) return MultipleAlignment(emptyList(), emptyList(), emptyList())
        if (sequences.size == 1) return MultipleAlignment(
            listOf(sequences[0].sequence), listOf(sequences[0].id), listOf(sequences[0].sequence)
        )

        // Build pairwise distance matrix
        val n = sequences.size
        val dist = Array(n) { DoubleArray(n) }
        for (i in 0 until n) for (j in i + 1 until n) {
            val aln = pairwiseAligner.align(sequences[i].sequence, sequences[j].sequence)
            dist[i][j] = 1.0 - aln.identity
            dist[j][i] = dist[i][j]
        }

        // UPGMA clustering for guide tree
        val order = upgmaOrder(dist, n)

        // Progressive alignment
        val aligned = sequences.map { it.sequence }.toMutableList()
        val alignedIds = sequences.map { it.id }.toMutableList()

        var currentAligned = mutableListOf(aligned[order[0]])
        for (k in 1 until order.size) {
            val next = aligned[order[k]]
            val profile = profileConsensus(currentAligned)
            val aln = pairwiseAligner.align(profile, next)
            currentAligned = propagateGaps(currentAligned, aln.alignedA)
            currentAligned.add(aln.alignedB)
        }

        return MultipleAlignment(
            sequences = sequences.map { it.sequence },
            ids = sequences.map { it.id },
            alignedSequences = currentAligned
        )
    }

    private fun upgmaOrder(dist: Array<DoubleArray>, n: Int): List<Int> {
        // Simple nearest-neighbor ordering
        val remaining = (0 until n).toMutableList()
        val order = mutableListOf<Int>()
        order.add(remaining.removeAt(0))
        while (remaining.isNotEmpty()) {
            val last = order.last()
            val nearest = remaining.minBy { dist[last][it] }!!
            order.add(nearest)
            remaining.remove(nearest)
        }
        return order
    }

    private fun profileConsensus(seqs: List<String>): String {
        if (seqs.isEmpty()) return ""
        val len = seqs[0].length
        return (0 until len).map { col ->
            val chars = seqs.map { it[col] }.filter { it != '-' }
            if (chars.isEmpty()) '-'
            else chars.groupingBy { it }.eachCount().maxBy { it.value }!!.key
        }.joinToString("")
    }

    private fun propagateGaps(seqs: List<String>, gappedProfile: String): MutableList<String> {
        if (seqs.isEmpty()) return mutableListOf()
        val result = seqs.map { StringBuilder(it) }.toMutableList()
        val orig = profileConsensus(seqs)
        var origIdx = 0
        gappedProfile.forEachIndexed { pos, c ->
            if (c == '-' && origIdx <= orig.length) {
                result.forEach { it.insert(origIdx, '-') }
                origIdx++
            } else origIdx++
        }
        return result.map { it.toString() }.toMutableList()
    }
}
