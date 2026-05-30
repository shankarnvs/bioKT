package biokt

// ─── BLAST-like Hit & HSP ───────────────────────────────────

data class BlastHSP(
    val queryStart: Int,
    val queryEnd: Int,
    val subjectStart: Int,
    val subjectEnd: Int,
    val alignment: PairwiseAlignment,
    val eValue: Double,
    val bitScore: Double
) {
    val identities: Int get() = alignment.alignedA.zip(alignment.alignedB).count { (a, b) -> a == b && a != '-' }
    val gaps: Int get() = alignment.gaps
    val length: Int get() = alignment.alignedA.length
    override fun toString() =
        "HSP(query=$queryStart-$queryEnd, score=%.1f, e-value=%.2e, identity=%.1f%%)".format(
            bitScore, eValue, alignment.identity * 100)
}

data class BlastHit(
    val id: String,
    val description: String,
    val length: Int,
    val hsps: List<BlastHSP>
) {
    val bestHsp: BlastHSP? get() = hsps.minBy { it.eValue }
    override fun toString() = "Hit(id='$id', hsps=${hsps.size}, best_e=${bestHsp?.eValue?.let { "%.2e".format(it) }})"
}

data class BlastResult(
    val queryId: String,
    val queryLength: Int,
    val hits: List<BlastHit>,
    val database: String = ""
) {
    val numHits: Int get() = hits.size
    fun bestHits(n: Int = 10): List<BlastHit> =
        hits.sortedBy { it.bestHsp?.eValue ?: Double.MAX_VALUE }.take(n)
    override fun toString() = "BlastResult(query='$queryId', hits=$numHits)"
}

// ─── Sequence Database ──────────────────────────────────────

class SequenceDatabase(val name: String = "db") {
    private val records = mutableMapOf<String, SeqRecord>()

    fun add(record: SeqRecord) { records[record.id] = record }
    fun addAll(recs: Iterable<SeqRecord>) = recs.forEach { add(it) }

    fun get(id: String): SeqRecord? = records[id]
    fun ids(): Set<String> = records.keys
    val size: Int get() = records.size

    // Load from file
    fun loadFasta(filename: String) = addAll(SeqIO.parse(filename, "fasta"))

    // k-mer index for fast search
    private val kmerIndex = mutableMapOf<String, MutableList<Pair<String, Int>>>()
    private var kmerSize = 11

    fun buildIndex(k: Int = 11) {
        kmerSize = k
        kmerIndex.clear()
        for ((id, record) in records) {
            val seq = record.sequence.sequence
            for (i in 0..seq.length - k) {
                val kmer = seq.substring(i, i + k)
                kmerIndex.getOrPut(kmer) { mutableListOf() }.add(id to i)
            }
        }
    }

    // ── BLASTn (DNA) ─────────────────────────────────────────

    fun blastn(
        query: DNASequence,
        eValueThreshold: Double = 0.001,
        maxHits: Int = 50,
        wordSize: Int = 11,
        matchScore: Int = 2,
        mismatchPenalty: Int = -3,
        gapOpen: Int = 5,
        gapExtend: Int = 2
    ): BlastResult {
        if (kmerIndex.isEmpty()) buildIndex(wordSize)
        return search(query, eValueThreshold, maxHits, wordSize,
            matchScore.toDouble(), mismatchPenalty.toDouble(),
            gapOpen.toDouble(), gapExtend.toDouble())
    }

    // ── BLASTp (protein) ─────────────────────────────────────

    fun blastp(
        query: ProteinSequence,
        eValueThreshold: Double = 0.001,
        maxHits: Int = 50,
        wordSize: Int = 3
    ): BlastResult {
        if (kmerIndex.isEmpty()) buildIndex(wordSize)
        return search(query, eValueThreshold, maxHits, wordSize,
            matchScore = 2.0, mismatchScore = -1.0,
            gapOpen = 11.0, gapExtend = 1.0,
            substitutionMatrix = ScoringMatrices.BLOSUM62)
    }

    private fun search(
        query: BioSequence,
        eValueThreshold: Double,
        maxHits: Int,
        wordSize: Int,
        matchScore: Double,
        mismatchScore: Double,
        gapOpen: Double,
        gapExtend: Double,
        substitutionMatrix: Map<Pair<Char, Char>, Int>? = null
    ): BlastResult {
        val aligner = PairwiseAligner(
            mode = PairwiseAligner.Mode.LOCAL,
            matchScore = matchScore,
            mismatchScore = mismatchScore,
            gapOpen = -gapOpen,
            gapExtend = -gapExtend,
            substitutionMatrix = substitutionMatrix
        )

        // Seed with k-mers
        val candidates = mutableMapOf<String, MutableList<Int>>()
        val qSeq = query.sequence
        for (i in 0..qSeq.length - wordSize) {
            val kmer = qSeq.substring(i, i + wordSize)
            kmerIndex[kmer]?.forEach { (id, pos) ->
                candidates.getOrPut(id) { mutableListOf() }.add(pos)
            }
        }

        val hits = mutableListOf<BlastHit>()
        for ((id, _) in candidates.entries.sortedByDescending { it.value.size }.take(maxHits * 2)) {
            val subject = records[id] ?: continue
            val aln = aligner.align(qSeq, subject.sequence.sequence)
            if (aln.score <= 0) continue

            val eValue = calculateEValue(aln.score, qSeq.length, subject.length, size)
            if (eValue > eValueThreshold) continue

            val bitScore = lambdaK(aln.score)
            val hsp = BlastHSP(
                queryStart = 0, queryEnd = qSeq.length,
                subjectStart = 0, subjectEnd = subject.length,
                alignment = aln, eValue = eValue, bitScore = bitScore
            )
            hits.add(BlastHit(id, subject.description, subject.length, listOf(hsp)))
        }

        return BlastResult(
            queryId = query.id,
            queryLength = qSeq.length,
            hits = hits.sortedBy { it.bestHsp?.eValue ?: Double.MAX_VALUE }.take(maxHits),
            database = name
        )
    }

    // Simplified E-value calculation
    private fun calculateEValue(score: Double, qLen: Int, dbLen: Int, dbSeqs: Int): Double {
        val lambda = 0.318; val K = 0.134
        val mn = qLen.toLong() * (dbLen.toLong() * dbSeqs)
        return K * mn * Math.exp(-lambda * score)
    }

    private fun lambdaK(rawScore: Double): Double {
        val lambda = 0.318; val K = 0.134; val ln2 = Math.log(2.0)
        return (lambda * rawScore - Math.log(K)) / ln2
    }
}
