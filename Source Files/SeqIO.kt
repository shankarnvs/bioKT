package biokt

import java.io.File
import java.io.BufferedReader
import java.io.PrintWriter

// ─── SeqRecord (analogous to BioPython's SeqRecord) ─────────

data class SeqRecord(
    val sequence: BioSequence,
    val id: String = "",
    val name: String = "",
    val description: String = "",
    val features: List<SeqFeature> = emptyList(),
    val annotations: Map<String, Any> = emptyMap(),
    val dbxrefs: List<String> = emptyList()
) {
    val length: Int get() = sequence.length
    override fun toString() = "SeqRecord(id='$id', length=$length)"
}

data class SeqFeature(
    val type: String,
    val location: FeatureLocation,
    val strand: Int = 1,          // 1 = forward, -1 = reverse
    val qualifiers: Map<String, List<String>> = emptyMap()
) {
    override fun toString() = "SeqFeature(type='$type', location=$location)"
}

data class FeatureLocation(
    val start: Int,
    val end: Int,
    val strand: Int = 1
) {
    val length: Int get() = end - start
    override fun toString() = "${start}..${end}"
}

// ─── SeqIO (analogous to BioPython's SeqIO) ─────────────────

object SeqIO {

    // ── Read from file ───────────────────────────────────────

    fun read(filename: String, format: String): SeqRecord {
        val records = parse(filename, format)
        if (records.isEmpty()) throw IllegalArgumentException("No records found in $filename")
        if (records.size > 1) throw IllegalArgumentException("More than one record in $filename — use parse()")
        return records.first()
    }

    fun parse(filename: String, format: String): List<SeqRecord> =
        parse(File(filename).bufferedReader(), format)

    fun parse(reader: BufferedReader, format: String): List<SeqRecord> {
        return when (format.toLowerCase()) {
            "fasta", "fa"  -> parseFasta(reader)
            "fastq"        -> parseFastq(reader)
            "genbank", "gb"-> parseGenbank(reader)
            "stockholm"    -> parseStockholm(reader)
            else           -> throw IllegalArgumentException("Unknown format: $format")
        }
    }

    fun parseString(text: String, format: String): List<SeqRecord> =
        parse(text.trimIndent().reader().buffered(), format)

    // ── Write to file ────────────────────────────────────────

    fun write(records: List<SeqRecord>, filename: String, format: String) {
        val writer = PrintWriter(File(filename))
        writer.use { write(records, it, format) }
    }

    fun write(records: List<SeqRecord>, writer: PrintWriter, format: String) {
        when (format.toLowerCase()) {
            "fasta"        -> records.forEach { writeFastaRecord(it, writer) }
            "fastq"        -> records.forEach { writeFastqRecord(it, writer) }
            else           -> throw IllegalArgumentException("Unknown format: $format")
        }
    }

    fun toString(records: List<SeqRecord>, format: String): String {
        val sb = StringBuilder()
        val pw = PrintWriter(object : java.io.Writer() {
            override fun write(buf: CharArray, off: Int, len: Int) = sb.append(String(buf, off, len)).let {}
            override fun flush() {}
            override fun close() {}
        })
        write(records, pw, format)
        pw.flush()
        return sb.toString()
    }

    // ── FASTA Parser ─────────────────────────────────────────

    private fun parseFasta(reader: BufferedReader): List<SeqRecord> {
        val records = mutableListOf<SeqRecord>()
        var currentId = ""
        var currentDesc = ""
        val seqBuilder = StringBuilder()

        fun flush() {
            if (currentId.isNotEmpty()) {
                val seq = seqBuilder.toString()
                val bioSeq = inferSequence(seq, currentId, currentDesc)
                records.add(SeqRecord(bioSeq, currentId, currentId, currentDesc))
                seqBuilder.clear()
            }
        }

        reader.forEachLine { line ->
            val trimmed = line.trim()
            when {
                trimmed.startsWith(">") -> {
                    flush()
                    val header = trimmed.drop(1)
                    val spaceIdx = header.indexOf(' ')
                    if (spaceIdx == -1) { currentId = header; currentDesc = "" }
                    else { currentId = header.substring(0, spaceIdx); currentDesc = header.substring(spaceIdx + 1) }
                }
                trimmed.isNotEmpty() && !trimmed.startsWith(";") -> seqBuilder.append(trimmed)
            }
        }
        flush()
        return records
    }

    private fun writeFastaRecord(record: SeqRecord, writer: PrintWriter) {
        val header = if (record.description.isNotEmpty())
            "${record.id} ${record.description}" else record.id
        writer.println(">$header")
        record.sequence.sequence.chunked(80).forEach { writer.println(it) }
    }

    // ── FASTQ Parser ─────────────────────────────────────────

    private fun parseFastq(reader: BufferedReader): List<SeqRecord> {
        val records = mutableListOf<SeqRecord>()
        val lines = reader.readLines()
        var i = 0
        while (i + 3 < lines.size) {
            val header = lines[i].drop(1)
            val seq = lines[i + 1]
            // lines[i+2] is "+"
            val quality = lines[i + 3]
            val (id, desc) = if (header.contains(' ')) header.split(' ', limit = 2).let { it[0] to it[1] }
                             else header to ""
            val bioSeq = inferSequence(seq, id, desc)
            records.add(SeqRecord(
                sequence = bioSeq, id = id, name = id, description = desc,
                annotations = mapOf("phred_quality" to quality.map { it.toInt() - 33 })
            ))
            i += 4
        }
        return records
    }

    private fun writeFastqRecord(record: SeqRecord, writer: PrintWriter) {
        val quality = (record.annotations["phred_quality"] as? List<*>)
            ?.joinToString("") { ((it as Int) + 33).toChar().toString() }
            ?: "I".repeat(record.length)
        val header = if (record.description.isNotEmpty())
            "${record.id} ${record.description}" else record.id
        writer.println("@$header")
        writer.println(record.sequence.sequence)
        writer.println("+")
        writer.println(quality)
    }

    // ── GenBank Parser (simplified) ──────────────────────────

    private fun parseGenbank(reader: BufferedReader): List<SeqRecord> {
        val records = mutableListOf<SeqRecord>()
        val lines = reader.readLines()
        var i = 0
        var id = ""; var desc = ""; val features = mutableListOf<SeqFeature>()
        val annotations = mutableMapOf<String, Any>()
        val seqBuilder = StringBuilder()

        while (i < lines.size) {
            val line = lines[i]
            when {
                line.startsWith("LOCUS") -> {
                    val parts = line.trim().split("\\s+".toRegex())
                    if (parts.size > 1) id = parts[1]
                }
                line.startsWith("DEFINITION") -> {
                    desc = line.removePrefix("DEFINITION").trim()
                }
                line.startsWith("ACCESSION") -> {
                    annotations["accession"] = line.removePrefix("ACCESSION").trim()
                }
                line.startsWith("ORGANISM") -> {
                    annotations["organism"] = line.removePrefix("ORGANISM").trim()
                }
                line.startsWith("     gene") || line.startsWith("     CDS") || line.startsWith("     exon") -> {
                    val type = line.trim().split("\\s+".toRegex())[0]
                    val locStr = line.trim().split("\\s+".toRegex()).getOrElse(1) { "1..1" }
                    val loc = parseLocation(locStr)
                    val qualMap = mutableMapOf<String, List<String>>()
                    var k = i + 1
                    while (k < lines.size && lines[k].startsWith("                     /")) {
                        val qual = lines[k].trim().removePrefix("/")
                        val eqIdx = qual.indexOf('=')
                        if (eqIdx >= 0) {
                            val qKey = qual.substring(0, eqIdx)
                            val qVal = qual.substring(eqIdx + 1).trim('"')
                            qualMap[qKey] = qualMap.getOrDefault(qKey, emptyList()) + qVal
                        }
                        k++
                    }
                    features.add(SeqFeature(type, loc, qualifiers = qualMap))
                    i = k - 1
                }
                line.startsWith("ORIGIN") -> {
                    i++
                    while (i < lines.size && !lines[i].startsWith("//")) {
                        val seq = lines[i].replace("\\s+".toRegex(), "").filter { it.isLetter() }
                        seqBuilder.append(seq)
                        i++
                    }
                }
                line.startsWith("//") -> {
                    val seq = seqBuilder.toString()
                    if (seq.isNotEmpty()) {
                        records.add(SeqRecord(
                            sequence = inferSequence(seq, id, desc),
                            id = id, name = id, description = desc,
                            features = features.toList(),
                            annotations = annotations.toMap()
                        ))
                    }
                    id = ""; desc = ""; features.clear(); annotations.clear(); seqBuilder.clear()
                }
            }
            i++
        }
        return records
    }

    private fun parseLocation(loc: String): FeatureLocation {
        val clean = loc.replace("complement(", "").replace("join(", "").trimEnd(')')
        return if (clean.contains("..")) {
            val parts = clean.split("..").map { it.replace("<", "").replace(">", "").toIntOrNull() ?: 0 }
            FeatureLocation(parts[0] - 1, parts[1], if (loc.contains("complement")) -1 else 1)
        } else {
            val pos = clean.toIntOrNull() ?: 0
            FeatureLocation(pos - 1, pos)
        }
    }

    // ── Stockholm Parser (for MSA files) ─────────────────────

    private fun parseStockholm(reader: BufferedReader): List<SeqRecord> {
        val seqMap = mutableMapOf<String, StringBuilder>()
        val idOrder = mutableListOf<String>()
        reader.forEachLine { line ->
            when {
                line.startsWith("//") || line.startsWith("#") || line.isBlank() -> {}
                else -> {
                    val parts = line.split("\\s+".toRegex())
                    if (parts.size >= 2) {
                        val seqId = parts[0]
                        val seqFrag = parts[1]
                        if (seqId !in seqMap) { seqMap[seqId] = StringBuilder(); idOrder.add(seqId) }
                        seqMap[seqId]!!.append(seqFrag)
                    }
                }
            }
        }
        return idOrder.map { id ->
            val seq = seqMap[id]!!.toString()
            SeqRecord(inferSequence(seq.replace("-", ""), id, ""), id, id, "")
        }
    }

    // ── Utility ──────────────────────────────────────────────

    private fun inferSequence(seq: String, id: String, desc: String): BioSequence {
        val upper = seq.toUpperCase()
        return when {
            upper.contains('U') -> RNASequence(upper, id, desc)
            upper.all { it in Alphabets.AMINO_ACIDS } &&
            upper.any { it in setOf('E','F','I','J','L','O','P','Q','X','Z') } ->
                ProteinSequence(upper, id, desc)
            else -> try { DNASequence(upper, id, desc) } catch (e: Exception) {
                try { RNASequence(upper, id, desc) } catch (e2: Exception) {
                    ProteinSequence(upper, id, desc)
                }
            }
        }
    }

    fun convert(inputFile: String, inputFormat: String, outputFile: String, outputFormat: String): Int {
        val records = parse(inputFile, inputFormat)
        write(records, outputFile, outputFormat)
        return records.size
    }
}
