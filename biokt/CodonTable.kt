package biokt

class CodonTable(
    val name: String,
    val id: Int,
    private val table: Map<String, Char>
) {
    fun translate(codon: String): Char =
        table[codon.toUpperCase()] ?: throw IllegalArgumentException("Unknown codon: $codon")

    fun isStopCodon(codon: String): Boolean = translate(codon) == '*'
    fun isStartCodon(codon: String): Boolean = codon.toUpperCase() == "ATG"
    fun stopCodons(): List<String> = table.entries.filter { it.value == '*' }.map { it.key }
    fun startCodons(): List<String> = listOf("ATG")
}

private fun buildCodonMap(vararg pairs: Pair<String, Char>): Map<String, Char> = mapOf(*pairs)

object CodonTables {

    val STANDARD: CodonTable by lazy {
        val m = mutableMapOf<String, Char>()
        m["TTT"] = 'F'; m["TTC"] = 'F'
        m["TTA"] = 'L'; m["TTG"] = 'L'
        m["CTT"] = 'L'; m["CTC"] = 'L'; m["CTA"] = 'L'; m["CTG"] = 'L'
        m["ATT"] = 'I'; m["ATC"] = 'I'; m["ATA"] = 'I'
        m["ATG"] = 'M'
        m["GTT"] = 'V'; m["GTC"] = 'V'; m["GTA"] = 'V'; m["GTG"] = 'V'
        m["TCT"] = 'S'; m["TCC"] = 'S'; m["TCA"] = 'S'; m["TCG"] = 'S'
        m["AGT"] = 'S'; m["AGC"] = 'S'
        m["CCT"] = 'P'; m["CCC"] = 'P'; m["CCA"] = 'P'; m["CCG"] = 'P'
        m["ACT"] = 'T'; m["ACC"] = 'T'; m["ACA"] = 'T'; m["ACG"] = 'T'
        m["GCT"] = 'A'; m["GCC"] = 'A'; m["GCA"] = 'A'; m["GCG"] = 'A'
        m["TAT"] = 'Y'; m["TAC"] = 'Y'
        m["TAA"] = '*'; m["TAG"] = '*'; m["TGA"] = '*'
        m["CAT"] = 'H'; m["CAC"] = 'H'
        m["CAA"] = 'Q'; m["CAG"] = 'Q'
        m["AAT"] = 'N'; m["AAC"] = 'N'
        m["AAA"] = 'K'; m["AAG"] = 'K'
        m["GAT"] = 'D'; m["GAC"] = 'D'
        m["GAA"] = 'E'; m["GAG"] = 'E'
        m["TGT"] = 'C'; m["TGC"] = 'C'
        m["TGG"] = 'W'
        m["CGT"] = 'R'; m["CGC"] = 'R'; m["CGA"] = 'R'; m["CGG"] = 'R'
        m["AGA"] = 'R'; m["AGG"] = 'R'
        m["GGT"] = 'G'; m["GGC"] = 'G'; m["GGA"] = 'G'; m["GGG"] = 'G'
        CodonTable("Standard", 1, m)
    }

    val VERTEBRATE_MITOCHONDRIAL: CodonTable by lazy {
        val m = STANDARD.let { st ->
            val mm = mutableMapOf<String, Char>()
            listOf("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG",
                   "ATT","ATC","ATG","GTT","GTC","GTA","GTG",
                   "TCT","TCC","TCA","TCG","AGT","AGC","CCT","CCC","CCA","CCG",
                   "ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG",
                   "TAT","TAC","TAA","TAG","CAT","CAC","CAA","CAG","AAT","AAC",
                   "AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGG",
                   "CGT","CGC","CGA","CGG","GGT","GGC","GGA","GGG"
            ).forEach { codon -> mm[codon] = st.translate(codon) }
            mm["AGA"] = '*'; mm["AGG"] = '*'
            mm["ATA"] = 'M'; mm["TGA"] = 'W'
            mm
        }
        CodonTable("Vertebrate Mitochondrial", 2, m)
    }

    val YEAST_MITOCHONDRIAL: CodonTable by lazy {
        val m = STANDARD.let { st ->
            val mm = mutableMapOf<String, Char>()
            listOf("TTT","TTC","TTA","TTG","ATT","ATC","ATG","GTT","GTC","GTA","GTG",
                   "TCT","TCC","TCA","TCG","AGT","AGC","CCT","CCC","CCA","CCG",
                   "ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG",
                   "TAT","TAC","TAA","TAG","TGA","CAT","CAC","CAA","CAG","AAT","AAC",
                   "AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGG",
                   "CGT","CGC","CGA","CGG","AGA","AGG","GGT","GGC","GGA","GGG"
            ).forEach { codon -> mm[codon] = st.translate(codon) }
            mm["CTT"] = 'T'; mm["CTC"] = 'T'; mm["CTA"] = 'T'; mm["CTG"] = 'T'
            mm["ATA"] = 'M'; mm["TGA"] = 'W'
            mm
        }
        CodonTable("Yeast Mitochondrial", 3, m)
    }

    fun byId(id: Int): CodonTable = when (id) {
        1 -> STANDARD
        2 -> VERTEBRATE_MITOCHONDRIAL
        3 -> YEAST_MITOCHONDRIAL
        else -> throw IllegalArgumentException("Unknown codon table ID: $id")
    }
}
