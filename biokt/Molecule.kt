package biokt

import kotlin.math.*

// ============================================================
// MOLECULE — Small molecule / drug representation
// SMILES parsing, molecular graph, physicochemical properties
// ============================================================

// ─── Atom ────────────────────────────────────────────────────

data class Atom(
    val index: Int,
    val symbol: String,
    val atomicNum: Int,
    val charge: Int = 0,
    val isotope: Int = 0,
    val isAromatic: Boolean = false,
    var implicitH: Int = 0,
    var explicitH: Int = 0
) {
    val totalH: Int get() = implicitH + explicitH
    val atomicMass: Double get() = ATOMIC_MASSES[symbol] ?: 0.0
    val electronegativity: Double get() = ELECTRONEGATIVITIES[symbol] ?: 0.0
    val vdwRadius: Double get() = VDW_RADII[symbol] ?: 1.7

    companion object {
        val ATOMIC_MASSES = mapOf(
            "H" to 1.008, "C" to 12.011, "N" to 14.007, "O" to 15.999,
            "F" to 18.998, "P" to 30.974, "S" to 32.06, "Cl" to 35.45,
            "Br" to 79.904, "I" to 126.904, "B" to 10.81, "Si" to 28.085,
            "Se" to 78.971, "Fe" to 55.845, "Zn" to 65.38, "Cu" to 63.546,
            "Mg" to 24.305, "Ca" to 40.078, "Na" to 22.990, "K" to 39.098
        )
        val ELECTRONEGATIVITIES = mapOf(
            "H" to 2.20, "C" to 2.55, "N" to 3.04, "O" to 3.44,
            "F" to 3.98, "P" to 2.19, "S" to 2.58, "Cl" to 3.16,
            "Br" to 2.96, "I" to 2.66, "B" to 2.04, "Si" to 1.90
        )
        val VDW_RADII = mapOf(
            "H" to 1.20, "C" to 1.70, "N" to 1.55, "O" to 1.52,
            "F" to 1.47, "P" to 1.80, "S" to 1.80, "Cl" to 1.75,
            "Br" to 1.85, "I" to 1.98
        )
        val VALENCES = mapOf(
            "H" to 1, "C" to 4, "N" to 3, "O" to 2, "F" to 1,
            "P" to 5, "S" to 6, "Cl" to 1, "Br" to 1, "I" to 1,
            "B" to 3, "Si" to 4, "Se" to 2
        )
    }
}

// ─── Bond ────────────────────────────────────────────────────

enum class BondType(val order: Double) {
    SINGLE(1.0), DOUBLE(2.0), TRIPLE(3.0), AROMATIC(1.5)
}

data class Bond(
    val atom1: Int,
    val atom2: Int,
    val type: BondType,
    val isRing: Boolean = false
) {
    val order: Double get() = type.order
}

// ─── Ring ────────────────────────────────────────────────────

data class Ring(
    val atomIndices: List<Int>,
    val isAromatic: Boolean
) {
    val size: Int get() = atomIndices.size
}

// ─── Molecule ────────────────────────────────────────────────

class Molecule(
    val smiles: String,
    val name: String = "",
    val id: String = ""
) {
    val atoms: MutableList<Atom> = mutableListOf()
    val bonds: MutableList<Bond> = mutableListOf()
    private val adjacency: MutableMap<Int, MutableList<Int>> = mutableMapOf()
    var rings: List<Ring> = emptyList()
        private set

    init {
        SmilesParser.parse(smiles, this)
        detectRings()
        assignImplicitHydrogens()
    }

    // ── Graph queries ────────────────────────────────────────

    fun neighborsOf(atomIdx: Int): List<Int> = adjacency[atomIdx] ?: emptyList()

    fun bondBetween(a: Int, b: Int): Bond? =
        bonds.find { (it.atom1 == a && it.atom2 == b) || (it.atom1 == b && it.atom2 == a) }

    fun addAtom(atom: Atom) {
        atoms.add(atom)
        adjacency[atom.index] = mutableListOf()
    }

    fun addBond(bond: Bond) {
        bonds.add(bond)
        adjacency.getOrPut(bond.atom1) { mutableListOf() }.add(bond.atom2)
        adjacency.getOrPut(bond.atom2) { mutableListOf() }.add(bond.atom1)
    }

    // ── Ring detection (SSSR via DFS) ────────────────────────

    private fun detectRings() {
        val found = mutableListOf<Ring>()
        val visited = BooleanArray(atoms.size)
        val path = mutableListOf<Int>()

        fun dfs(node: Int, parent: Int) {
            visited[node] = true
            path.add(node)
            for (nb in neighborsOf(node)) {
                if (nb == parent) continue
                if (visited[nb]) {
                    val ringStart = path.indexOf(nb)
                    if (ringStart >= 0) {
                        val ringAtoms = path.subList(ringStart, path.size).toList()
                        val aromatic = ringAtoms.all { atoms[it].isAromatic }
                        found.add(Ring(ringAtoms, aromatic))
                    }
                } else {
                    dfs(nb, node)
                }
            }
            path.removeAt(path.size - 1)
        }

        for (i in atoms.indices) if (!visited[i]) dfs(i, -1)
        rings = found.distinctBy { it.atomIndices.toSet() }
    }

    // ── Implicit H assignment ────────────────────────────────

    private fun assignImplicitHydrogens() {
        for (atom in atoms) {
            if (atom.symbol == "H") continue
            val valence = Atom.VALENCES[atom.symbol] ?: continue
            val bondOrder = bonds.filter { it.atom1 == atom.index || it.atom2 == atom.index }
                .map { it.order.toInt() }.sum()
            val implicit = maxOf(0, valence - bondOrder - atom.charge + atom.explicitH)
            atoms[atom.index] = atom.copy(implicitH = implicit)
        }
    }

    // ── Molecular formula ────────────────────────────────────

    fun molecularFormula(): String {
        val counts = mutableMapOf<String, Int>()
        for (atom in atoms) {
            counts[atom.symbol] = (counts[atom.symbol] ?: 0) + 1
            if (atom.totalH > 0)
                counts["H"] = (counts["H"] ?: 0) + atom.totalH
        }
        return buildString {
            // Hill order: C first, then H, then alphabetical
            listOf("C", "H").forEach { s -> counts[s]?.let { append(s); if (it > 1) append(it) } }
            counts.keys.filter { it != "C" && it != "H" }.sorted().forEach { s ->
                val n = counts[s]!!; append(s); if (n > 1) append(n)
            }
        }
    }

    fun exactMass(): Double {
        var mass = atoms.map { it.atomicMass }.sum()
        mass += atoms.map { it.totalH }.sum() * (Atom.ATOMIC_MASSES["H"] ?: 1.008)
        return mass
    }

    fun numAtoms(includeH: Boolean = false): Int =
        if (includeH) atoms.size + atoms.map { it.totalH }.sum() else atoms.filter { it.symbol != "H" }.size

    fun numBonds(): Int = bonds.size
    fun numRings(): Int = rings.size
    fun numAromaticRings(): Int = rings.count { it.isAromatic }

    override fun toString() = "Molecule(name='$name', formula=${molecularFormula()}, atoms=${atoms.size})"
}

// ─── SMILES Parser ───────────────────────────────────────────

object SmilesParser {

    private val ORGANIC_SUBSET = setOf("B","C","N","O","P","S","F","Cl","Br","I")
    private val TWO_CHAR = setOf("Cl","Br","Si","Se","Na","Mg","Ca","Fe","Zn","Cu")

    fun parse(smiles: String, mol: Molecule) {
        val s = smiles.trim()
        var i = 0
        val atomStack = java.util.LinkedList<Int>()
        val ringOpens = mutableMapOf<Int, Pair<Int, BondType>>() // ringNum -> (atomIdx, bondType)
        var pendingBond: BondType? = null
        var atomIdx = 0

        fun addAtomNode(symbol: String, aromatic: Boolean, charge: Int = 0, isotope: Int = 0, explicitH: Int = 0): Int {
            val aNum = ATOMIC_NUMBERS[symbol] ?: 0
            val atom = Atom(atomIdx, symbol, aNum, charge, isotope, aromatic, explicitH = explicitH)
            mol.addAtom(atom)
            return atomIdx++
        }

        fun connectToPrev(newIdx: Int, bond: BondType) {
            if (atomStack.isNotEmpty()) {
                mol.addBond(Bond(atomStack.last, newIdx, bond))
            }
        }

        while (i < s.length) {
            when {
                s[i] == '(' -> {
                    atomStack.addLast(atomStack.last)
                    pendingBond = null
                    i++
                }
                s[i] == ')' -> {
                    if (atomStack.isEmpty()) null else atomStack.removeLast()
                    pendingBond = null
                    i++
                }
                s[i] == '=' -> { pendingBond = BondType.DOUBLE; i++ }
                s[i] == '#' -> { pendingBond = BondType.TRIPLE; i++ }
                s[i] == ':' -> { pendingBond = BondType.AROMATIC; i++ }
                s[i] == '-' -> { pendingBond = BondType.SINGLE; i++ }
                s[i] == '.' -> { atomStack.clear(); pendingBond = null; i++ }

                s[i] == '[' -> {
                    val end = s.indexOf(']', i)
                    val inner = s.substring(i + 1, end)
                    val (sym, charge, hCount, isotope, aromatic) = parseBracketAtom(inner)
                    val idx = addAtomNode(sym, aromatic, charge, isotope, hCount)
                    val bond = pendingBond ?: if (aromatic && atomStack.isNotEmpty() && mol.atoms.getOrNull(atomStack.last)?.isAromatic == true) BondType.AROMATIC else BondType.SINGLE
                    connectToPrev(idx, bond)
                    atomStack.addLast(idx)
                    pendingBond = null
                    i = end + 1
                }

                s[i] == '%' && i + 2 < s.length -> {
                    val rnum = s.substring(i + 1, i + 3).toIntOrNull() ?: 0
                    handleRing(rnum, atomStack, ringOpens, mol, pendingBond)
                    pendingBond = null; i += 3
                }

                s[i].isDigit() -> {
                    val rnum = s[i].toString().toInt()
                    handleRing(rnum, atomStack, ringOpens, mol, pendingBond)
                    pendingBond = null; i++
                }

                else -> {
                    // Try two-char symbol
                    val twoChar = if (i + 1 < s.length) s.substring(i, i + 2) else ""
                    val (sym, len, aromatic) = when {
                        twoChar in TWO_CHAR -> Triple(twoChar, 2, false)
                        s[i].isUpperCase() -> Triple(s[i].toString(), 1, false)
                        s[i].isLowerCase() -> Triple(s[i].toUpperCase().toString(), 1, true)
                        else -> Triple("?", 1, false)
                    }
                    val idx = addAtomNode(sym, aromatic)
                    val bond = pendingBond ?: if (aromatic && atomStack.isNotEmpty() && mol.atoms.getOrNull(atomStack.last)?.isAromatic == true) BondType.AROMATIC else BondType.SINGLE
                    connectToPrev(idx, bond)
                    atomStack.addLast(idx)
                    pendingBond = null
                    i += len
                }
            }
        }
    }

    data class BracketAtomInfo(val symbol: String, val charge: Int, val hCount: Int, val isotope: Int, val aromatic: Boolean)

    private fun parseBracketAtom(inner: String): BracketAtomInfo {
        var s = inner
        var isotope = 0
        var charge = 0
        var hCount = 0

        // Isotope
        val isotopeMatch = Regex("^(\\d+)").find(s)
        if (isotopeMatch != null) { isotope = isotopeMatch.value.toInt(); s = s.drop(isotopeMatch.value.length) }

        // Symbol
        val symMatch = Regex("^([A-Z][a-z]?|[a-z])").find(s) ?: return BracketAtomInfo("C", 0, 0, 0, false)
        val rawSym = symMatch.value
        val aromatic = rawSym[0].isLowerCase()
        val symbol = rawSym[0].toUpperCase().toString() + rawSym.drop(1)
        s = s.drop(rawSym.length)

        // H count
        if (s.startsWith("H")) {
            s = s.drop(1)
            val hMatch = Regex("^(\\d+)").find(s)
            hCount = if (hMatch != null) { s = s.drop(hMatch.value.length); hMatch.value.toInt() } else 1
        }

        // Charge
        val chargeMatch = Regex("([+-])(\\d*)").find(s)
        if (chargeMatch != null) {
            val sign = if (chargeMatch.groupValues[1] == "+") 1 else -1
            val mag = chargeMatch.groupValues[2].toIntOrNull() ?: 1
            charge = sign * mag
        }

        return BracketAtomInfo(symbol, charge, hCount, isotope, aromatic)
    }

    private fun handleRing(rnum: Int, stack: java.util.LinkedList<Int>, opens: MutableMap<Int, Pair<Int, BondType>>, mol: Molecule, pending: BondType?) {
        val current = stack.lastOrNull() ?: return
        if (opens.containsKey(rnum)) {
            val (openAtom, openBond) = opens.remove(rnum)!!
            val bond = pending ?: openBond
            mol.addBond(Bond(openAtom, current, bond, isRing = true))
        } else {
            opens[rnum] = current to (pending ?: BondType.SINGLE)
        }
    }

    val ATOMIC_NUMBERS = mapOf(
        "H" to 1, "He" to 2, "Li" to 3, "Be" to 4, "B" to 5, "C" to 6,
        "N" to 7, "O" to 8, "F" to 9, "Ne" to 10, "Na" to 11, "Mg" to 12,
        "Al" to 13, "Si" to 14, "P" to 15, "S" to 16, "Cl" to 17, "Ar" to 18,
        "K" to 19, "Ca" to 20, "Fe" to 26, "Cu" to 29, "Zn" to 30,
        "Br" to 35, "Se" to 34, "I" to 53
    )
}
