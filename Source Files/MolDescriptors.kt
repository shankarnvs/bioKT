package biokt

import kotlin.math.*

// ============================================================
// MOLECULAR DESCRIPTORS
// Lipinski / Veber / ADMET / Fingerprints / Drug-likeness
// Mirrors and extends pyBioMed's Chem module
// ============================================================

// ─── Descriptor Result ───────────────────────────────────────

data class MolecularDescriptors(
    // Basic counts
    val molecularFormula: String,
    val exactMass: Double,
    val numAtoms: Int,
    val numHeavyAtoms: Int,
    val numBonds: Int,
    val numRings: Int,
    val numAromaticRings: Int,
    val numRotatableBonds: Int,
    val numStereocenters: Int,

    // Lipinski Rule of Five
    val molecularWeight: Double,
    val logP: Double,          // Wildman-Crippen estimate
    val hBondDonors: Int,
    val hBondAcceptors: Int,
    val lipinskiPasses: Boolean,

    // Veber rules (oral bioavailability)
    val topologicalPolarSurfaceArea: Double,   // TPSA (Å²)
    val veberPasses: Boolean,

    // Extended drug-likeness
    val molarRefractivity: Double,
    val fractionCSP3: Double,       // fraction of sp3 carbons
    val numAromaticAtoms: Int,

    // Ghose filter
    val ghosePasses: Boolean,

    // Egan filter
    val eganPasses: Boolean,

    // Muegge filter
    val mueggeScore: Int,   // number of rules passed (max 9)

    // ADMET estimates
    val admet: ADMETProfile
) {
    fun summary(): String = buildString {
        append("=== Molecular Descriptors ===\n")
        append("Formula     : $molecularFormula\n")
        append("Exact mass  : ${"%.4f".format(exactMass)} Da\n")
        append("MW          : ${"%.2f".format(molecularWeight)} Da\n")
        append("LogP        : ${"%.2f".format(logP)}\n")
        append("HBD / HBA   : $hBondDonors / $hBondAcceptors\n")
        append("TPSA        : ${"%.1f".format(topologicalPolarSurfaceArea)} Å²\n")
        append("Rot. bonds  : $numRotatableBonds\n")
        append("Rings       : $numRings (aromatic: $numAromaticRings)\n")
        append("fsp3        : ${"%.2f".format(fractionCSP3)}\n")
        append("Lipinski    : ${if (lipinskiPasses) "PASS ✓" else "FAIL ✗"}\n")
        append("Veber       : ${if (veberPasses) "PASS ✓" else "FAIL ✗"}\n")
        append("Ghose       : ${if (ghosePasses) "PASS ✓" else "FAIL ✗"}\n")
        append("\n--- ADMET ---\n")
        append(admet.summary())
    }
}

data class ADMETProfile(
    // Absorption
    val oralBioavailability: String,   // "High" / "Low" / "Moderate"
    val caco2Permeability: String,     // "High" / "Low"
    val pgpSubstrate: Boolean,
    val hia: String,                   // Human Intestinal Absorption

    // Distribution
    val bbbPenetration: Boolean,       // Blood-Brain Barrier
    val vdss: Double,                  // Volume of distribution (L/kg)
    val plasmaProteinBinding: String,  // "High >90%" / "Moderate" / "Low"

    // Metabolism
    val cyp3a4Substrate: Boolean,
    val cyp2d6Substrate: Boolean,
    val cyp2c9Substrate: Boolean,
    val cyp3a4Inhibitor: Boolean,
    val cyp2d6Inhibitor: Boolean,

    // Excretion
    val halfLife: String,              // "Short <2h" / "Moderate" / "Long >8h"
    val renalClearance: String,

    // Toxicity
    val hergInhibition: Boolean,       // Cardiotoxicity risk
    val amesTest: String,              // "Mutagenic" / "Non-mutagenic"
    val oralToxicity: String,          // LD50 class
    val skinSensitization: Boolean
) {
    fun summary() = buildString {
        append("Oral bioavailability : $oralBioavailability\n")
        append("BBB penetration      : ${if (bbbPenetration) "Yes" else "No"}\n")
        append("CYP3A4 substrate     : ${if (cyp3a4Substrate) "Yes" else "No"}\n")
        append("CYP2D6 substrate     : ${if (cyp2d6Substrate) "Yes" else "No"}\n")
        append("hERG inhibition      : ${if (hergInhibition) "Risk" else "Low risk"}\n")
        append("Ames test            : $amesTest\n")
        append("Oral toxicity        : $oralToxicity\n")
    }
}

// ─── Descriptor Calculator ───────────────────────────────────

object MolDescriptors {


    // ── Molecular Weight ─────────────────────────────────────

    fun molecularWeight(mol: Molecule): Double {
        val avgMasses = mapOf(
            "H" to 1.008, "C" to 12.011, "N" to 14.007, "O" to 15.999,
            "F" to 18.998, "P" to 30.974, "S" to 32.065, "Cl" to 35.453,
            "Br" to 79.904, "I" to 126.905, "B" to 10.811, "Si" to 28.086,
            "Se" to 78.971, "Fe" to 55.845, "Zn" to 65.38
        )
        var mw = 0.0
        for (atom in mol.atoms) mw += avgMasses[atom.symbol] ?: 12.0
        mw += mol.atoms.map { it.totalH * 1.008 }.sum()
        return mw
    }

    // ── Wildman-Crippen LogP ──────────────────────────────────
    // Contribution-based estimate using atom types

    fun wilmanCrippenLogP(mol: Molecule): Double {
        // Simplified atom contribution table (Wildman & Crippen 1999)
        val contributions = mapOf(
            "C_ar"  to  0.1441,  // aromatic C
            "C_sp3" to  0.1234,  // aliphatic C
            "C_sp2" to  0.0000,  // sp2 C (non-aromatic)
            "N_ar"  to -0.7736,  // aromatic N
            "N_am"  to -1.0190,  // amide N
            "N_sp3" to -0.9275,  // aliphatic N
            "O_ar"  to -0.2405,  // aromatic O
            "O_sp3" to -0.2677,  // aliphatic O
            "O_co2" to -0.0200,  // carboxylate O
            "F"     to  0.4202,
            "Cl"    to  0.7338,
            "Br"    to  0.8456,
            "I"     to  0.8857,
            "S"     to  0.6895,
            "P"     to  0.8612,
            "H"     to  0.1234
        )

        var logP = 0.0
        for (atom in mol.atoms) {
            val neighbors = mol.neighborsOf(atom.index)
            val hasO = neighbors.any { mol.atoms[it].symbol == "O" }
            val hasN = neighbors.any { mol.atoms[it].symbol == "N" }

            val contrib = when (atom.symbol) {
                "C"  -> when {
                    atom.isAromatic -> contributions["C_ar"]!!
                    mol.bonds.filter { it.atom1 == atom.index || it.atom2 == atom.index }
                        .any { it.type == BondType.DOUBLE } -> contributions["C_sp2"]!!
                    else -> contributions["C_sp3"]!!
                }
                "N"  -> when {
                    atom.isAromatic -> contributions["N_ar"]!!
                    hasO -> contributions["N_am"]!!
                    else -> contributions["N_sp3"]!!
                }
                "O"  -> when {
                    atom.isAromatic -> contributions["O_ar"]!!
                    neighbors.any { n -> mol.bonds.filter { it.atom1 == n || it.atom2 == n }
                        .any { it.type == BondType.DOUBLE } } -> contributions["O_co2"]!!
                    else -> contributions["O_sp3"]!!
                }
                "F"  -> contributions["F"]!!
                "Cl" -> contributions["Cl"]!!
                "Br" -> contributions["Br"]!!
                "I"  -> contributions["I"]!!
                "S"  -> contributions["S"]!!
                "P"  -> contributions["P"]!!
                "H"  -> 0.0
                else -> 0.0
            }
            logP += contrib
        }
        // Add H contributions
        logP += mol.atoms.map { it.totalH }.sum() * (contributions["H"] ?: 0.0)
        return logP
    }

    // ── H-Bond Donors / Acceptors ────────────────────────────

    fun hBondDonors(mol: Molecule): Int {
        // O-H and N-H groups
        return mol.atoms.count { atom ->
            (atom.symbol == "O" || atom.symbol == "N") && atom.totalH > 0
        }
    }

    fun hBondAcceptors(mol: Molecule): Int {
        // N and O atoms (Lipinski definition)
        return mol.atoms.count { it.symbol == "N" || it.symbol == "O" }
    }

    // ── TPSA (Ertl 2000) ─────────────────────────────────────

    fun tpsa(mol: Molecule): Double {
        // Fragment contributions (simplified Ertl table)
        var tpsa = 0.0
        for (atom in mol.atoms) {
            tpsa += when (atom.symbol) {
                "N" -> when {
                    atom.isAromatic && atom.totalH == 0 -> 12.89
                    atom.isAromatic && atom.totalH == 1 -> 13.59
                    !atom.isAromatic && atom.totalH == 0 -> 3.24
                    !atom.isAromatic && atom.totalH == 1 -> 23.85
                    !atom.isAromatic && atom.totalH == 2 -> 26.02
                    else -> 0.0
                }
                "O" -> when {
                    atom.isAromatic -> 13.14
                    atom.totalH == 1 -> 20.23   // O-H
                    else -> 17.07               // ether/carbonyl
                }
                "S" -> when {
                    atom.totalH == 1 -> 34.14
                    else -> 28.24
                }
                "P" -> when {
                    atom.totalH == 0 -> 9.81
                    else -> 23.47
                }
                else -> 0.0
            }
        }
        return tpsa
    }

    // ── Rotatable Bonds ──────────────────────────────────────

    fun rotatableBonds(mol: Molecule): Int {
        return mol.bonds.count { bond ->
            if (bond.type != BondType.SINGLE) return@count false
            if (bond.isRing) return@count false
            val a1 = mol.atoms[bond.atom1]; val a2 = mol.atoms[bond.atom2]
            if (a1.symbol == "H" || a2.symbol == "H") return@count false
            // Terminal atoms (degree 1) don't rotate meaningfully
            val deg1 = mol.neighborsOf(bond.atom1).size
            val deg2 = mol.neighborsOf(bond.atom2).size
            deg1 > 1 && deg2 > 1
        }
    }

    // ── Fraction sp3 carbons ─────────────────────────────────

    fun fractionCSP3(mol: Molecule): Double {
        val carbons = mol.atoms.filter { it.symbol == "C" }
        if (carbons.isEmpty()) return 0.0
        val sp3 = carbons.count { atom ->
            !atom.isAromatic && mol.bonds
                .filter { it.atom1 == atom.index || it.atom2 == atom.index }
                .none { it.type == BondType.DOUBLE || it.type == BondType.TRIPLE }
        }
        return sp3.toDouble() / carbons.size
    }

    // ── Molar Refractivity (Wildman-Crippen) ─────────────────

    fun molarRefractivity(mol: Molecule): Double {
        val mw = molecularWeight(mol)
        val logp = wilmanCrippenLogP(mol)
        // Approximate: MR ≈ 0.029MW + 0.985 * some factor
        return 0.029 * mw + 1.512 * mol.atoms.count { it.symbol == "C" || it.symbol == "N" } +
               0.686 * mol.atoms.count { it.symbol == "O" || it.symbol == "S" }
    }

    // ── Stereocenters ────────────────────────────────────────

    fun stereocenters(mol: Molecule): Int {
        return mol.atoms.count { atom ->
            if (atom.symbol != "C") return@count false
            val neighbors = mol.neighborsOf(atom.index)
            if (neighbors.size != 4) return@count false
            // Check all 4 substituents are different (simplified)
            val syms = neighbors.map { mol.atoms[it].symbol }
            syms.size == syms.toSet().size
        }
    }

    // ── Drug-Likeness Filters ────────────────────────────────

    fun lipinskiRo5(mw: Double, logP: Double, hbd: Int, hba: Int): Boolean =
        mw <= 500 && logP <= 5 && hbd <= 5 && hba <= 10

    fun veberRules(tpsa: Double, rotBonds: Int): Boolean =
        tpsa <= 140 && rotBonds <= 10

    fun ghoseFilter(mw: Double, logP: Double, mr: Double, numAtoms: Int): Boolean =
        mw in 160.0..480.0 && logP in -0.4..5.6 && mr in 40.0..130.0 && numAtoms in 20..70

    fun eganFilter(logP: Double, tpsa: Double): Boolean =
        logP in -1.0..6.0 && tpsa <= 150.0

    fun mueggeFilter(mol: Molecule, mw: Double, logP: Double, tpsa: Double, hbd: Int, hba: Int, rotB: Int): Int {
        var score = 0
        if (mw in 200.0..600.0) score++
        if (logP in -2.0..5.0) score++
        if (tpsa <= 150) score++
        if (hbd <= 5) score++
        if (hba <= 10) score++
        if (rotB <= 15) score++
        if (mol.rings.size <= 7) score++
        if (mol.rings.size >= 1) score++
        if (mol.atoms.none { it.symbol in setOf("Fe","Hg","Pb","Cd") }) score++
        return score
    }

    // ── ADMET Estimation ─────────────────────────────────────

    fun estimateADMET(mol: Molecule, mw: Double, logP: Double, tpsa: Double, hbd: Int, hba: Int, mr: Double): ADMETProfile {
        // Rule-based ADMET estimates (Egan, Veber, various literature rules)

        val oralBio = when {
            tpsa < 60 && logP in 0.0..3.0 -> "High"
            tpsa > 140 || logP > 5.0 || logP < -2.0 -> "Low"
            else -> "Moderate"
        }

        val caco2 = if (tpsa < 60 && mw < 500 && logP > 0) "High" else "Low"

        val pgp = logP > 3.0 && mw > 400 && mol.rings.size >= 2

        val hia = when {
            tpsa < 75 && logP in -1.0..4.0 -> "High (>80%)"
            tpsa > 120 -> "Low (<20%)"
            else -> "Moderate (20-80%)"
        }

        val bbb = logP in 1.0..3.5 && mw < 450 && tpsa < 90 && hbd <= 3 && hba <= 7

        val vdss = when {
            logP > 3.0 -> 2.5 + (logP - 3.0) * 0.5
            logP < 0.0 -> 0.2
            else -> 0.7 + logP * 0.6
        }

        val ppb = when {
            logP > 3.0 -> "High >90%"
            logP > 1.0 -> "Moderate 50-90%"
            else -> "Low <50%"
        }

        val cyp3a4Sub = mw > 400 && mol.rings.count { it.isAromatic } >= 1
        val cyp2d6Sub = mol.atoms.any { it.symbol == "N" && !it.isAromatic } && mw < 500
        val cyp2c9Sub = mol.atoms.any { it.symbol == "O" } && logP in 1.0..5.0
        val cyp3a4Inh = mol.rings.count { it.isAromatic } >= 2 && mw > 350
        val cyp2d6Inh = mol.atoms.any { it.symbol == "N" } && mol.rings.count { it.isAromatic } >= 1

        val halfLife = when {
            mw < 200 -> "Short <2h"
            mw > 500 -> "Long >8h"
            else -> "Moderate 2-8h"
        }

        val renalClear = if (logP < 1.0 && mw < 300) "High" else "Low/Moderate"

        // hERG inhibition risk (lipophilic + basic nitrogen = risk)
        val herg = logP > 2.0 && mol.atoms.any { it.symbol == "N" && !it.isAromatic }

        // Ames mutagenicity (simplified — aromatic amines, nitroso groups)
        val ames = if (mol.atoms.any { it.symbol == "N" && it.isAromatic } &&
                       mol.rings.count { it.isAromatic } >= 2) "Potentially mutagenic" else "Non-mutagenic"

        // Oral toxicity (Cramer class rules — simplified)
        val oralTox = when {
            mw < 200 && logP < 1.0 -> "Class I (Low risk)"
            mol.rings.count { it.isAromatic } >= 3 -> "Class III (High risk)"
            else -> "Class II (Moderate risk)"
        }

        val skinSens = mol.atoms.any { it.symbol == "S" } && logP > 2.0

        return ADMETProfile(oralBio, caco2, pgp, hia, bbb, vdss, ppb,
            cyp3a4Sub, cyp2d6Sub, cyp2c9Sub, cyp3a4Inh, cyp2d6Inh,
            halfLife, renalClear, herg, ames, oralTox, skinSens)
    }

    // ── Structural Alerts (PAINS, Brenk) ────────────────────

    data class StructuralAlert(val name: String, val description: String, val severity: String)

    fun checkStructuralAlerts(mol: Molecule): List<StructuralAlert> {
        val alerts = mutableListOf<StructuralAlert>()

        // Reactive groups
        if (mol.bonds.any { it.type == BondType.DOUBLE } &&
            mol.atoms.any { it.symbol == "O" }) {
            val carbonylC = mol.atoms.filter { it.symbol == "C" }
                .any { c -> mol.bonds.filter { it.atom1 == c.index || it.atom2 == c.index }
                    .any { b -> b.type == BondType.DOUBLE &&
                        mol.atoms[if (b.atom1 == c.index) b.atom2 else b.atom1].symbol == "O" } }
            if (carbonylC && mol.bonds.any { b ->
                    b.type == BondType.DOUBLE &&
                    mol.atoms[b.atom1].symbol == "C" && mol.atoms[b.atom2].symbol == "C" }) {
                alerts.add(StructuralAlert("Michael Acceptor", "α,β-unsaturated carbonyl — can alkylate nucleophiles", "High"))
            }
        }

        // Too many halogens
        val halogens = mol.atoms.count { it.symbol in setOf("F", "Cl", "Br", "I") }
        if (halogens >= 4)
            alerts.add(StructuralAlert("Polyhalogenated", "≥4 halogen atoms — potential toxicity", "Moderate"))

        // High MW / complexity
        if (mol.atoms.filter { it.symbol != "H" }.size > 50)
            alerts.add(StructuralAlert("High Complexity", "More than 50 heavy atoms — poor oral absorption likely", "Low"))

        // Nitro groups (aromatic nitro = mutagenicity risk)
        val nitrogenCount = mol.atoms.count { it.symbol == "N" }
        val oxygenCount = mol.atoms.count { it.symbol == "O" }
        if (nitrogenCount >= 1 && oxygenCount >= 2 && mol.rings.count { it.isAromatic } >= 1)
            alerts.add(StructuralAlert("Potential Nitro Group", "Aromatic nitro compounds — mutagenicity risk", "High"))

        return alerts
    }

    // ── Full descriptor calculation ──────────────────────────

    fun calculate(mol: Molecule): MolecularDescriptors {
        val mw   = molecularWeight(mol)
        val logp = wilmanCrippenLogP(mol)
        val hbd  = hBondDonors(mol)
        val hba  = hBondAcceptors(mol)
        val tpsa = tpsa(mol)
        val rotB = rotatableBonds(mol)
        val fsp3 = fractionCSP3(mol)
        val mr   = molarRefractivity(mol)
        val nRings   = mol.rings.size
        val nArRings = mol.rings.count { it.isAromatic }
        val nAtoms   = mol.atoms.filter { it.symbol != "H" }.size
        val nHAtoms  = mol.atoms.size + mol.atoms.map { it.totalH }.sum()
        val mass     = mol.atoms.map { it.atomicMass }.sum() + mol.atoms.map { it.totalH }.sum() * 1.008

        return MolecularDescriptors(
            molecularFormula    = mol.molecularFormula(),
            exactMass           = mass,
            numAtoms            = nHAtoms,
            numHeavyAtoms       = nAtoms,
            numBonds            = mol.bonds.size,
            numRings            = nRings,
            numAromaticRings    = nArRings,
            numRotatableBonds   = rotB,
            numStereocenters    = stereocenters(mol),
            molecularWeight     = mw,
            logP                = logp,
            hBondDonors         = hbd,
            hBondAcceptors      = hba,
            lipinskiPasses      = lipinskiRo5(mw, logp, hbd, hba),
            topologicalPolarSurfaceArea = tpsa,
            veberPasses         = veberRules(tpsa, rotB),
            molarRefractivity   = mr,
            fractionCSP3        = fsp3,
            numAromaticAtoms    = mol.atoms.count { it.isAromatic },
            ghosePasses         = ghoseFilter(mw, logp, mr, nHAtoms),
            eganPasses          = eganFilter(logp, tpsa),
            mueggeScore         = mueggeFilter(mol, mw, logp, tpsa, hbd, hba, rotB),
            admet               = estimateADMET(mol, mw, logp, tpsa, hbd, hba, mr)
        )
    }

}

// ─── Molecular Fingerprints ──────────────────────────────────

object Fingerprints {

    // Morgan/ECFP-like circular fingerprint (bit vector)
    fun morgan(mol: Molecule, radius: Int = 2, nBits: Int = 1024): BooleanArray {
        val fp = BooleanArray(nBits)
        for (atom in mol.atoms) {
            var hash = atomHash(atom)
            fp[Math.abs(hash) % nBits] = true
            for (r in 1..radius) {
                val neighbors = mol.neighborsOf(atom.index)
                    .map { mol.atoms[it] }
                    .sortedBy { it.atomicNum }
                hash = hash * 31 + neighbors.fold(0) { acc, a -> acc * 31 + atomHash(a) }
                fp[Math.abs(hash) % nBits] = true
            }
        }
        return fp
    }

    // MACCS keys (simplified 166-bit)
    fun maccs(mol: Molecule): BooleanArray {
        val fp = BooleanArray(166)
        // A selection of the 166 MACCS keys
        fp[1]  = mol.atoms.any { it.symbol == "Cl" || it.symbol == "Br" }
        fp[4]  = mol.atoms.any { it.symbol == "N" && it.isAromatic }
        fp[6]  = mol.atoms.any { it.symbol == "O" && !it.isAromatic }
        fp[7]  = mol.atoms.any { it.symbol == "N" && !it.isAromatic }
        fp[14] = mol.rings.size > 0
        fp[16] = mol.rings.count { it.isAromatic } > 0
        fp[18] = mol.atoms.any { it.symbol == "S" }
        fp[24] = mol.atoms.any { it.symbol == "F" }
        fp[27] = mol.rings.count { it.isAromatic } >= 2
        fp[32] = mol.bonds.any { it.type == BondType.TRIPLE }
        fp[44] = mol.atoms.any { it.symbol == "I" }
        fp[65] = mol.rings.count { it.isAromatic } >= 3
        fp[91] = mol.atoms.count { it.symbol == "O" } >= 2
        fp[120]= mol.atoms.count { it.symbol == "N" } >= 2
        fp[160]= mol.atoms.any { it.symbol == "P" }
        fp[162]= mol.bonds.any { it.type == BondType.DOUBLE }
        return fp
    }

    // Topological (RDKit-style path fingerprint)
    fun topological(mol: Molecule, nBits: Int = 2048, maxPath: Int = 7): BooleanArray {
        val fp = BooleanArray(nBits)
        // Enumerate all paths up to maxPath length via DFS
        fun dfs(path: List<Int>) {
            if (path.size > 1) {
                val hash = path.map { mol.atoms[it].atomicNum }.hashCode()
                fp[Math.abs(hash) % nBits] = true
            }
            if (path.size < maxPath) {
                val last = path.last()
                for (nb in mol.neighborsOf(last)) {
                    if (nb !in path) dfs(path + nb)
                }
            }
        }
        for (i in mol.atoms.indices) dfs(listOf(i))
        return fp
    }

    // Tanimoto similarity between two fingerprints
    fun tanimoto(fp1: BooleanArray, fp2: BooleanArray): Double {
        require(fp1.size == fp2.size) { "Fingerprints must be same length" }
        val a = fp1.count { it }
        val b = fp2.count { it }
        val c = fp1.indices.count { fp1[it] && fp2[it] }
        return if (a + b - c == 0) 0.0 else c.toDouble() / (a + b - c)
    }

    // Dice similarity
    fun dice(fp1: BooleanArray, fp2: BooleanArray): Double {
        require(fp1.size == fp2.size)
        val a = fp1.count { it }
        val b = fp2.count { it }
        val c = fp1.indices.count { fp1[it] && fp2[it] }
        return if (a + b == 0) 0.0 else 2.0 * c / (a + b)
    }

    // Cosine similarity
    fun cosine(fp1: BooleanArray, fp2: BooleanArray): Double {
        require(fp1.size == fp2.size)
        val dot = fp1.indices.count { fp1[it] && fp2[it] }.toDouble()
        val n1 = sqrt(fp1.count { it }.toDouble())
        val n2 = sqrt(fp2.count { it }.toDouble())
        return if (n1 == 0.0 || n2 == 0.0) 0.0 else dot / (n1 * n2)
    }

    private fun atomHash(atom: Atom): Int =
        atom.atomicNum * 31 + atom.charge * 17 + (if (atom.isAromatic) 7 else 0) + atom.totalH * 3
}

// ─── Virtual Screening ───────────────────────────────────────

data class VirtualScreeningResult(
    val molecule: Molecule,
    val tanimoto: Double,
    val descriptors: MolecularDescriptors,
    val passesLipinski: Boolean,
    val rank: Int
)

object VirtualScreening {

    fun screenLibrary(
        query: Molecule,
        library: List<Molecule>,
        fingerprintType: String = "morgan",
        tanimotoThreshold: Double = 0.4,
        applyLipinski: Boolean = true,
        topN: Int = 20
    ): List<VirtualScreeningResult> {
        val queryFp = when (fingerprintType) {
            "morgan" -> Fingerprints.morgan(query)
            "maccs"  -> Fingerprints.maccs(query)
            "topo"   -> Fingerprints.topological(query)
            else     -> Fingerprints.morgan(query)
        }

        val results = library.mapIndexedNotNull { i, mol ->
            val fp = when (fingerprintType) {
                "morgan" -> Fingerprints.morgan(mol)
                "maccs"  -> Fingerprints.maccs(mol)
                "topo"   -> Fingerprints.topological(mol)
                else     -> Fingerprints.morgan(mol)
            }
            val sim = Fingerprints.tanimoto(queryFp, fp)
            if (sim < tanimotoThreshold) return@mapIndexedNotNull null
            val desc = MolDescriptors.calculate(mol)
            if (applyLipinski && !desc.lipinskiPasses) return@mapIndexedNotNull null
            Triple(mol, sim, desc)
        }.sortedByDescending { it.second }.take(topN)

        return results.mapIndexed { i, (mol, sim, desc) ->
            VirtualScreeningResult(mol, sim, desc, desc.lipinskiPasses, i + 1)
        }
    }

    fun diversityPicking(molecules: List<Molecule>, n: Int, fingerprintType: String = "morgan"): List<Molecule> {
        if (molecules.size <= n) return molecules
        val fps = molecules.map { Fingerprints.morgan(it) }
        val selected = mutableListOf(0)
        val remaining = (1 until molecules.size).toMutableList()

        while (selected.size < n && remaining.isNotEmpty()) {
            // MaxMin diversity: pick molecule most dissimilar to already selected
            val next = remaining.maxBy { idx ->
                selected.map { sel -> 1.0 - Fingerprints.tanimoto(fps[idx], fps[sel]) }.min()!!
            }!!
            selected.add(next)
            remaining.remove(next)
        }
        return selected.map { molecules[it] }
    }
}
