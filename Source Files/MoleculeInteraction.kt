package biokt

import kotlin.math.*

internal fun Double.fmt(decimals: Int) = "%.${decimals}f".format(this)

// ============================================================
// MOLECULE–SEQUENCE INTERACTIONS
// Drug/ligand binding to DNA, RNA, and Protein
// Docking scores, binding affinity, pharmacophore
// ============================================================

// ─── Binding Result ──────────────────────────────────────────

data class BindingResult(
    val molecule: Molecule,
    val target: String,          // target ID
    val bindingScore: Double,    // lower = better (like docking score, kcal/mol)
    val estimatedKd: Double,     // dissociation constant (μM)
    val bindingSite: BindingSite?,
    val interactions: List<MolecularInteraction>,
    val confidence: Double
) {
    val isLikelyBinder: Boolean get() = bindingScore < -5.0 && confidence > 0.4
    override fun toString() =
        "Binding(score=%.2f kcal/mol, Kd=%.2f μM, likely=$isLikelyBinder)".format(bindingScore, estimatedKd)
}

data class BindingSite(
    val residues: List<Int>,   // residue indices involved
    val center: Triple<Double, Double, Double>,
    val volume: Double,        // Å³
    val hydrophobicity: Double
)

data class MolecularInteraction(
    val type: InteractionType,
    val strength: Double,         // kcal/mol contribution
    val description: String
)

enum class InteractionType {
    HYDROGEN_BOND, HYDROPHOBIC, ELECTROSTATIC, PI_STACKING,
    VAN_DER_WAALS, HALOGEN_BOND, METAL_COORDINATION, CATION_PI
}

// ─── Drug–Protein Binding ────────────────────────────────────

object DrugProteinBinding {

    // ── Empirical scoring function (ChemScore-like) ──────────

    fun score(drug: Molecule, protein: ProteinSequence): BindingResult {
        val desc = MolDescriptors.calculate(drug)
        val interactions = mutableListOf<MolecularInteraction>()
        var totalScore = 0.0

        // Hydrogen bonding contribution
        val hbScore = estimateHBondScore(drug, protein)
        totalScore += hbScore
        if (hbScore < -1.0) interactions.add(MolecularInteraction(
            InteractionType.HYDROGEN_BOND, hbScore,
            "${desc.hBondDonors + desc.hBondAcceptors} H-bond donors/acceptors"))

        // Hydrophobic contribution
        val hydrophobicScore = estimateHydrophobicScore(drug, protein, desc)
        totalScore += hydrophobicScore
        if (hydrophobicScore < -1.0) interactions.add(MolecularInteraction(
            InteractionType.HYDROPHOBIC, hydrophobicScore,
            "LogP=${"%.2f".format(desc.logP)}, hydrophobic contacts"))

        // Electrostatic contribution
        val elecScore = estimateElectrostaticScore(drug, protein)
        totalScore += elecScore
        if (abs(elecScore) > 0.5) interactions.add(MolecularInteraction(
            InteractionType.ELECTROSTATIC, elecScore,
            "Ionic / polar interactions"))

        // Aromatic stacking
        val aromaticScore = estimateAromaticScore(drug, protein)
        totalScore += aromaticScore
        if (aromaticScore < -0.5) interactions.add(MolecularInteraction(
            InteractionType.PI_STACKING, aromaticScore,
            "${drug.numAromaticRings()} aromatic rings, π-stacking"))

        // Desolvation penalty
        val desolv = desc.topologicalPolarSurfaceArea * 0.02
        totalScore += desolv
        interactions.add(MolecularInteraction(InteractionType.VAN_DER_WAALS, desolv,
            "Desolvation penalty (TPSA=${"%.1f".format(desc.topologicalPolarSurfaceArea)} Å²)"))

        // Rotatable bond entropy penalty
        val rotPenalty = desc.numRotatableBonds * 0.3
        totalScore += rotPenalty

        // Convert score to Kd using linear approximation (ΔG = RT ln Kd)
        val rt = 0.001987 * 310.15  // RT at 37°C in kcal/mol
        val kd = exp(totalScore / rt) * 1e6  // μM

        // Binding site prediction from surface-exposed residues
        val bindingSite = predictBindingSite(protein)

        val confidence = minOf(1.0, maxOf(0.0, (-totalScore - 2.0) / 8.0))

        return BindingResult(drug, protein.id, totalScore, kd, bindingSite, interactions, confidence)
    }

    private fun estimateHBondScore(drug: Molecule, protein: ProteinSequence): Double {
        val drugHBD = MolDescriptors.hBondDonors(drug)
        val drugHBA = MolDescriptors.hBondAcceptors(drug)
        val proteinHBD = protein.sequence.count { it in setOf('S','T','Y','K','R','H','N','Q','W') }
        val proteinHBA = protein.sequence.count { it in setOf('D','E','S','T','Y','N','Q','H') }

        val donorPairs  = minOf(drugHBD, proteinHBA)
        val acceptorPairs = minOf(drugHBA, proteinHBD)
        return -(donorPairs + acceptorPairs) * 0.9  // ~0.9 kcal/mol per H-bond
    }

    private fun estimateHydrophobicScore(drug: Molecule, protein: ProteinSequence, desc: MolecularDescriptors): Double {
        val hydrophobicProtein = protein.sequence.count { it in setOf('A','V','I','L','M','F','W','P','Y') }
        val hydrophobicFraction = hydrophobicProtein.toDouble() / protein.length
        return if (desc.logP > 0) -desc.logP * hydrophobicFraction * 0.8 else 0.0
    }

    private fun estimateElectrostaticScore(drug: Molecule, protein: ProteinSequence): Double {
        val chargedAtoms = drug.atoms.count { it.charge != 0 }
        val chargedResidues = protein.sequence.count { it in setOf('D','E','K','R','H') }
        val chargeProduct = drug.atoms.map { it.charge.toDouble() }.sum() *
                           ProteinDescriptors.chargeAtPH(protein, 7.0)
        return if (chargeProduct < 0) chargedAtoms * -0.5 else chargedAtoms * 0.2
    }

    private fun estimateAromaticScore(drug: Molecule, protein: ProteinSequence): Double {
        val aromatic = drug.numAromaticRings()
        val aromaticResidues = protein.sequence.count { it in setOf('F','Y','W','H') }
        return -(aromatic * aromaticResidues * 0.3).coerceAtMost(3.0)
    }

    private fun predictBindingSite(protein: ProteinSequence): BindingSite {
        val hydrophobicRes = protein.sequence.indices.filter {
            protein[it] in setOf('A','V','I','L','M','F','W','P')
        }.take(10)
        return BindingSite(
            residues = hydrophobicRes,
            center = Triple(0.0, 0.0, 0.0),
            volume = protein.length * 2.5,
            hydrophobicity = ProteinDescriptors.hydrophobicRatio(protein)
        )
    }

    // ── Pharmacophore Model ──────────────────────────────────

    data class Pharmacophore(
        val features: List<PharmacophoreFeature>
    )

    data class PharmacophoreFeature(
        val type: String,    // "HBD", "HBA", "Hydrophobic", "Aromatic", "Positive", "Negative"
        val atomIndices: List<Int>
    )

    fun generatePharmacophore(drug: Molecule): Pharmacophore {
        val features = mutableListOf<PharmacophoreFeature>()

        // H-bond donors (N-H, O-H)
        val hbd = drug.atoms.indices.filter { i ->
            val a = drug.atoms[i]
            (a.symbol == "N" || a.symbol == "O") && a.totalH > 0
        }
        if (hbd.isNotEmpty()) features.add(PharmacophoreFeature("HBD", hbd))

        // H-bond acceptors (N, O)
        val hba = drug.atoms.indices.filter { drug.atoms[it].symbol in setOf("N", "O") }
        if (hba.isNotEmpty()) features.add(PharmacophoreFeature("HBA", hba))

        // Hydrophobic (aliphatic C)
        val hydro = drug.atoms.indices.filter { i ->
            drug.atoms[i].symbol == "C" && !drug.atoms[i].isAromatic
        }
        if (hydro.isNotEmpty()) features.add(PharmacophoreFeature("Hydrophobic", hydro))

        // Aromatic rings
        val aromatic = drug.atoms.indices.filter { drug.atoms[it].isAromatic }
        if (aromatic.isNotEmpty()) features.add(PharmacophoreFeature("Aromatic", aromatic))

        // Positive charge (N+)
        val pos = drug.atoms.indices.filter { drug.atoms[it].symbol == "N" && drug.atoms[it].charge > 0 }
        if (pos.isNotEmpty()) features.add(PharmacophoreFeature("Positive", pos))

        // Negative charge (O-, carboxylate)
        val neg = drug.atoms.indices.filter { drug.atoms[it].symbol == "O" && drug.atoms[it].charge < 0 }
        if (neg.isNotEmpty()) features.add(PharmacophoreFeature("Negative", neg))

        return Pharmacophore(features)
    }

    fun pharmacophoreMatch(query: Pharmacophore, candidate: Molecule): Double {
        val candPharm = generatePharmacophore(candidate)
        var matches = 0
        for (feat in query.features) {
            if (candPharm.features.any { it.type == feat.type }) matches++
        }
        return if (query.features.isEmpty()) 0.0 else matches.toDouble() / query.features.size
    }

    // ── QSAR model (linear regression on descriptors) ────────

    data class QSARModel(
        val coefficients: Map<String, Double>,
        val intercept: Double,
        val r2: Double
    ) {
        fun predict(drug: Molecule): Double {
            val desc = MolDescriptors.calculate(drug)
            val features = mapOf(
                "mw" to desc.molecularWeight, "logP" to desc.logP,
                "hbd" to desc.hBondDonors.toDouble(), "hba" to desc.hBondAcceptors.toDouble(),
                "tpsa" to desc.topologicalPolarSurfaceArea, "rotB" to desc.numRotatableBonds.toDouble(),
                "rings" to desc.numRings.toDouble(), "mr" to desc.molarRefractivity,
                "fsp3" to desc.fractionCSP3
            )
            return intercept + coefficients.entries.map { (k, v) -> (features[k] ?: 0.0) * v }.sum()
        }
    }

    fun trainQSAR(drugs: List<Molecule>, activities: List<Double>): QSARModel {
        require(drugs.size == activities.size) { "Drugs and activities must be same size" }
        // Simple multiple linear regression via least squares (simplified)
        val descriptors = drugs.map { drug ->
            val d = MolDescriptors.calculate(drug)
            listOf(d.molecularWeight, d.logP, d.hBondDonors.toDouble(),
                   d.hBondAcceptors.toDouble(), d.topologicalPolarSurfaceArea,
                   d.numRotatableBonds.toDouble(), d.numRings.toDouble(),
                   d.molarRefractivity, d.fractionCSP3)
        }
        val keys = listOf("mw","logP","hbd","hba","tpsa","rotB","rings","mr","fsp3")
        val means = keys.indices.map { j -> descriptors.map { it[j] }.average() }
        val stds  = keys.indices.map { j ->
            val m = means[j]
            sqrt(descriptors.map { (it[j] - m).pow(2) }.average())
        }
        val actMean = activities.average()

        // Pearson correlation for each feature as coefficient proxy
        val coeffs = keys.indices.associate { j ->
            val std = if (stds[j] == 0.0) 1.0 else stds[j]
            val cov = activities.indices.map { i -> (descriptors[i][j] - means[j]) * (activities[i] - actMean) }.sum() / drugs.size
            keys[j] to cov / (std * sqrt(activities.map { (it - actMean).pow(2) }.average()))
        }

        val predicted = drugs.mapIndexed { i, drug ->
            coeffs.entries.map { (k, v) ->
                (descriptors[i][keys.indexOf(k)] - means[keys.indexOf(k)]) *
                v / (if (stds[keys.indexOf(k)] == 0.0) 1.0 else stds[keys.indexOf(k)])
            }.sum() + actMean
        }
        val ss_res = activities.indices.map { (activities[it] - predicted[it]).pow(2) }.sum()
        val ss_tot = activities.map { (it - actMean).pow(2) }.sum()
        val r2 = if (ss_tot == 0.0) 0.0 else 1.0 - ss_res / ss_tot

        return QSARModel(coeffs, actMean, r2)
    }
}

// ─── Drug–DNA Binding ────────────────────────────────────────

object DrugDNABinding {

    // ── Intercalation prediction ─────────────────────────────

    fun intercalationScore(drug: Molecule, dna: DNASequence): Double {
        var score = 0.0

        // Planar aromatic system required for intercalation
        val aromaticRings = drug.numAromaticRings()
        val planarScore = when {
            aromaticRings >= 3 -> 1.0
            aromaticRings == 2 -> 0.7
            aromaticRings == 1 -> 0.3
            else -> 0.0
        }
        score += planarScore * 3.0

        // Positive charge favors intercalation (electrostatic with DNA)
        val totalCharge = drug.atoms.map { it.charge.toDouble() }.sum()
        if (totalCharge > 0) score += totalCharge * 0.5

        // Molecular size (too large = clashes)
        val mw = MolDescriptors.molecularWeight(drug)
        val sizePenalty = if (mw > 600) (mw - 600) * 0.01 else 0.0
        score -= sizePenalty

        // GC content affects intercalation preference
        val gcFavor = dna.gcContent() / 100.0 * 0.5
        score += gcFavor

        return score
    }

    // ── Minor groove binding ─────────────────────────────────

    fun minorGrooveScore(drug: Molecule, dna: DNASequence): Double {
        var score = 0.0
        val desc = MolDescriptors.calculate(drug)

        // Crescent shape (low rotatable bonds, compact)
        val shapeScore = maxOf(0.0, 1.0 - desc.numRotatableBonds * 0.1)
        score += shapeScore * 2.0

        // AT preference for minor groove
        val atContent = (100.0 - dna.gcContent()) / 100.0
        score += atContent * 1.5

        // H-bond donors (interact with O2 of AT base pairs)
        score += desc.hBondDonors * 0.5

        // Appropriate width (MW 300-600 preferred)
        val mw = desc.molecularWeight
        score += if (mw in 300.0..600.0) 1.0 else 0.0

        return score
    }

    // ── Covalent binding prediction (alkylating agents) ──────

data class CovalentBindingRisk(
        val isElectrophilic: Boolean,
        val electrophilicAtomCount: Int,
        val structuralAlerts: List<MolDescriptors.StructuralAlert>,
        val riskLevel: String
    )

    fun covalentBindingRisk(drug: Molecule): CovalentBindingRisk {
        val alerts = MolDescriptors.checkStructuralAlerts(drug)
        val highAlerts = alerts.filter { it.severity == "High" }

        val electrophilicAtoms = drug.atoms.filter { atom ->
            // Atoms adjacent to electron-withdrawing groups
            val neighbors = drug.neighborsOf(atom.index)
            atom.symbol == "C" && neighbors.any { n ->
                drug.atoms[n].symbol in setOf("O", "N", "S") &&
                drug.bonds.filter { it.atom1 == n || it.atom2 == n }
                    .any { it.type == BondType.DOUBLE }
            }
        }

        return CovalentBindingRisk(
            isElectrophilic = electrophilicAtoms.isNotEmpty(),
            electrophilicAtomCount = electrophilicAtoms.size,
            structuralAlerts = highAlerts,
            riskLevel = when {
                highAlerts.size >= 2 -> "High"
                highAlerts.isNotEmpty() || electrophilicAtoms.size >= 2 -> "Moderate"
                else -> "Low"
            }
        )
    }

        // ── Full drug–DNA binding profile ────────────────────────

    fun profile(drug: Molecule, dna: DNASequence): DrugDNAProfile {
        val intercalation = intercalationScore(drug, dna)
        val minorGroove = minorGrooveScore(drug, dna)
        val covalentRisk = covalentBindingRisk(drug)
        val desc = MolDescriptors.calculate(drug)

        val bindingMode = when {
            intercalation > 4.0 -> "Intercalation"
            minorGroove > 3.0   -> "Minor groove binding"
            covalentRisk.riskLevel == "High" -> "Covalent alkylation"
            else -> "Electrostatic / non-specific"
        }

        // Approximate ΔG from intercalation score
        val deltaG = -(maxOf(intercalation, minorGroove)) * 0.8

        return DrugDNAProfile(
            drug, dna, bindingMode, intercalation, minorGroove,
            covalentRisk, deltaG, desc
        )
    }

    data class DrugDNAProfile(
        val drug: Molecule,
        val dna: DNASequence,
        val bindingMode: String,
        val intercalationScore: Double,
        val minorGrooveScore: Double,
        val covalentRisk: CovalentBindingRisk,
        val deltaG: Double,
        val descriptors: MolecularDescriptors
    ) {
        override fun toString() = buildString {
            append("Drug-DNA Binding Profile\n")
            append("  Molecule     : ${drug.molecularFormula()} (MW=${"%.2f".format(descriptors.molecularWeight)})\n")
            append("  Target DNA   : ${dna.length} bp, GC=${"%.1f".format(dna.gcContent())}%\n")
            append("  Binding mode : $bindingMode\n")
            append("  ΔG estimate  : ${"%.2f".format(deltaG)} kcal/mol\n")
            append("  Intercalation: ${"%.2f".format(intercalationScore)}\n")
            append("  Minor groove : ${"%.2f".format(minorGrooveScore)}\n")
            append("  Covalent risk: ${covalentRisk.riskLevel}\n")
        }
    }
}

// ─── Drug–RNA Binding ────────────────────────────────────────

object DrugRNABinding {

    // RNA is more structurally diverse than DNA — score based on
    // structure and ligand properties

    fun score(drug: Molecule, rna: RNASequence): Double {
        val structure = RNAInteraction.fold(rna)
        val desc = MolDescriptors.calculate(drug)

        // Structured RNA (more base pairs) creates binding pockets
        val structureScore = structure.numBasePairs.toDouble() / rna.length

        // Positive charge favors RNA binding (electrostatic with phosphate)
        val chargeScore = drug.atoms.count { it.charge > 0 } * 0.5

        // Aromatic rings for intercalation into double-stranded regions
        val aromaticScore = drug.numAromaticRings() * 0.4

        // Aminoglycoside-like features (multiple amines)
        val amineScore = drug.atoms.count { it.symbol == "N" && !it.isAromatic } * 0.3

        return structureScore * 2.0 + chargeScore + aromaticScore + amineScore - desc.logP * 0.2
    }
}

// ─── Metal–Ligand Coordination ───────────────────────────────

object MetalLigandBinding {

    val METAL_SYMBOLS = setOf("Fe","Zn","Cu","Mg","Ca","Mn","Co","Ni","Mo","Se")

    fun hasMetal(mol: Molecule): Boolean = mol.atoms.any { it.symbol in METAL_SYMBOLS }

    fun coordinationScore(mol: Molecule): Double {
        val metals = mol.atoms.filter { it.symbol in METAL_SYMBOLS }
        if (metals.isEmpty()) return 0.0

        var score = 0.0
        for (metal in metals) {
            val ligands = mol.neighborsOf(metal.index).map { mol.atoms[it] }
            // Nitrogen and oxygen are classic ligands
            score += ligands.count { it.symbol in setOf("N", "O", "S") } * 1.5
        }
        return score
    }

    fun metalBindingResidues(protein: ProteinSequence): List<Int> {
        return protein.sequence.indices.filter { i ->
            protein[i] in setOf('H', 'C', 'D', 'E', 'M')
        }
    }
}

// ─── Helper extension ────────────────────────────────────────

