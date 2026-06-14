package biokt

import kotlin.math.*

// ============================================================
// ML PREDICTOR — Inference-only machine learning for BioKt
//
// Architecture:
//   • No training happens here — models are loaded from pretrained
//     weight files (JSON) or from built-in calibrated priors.
//   • Four algorithm families: NaiveBayes, kNN, DecisionTree, LinearClassifier
//   • Three task domains: Sequence, Protein, Drug
//   • Each predictor can be used standalone OR stacked in an Ensemble.
//   • GPTClient (see GPTClient.kt) provides LLM-based interpretation.
// ============================================================

// ─────────────────────────────────────────────────────────────
// CORE DATA TYPES
// ─────────────────────────────────────────────────────────────

/** A prediction result from any classifier */
data class MLPrediction(
    val label: String,
    val confidence: Double,           // 0.0 – 1.0
    val probabilities: Map<String, Double>,
    val features: Map<String, Double>,
    val modelName: String,
    val algorithm: String
) {
    val isHighConfidence: Boolean get() = confidence >= 0.75
    val isMediumConfidence: Boolean get() = confidence in 0.5..0.75
    val isLowConfidence: Boolean get() = confidence < 0.5

    fun summary(): String = buildString {
        append("[$modelName / $algorithm]\n")
        append("  Prediction : $label\n")
        append("  Confidence : ${"%.1f".format(confidence * 100)}%\n")
        append("  Probabilities:\n")
        probabilities.entries.sortedByDescending { it.value }
            .forEach { (lbl, p) -> append("    $lbl : ${"%.3f".format(p)}\n") }
    }
}

/** A collection of predictions from multiple models on the same input */
data class EnsemblePrediction(
    val finalLabel: String,
    val finalConfidence: Double,
    val votes: Map<String, Int>,
    val individualPredictions: List<MLPrediction>
) {
    fun summary(): String = buildString {
        append("=== Ensemble Prediction ===\n")
        append("Final:  $finalLabel (${"%.1f".format(finalConfidence * 100)}% confidence)\n")
        append("Votes:  ${votes.entries.joinToString(", ") { "${it.key}:${it.value}" }}\n")
        append("\nIndividual models:\n")
        individualPredictions.forEach { append(it.summary()) }
    }
}

// ─────────────────────────────────────────────────────────────
// MODEL WEIGHT STORE
// Load pretrained weights from JSON files or use built-in priors
// ─────────────────────────────────────────────────────────────

object ModelStore {

    /** Load weights from a JSON file produced by external training */
    fun loadFromJson(path: String): Map<String, Any> {
        val text = java.io.File(path).readText()
        return parseSimpleJson(text)
    }

    /** Save model weights to JSON */
    fun saveToJson(weights: Map<String, Any>, path: String) {
        java.io.File(path).writeText(toJson(weights))
    }

    private fun toJson(obj: Any?, indent: Int = 0): String {
        val pad = "  ".repeat(indent)
        return when (obj) {
            null          -> "null"
            is String     -> "\"${obj.replace("\"", "\\\"")}\""
            is Double     -> if (obj.isFinite()) obj.toString() else "0.0"
            is Float      -> if (obj.isFinite()) obj.toString() else "0.0"
            is Int        -> obj.toString()
            is Long       -> obj.toString()
            is Boolean    -> obj.toString()
            is Map<*, *>  -> "{\n" + obj.entries.joinToString(",\n") { (k, v) ->
                "$pad  \"$k\": ${toJson(v, indent + 1)}"
            } + "\n$pad}"
            is List<*>    -> "[${obj.joinToString(", ") { toJson(it, indent) }}]"
            is DoubleArray -> "[${obj.joinToString(", ")}]"
            is IntArray    -> "[${obj.joinToString(", ")}]"
            else          -> "\"$obj\""
        }
    }

    /** Very small JSON parser — handles the flat structures BioKt uses */
    fun parseSimpleJson(json: String): Map<String, Any> {
        val result = mutableMapOf<String, Any>()
        val clean  = json.trim().removePrefix("{").removeSuffix("}")
        var i = 0
        while (i < clean.length) {
            // skip whitespace and commas
            while (i < clean.length && (clean[i] == ' ' || clean[i] == '\n' || clean[i] == ',')) i++
            if (i >= clean.length) break
            // key
            if (clean[i] != '"') { i++; continue }
            val keyEnd = clean.indexOf('"', i + 1)
            val key    = clean.substring(i + 1, keyEnd)
            i = keyEnd + 1
            // colon
            while (i < clean.length && clean[i] != ':') i++
            i++
            // skip whitespace
            while (i < clean.length && clean[i] == ' ') i++
            // value
            val (value, nextI) = parseJsonValue(clean, i)
            result[key] = value
            i = nextI
        }
        return result
    }

    private fun parseJsonValue(s: String, start: Int): Pair<Any, Int> {
        var i = start
        return when {
            s[i] == '"' -> {
                val end = s.indexOf('"', i + 1)
                s.substring(i + 1, end) to end + 1
            }
            s[i] == '[' -> {
                val end = s.indexOf(']', i)
                val arr = s.substring(i + 1, end).split(",")
                    .map { it.trim() }
                    .filter { it.isNotEmpty() }
                    .mapNotNull { it.toDoubleOrNull() ?: it.removeSurrounding("\"") }
                arr to end + 1
            }
            s[i] == '{' -> {
                // nested object — find matching }
                var depth = 0; var j = i
                while (j < s.length) { if (s[j]=='{') depth++ else if (s[j]=='}') { depth--; if (depth==0) break }; j++ }
                parseSimpleJson(s.substring(i, j + 1)) to j + 1
            }
            s.startsWith("true",  i) -> true  to i + 4
            s.startsWith("false", i) -> false to i + 5
            s.startsWith("null",  i) -> ""    to i + 4
            else -> {
                var j = i
                while (j < s.length && s[j] != ',' && s[j] != '}' && s[j] != ']' && s[j] != '\n') j++
                val num = s.substring(i, j).trim()
                (num.toDoubleOrNull() ?: num) to j
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
// ALGORITHM IMPLEMENTATIONS (pure Kotlin, no dependencies)
// ─────────────────────────────────────────────────────────────

/** Gaussian Naive Bayes — uses per-class mean/variance per feature */
class NaiveBayesClassifier(
    val name: String,
    private val classPriors: Map<String, Double>,
    private val means: Map<String, DoubleArray>,       // class -> feature means
    private val variances: Map<String, DoubleArray>,   // class -> feature variances
    private val featureNames: List<String>
) {
    fun predict(features: DoubleArray): MLPrediction {
        require(features.size == featureNames.size)
        val logLikelihoods = classPriors.keys.associateWith { cls ->
            var ll = ln(classPriors[cls]!! + 1e-10)
            val mu  = means[cls]!!
            val sig = variances[cls]!!
            for (i in features.indices) {
                val v = sig[i].coerceAtLeast(1e-9)
                ll += -0.5 * ln(2 * PI * v) - (features[i] - mu[i]).pow(2) / (2 * v)
            }
            ll
        }
        val maxLL     = logLikelihoods.values.max()!!
        val expScores = logLikelihoods.mapValues { exp(it.value - maxLL) }
        val total     = expScores.values.sum()
        val probs     = expScores.mapValues { it.value / total }
        val best      = probs.maxBy { it.value }!!
        return MLPrediction(best.key, best.value, probs,
            featureNames.zip(features.toList()).toMap(),
            name, "Naive Bayes")
    }

    companion object {
        /** Build from a weight map (loaded from JSON or built-in) */
        fun fromWeights(name: String, weights: Map<String, Any>): NaiveBayesClassifier {
            @Suppress("UNCHECKED_CAST")
            val priors  = (weights["priors"]  as Map<String, Any>).mapValues { (it.value as Number).toDouble() }
            @Suppress("UNCHECKED_CAST")
            val meansMap= (weights["means"]   as Map<String, Any>).mapValues { (k, v) ->
                ((v as List<Any>).map { (it as Number).toDouble() }).toDoubleArray()
            }
            @Suppress("UNCHECKED_CAST")
            val varMap  = (weights["variances"] as Map<String, Any>).mapValues { (k, v) ->
                ((v as List<Any>).map { (it as Number).toDouble() }).toDoubleArray()
            }
            @Suppress("UNCHECKED_CAST")
            val fNames  = (weights["features"] as List<Any>).map { it.toString() }
            return NaiveBayesClassifier(name, priors, meansMap, varMap, fNames)
        }
    }
}

/** k-Nearest Neighbours — stores training vectors, computes Euclidean distance */
class KNNClassifier(
    val name: String,
    private val k: Int,
    private val trainingPoints: List<Pair<DoubleArray, String>>,
    private val featureNames: List<String>
) {
    fun predict(features: DoubleArray): MLPrediction {
        val distances = trainingPoints.map { (vec, label) ->
            val d = sqrt(vec.indices.map { (vec[it] - features[it]).pow(2) }.sum())
            d to label
        }.sortedBy { it.first }

        val neighbours = distances.take(k)
        val votes = mutableMapOf<String, Int>()
        val weightedVotes = mutableMapOf<String, Double>()
        neighbours.forEach { (dist, lbl) ->
            votes[lbl] = (votes[lbl] ?: 0) + 1
            val w = 1.0 / (dist + 1e-6)
            weightedVotes[lbl] = (weightedVotes[lbl] ?: 0.0) + w
        }

        val total = weightedVotes.values.sum()
        val probs = weightedVotes.mapValues { it.value / total }
        val best  = probs.maxBy { it.value }!!
        return MLPrediction(best.key, best.value, probs,
            featureNames.zip(features.toList()).toMap(),
            name, "kNN(k=$k)")
    }

    companion object {
        fun fromWeights(name: String, weights: Map<String, Any>, k: Int = 5): KNNClassifier {
            @Suppress("UNCHECKED_CAST")
            val points = (weights["points"] as List<Any>).map { pt ->
                val m = pt as Map<String, Any>
                val vec = (m["features"] as List<Any>).map { (it as Number).toDouble() }.toDoubleArray()
                vec to m["label"].toString()
            }
            @Suppress("UNCHECKED_CAST")
            val fNames = (weights["features"] as List<Any>).map { it.toString() }
            return KNNClassifier(name, k, points, fNames)
        }
    }
}

/** Decision Tree node */
data class DTNode(
    val featureIndex: Int = -1,
    val threshold: Double = 0.0,
    val label: String = "",
    val confidence: Double = 1.0,
    val left: DTNode? = null,
    val right: DTNode? = null
) {
    val isLeaf: Boolean get() = left == null && right == null
}

/** Decision Tree Classifier */
class DecisionTreeClassifier(
    val name: String,
    private val root: DTNode,
    private val featureNames: List<String>,
    private val classes: List<String>
) {
    fun predict(features: DoubleArray): MLPrediction {
        var node = root
        val path = mutableListOf<String>()
        while (!node.isLeaf) {
            val v = features[node.featureIndex]
            path.add("${featureNames[node.featureIndex]}${if (v <= node.threshold) "≤" else ">"}${node.threshold}")
            node = if (v <= node.threshold) node.left!! else node.right!!
        }
        val probs = classes.associateWith { if (it == node.label) node.confidence else (1 - node.confidence) / (classes.size - 1).coerceAtLeast(1) }
        return MLPrediction(node.label, node.confidence, probs,
            featureNames.zip(features.toList()).toMap(), name, "Decision Tree")
    }

    companion object {
        fun fromWeights(name: String, weights: Map<String, Any>): DecisionTreeClassifier {
            @Suppress("UNCHECKED_CAST")
            val fNames  = (weights["features"] as List<Any>).map { it.toString() }
            @Suppress("UNCHECKED_CAST")
            val classes = (weights["classes"]  as List<Any>).map { it.toString() }
            @Suppress("UNCHECKED_CAST")
            val tree    = weights["tree"] as Map<String, Any>
            return DecisionTreeClassifier(name, buildNode(tree), fNames, classes)
        }

        @Suppress("UNCHECKED_CAST")
        private fun buildNode(m: Map<String, Any>): DTNode {
            val isLeaf = m["leaf"] as? Boolean ?: false
            return if (isLeaf) {
                DTNode(label = m["label"].toString(),
                    confidence = (m["confidence"] as? Number)?.toDouble() ?: 1.0)
            } else {
                DTNode(
                    featureIndex = (m["feature"] as Number).toInt(),
                    threshold    = (m["threshold"] as Number).toDouble(),
                    left         = buildNode(m["left"]  as Map<String, Any>),
                    right        = buildNode(m["right"] as Map<String, Any>)
                )
            }
        }
    }
}

/** Random Forest — ensemble of Decision Trees with majority vote */
class RandomForestClassifier(
    val name: String,
    private val trees: List<DecisionTreeClassifier>,
    private val classes: List<String>
) {
    fun predict(features: DoubleArray): MLPrediction {
        val votes = mutableMapOf<String, Int>()
        val probAcc = classes.associateWith { 0.0 }.toMutableMap()

        trees.forEach { tree ->
            val pred = tree.predict(features)
            votes[pred.label] = (votes[pred.label] ?: 0) + 1
            pred.probabilities.forEach { (lbl, p) -> probAcc[lbl] = (probAcc[lbl] ?: 0.0) + p }
        }
        val n = trees.size.toDouble()
        val probs = probAcc.mapValues { it.value / n }
        val best  = probs.maxBy { it.value }!!
        return MLPrediction(best.key, best.value, probs,
            mapOf("trees" to n), name, "Random Forest(${trees.size} trees)")
    }

    companion object {
        fun fromWeights(name: String, weights: Map<String, Any>): RandomForestClassifier {
            @Suppress("UNCHECKED_CAST")
            val classes  = (weights["classes"] as List<Any>).map { it.toString() }
            @Suppress("UNCHECKED_CAST")
            val treeList = (weights["trees"] as List<Any>).mapIndexed { i, t ->
                @Suppress("UNCHECKED_CAST")
                val tMap = t as Map<String, Any>
                @Suppress("UNCHECKED_CAST")
                val fNames = (weights["features"] as List<Any>).map { it.toString() }
                @Suppress("UNCHECKED_CAST")
                DecisionTreeClassifier.fromWeights("tree_$i",
                    mapOf("features" to fNames, "classes" to classes, "tree" to tMap))
            }
            return RandomForestClassifier(name, treeList, classes)
        }
    }
}

/** Linear Classifier (Logistic Regression) — w·x + b */
class LinearClassifier(
    val name: String,
    private val weights: Map<String, DoubleArray>,   // class -> weight vector
    private val biases:  Map<String, Double>,
    private val featureNames: List<String>
) {
    fun predict(features: DoubleArray): MLPrediction {
        val scores = weights.keys.associateWith { cls ->
            val w = weights[cls]!!
            val b = biases[cls] ?: 0.0
            w.indices.map { w[it] * features[it] }.sum() + b
        }
        // Softmax
        val maxS   = scores.values.max()!!
        val expScores = scores.mapValues { exp(it.value - maxS) }
        val total  = expScores.values.sum()
        val probs  = expScores.mapValues { it.value / total }
        val best   = probs.maxBy { it.value }!!
        return MLPrediction(best.key, best.value, probs,
            featureNames.zip(features.toList()).toMap(), name, "Logistic Regression")
    }

    companion object {
        fun fromWeights(name: String, weights: Map<String, Any>): LinearClassifier {
            @Suppress("UNCHECKED_CAST")
            val wMap = (weights["weights"] as Map<String, Any>).mapValues { (_, v) ->
                ((v as List<Any>).map { (it as Number).toDouble() }).toDoubleArray()
            }
            @Suppress("UNCHECKED_CAST")
            val bMap = (weights["biases"] as Map<String, Any>).mapValues { (_, v) -> (v as Number).toDouble() }
            @Suppress("UNCHECKED_CAST")
            val fNames = (weights["features"] as List<Any>).map { it.toString() }
            return LinearClassifier(name, wMap, bMap, fNames)
        }
    }
}

// ─────────────────────────────────────────────────────────────
// BUILT-IN PRETRAINED MODELS
// Calibrated priors based on published literature statistics.
// These are NOT deep learning models — they are interpretable
// rule-based classifiers with literature-derived parameters.
// Reference values from PROSITE, UniProtKB, GenBank statistics.
// ─────────────────────────────────────────────────────────────

object PretrainedModels {

    // ── Sequence Models ───────────────────────────────────────

    /** Coding vs non-coding DNA (Fickett 1982 + hexamer approach) */
    val codingSequenceClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "CDS Classifier (Fickett-inspired)",
            classPriors  = mapOf("coding" to 0.55, "non-coding" to 0.45),
            means = mapOf(
                "coding"     to doubleArrayOf(0.62, 0.52, 1.82, 0.72, 0.31, 0.42, 0.88),
                "non-coding" to doubleArrayOf(0.48, 0.38, 1.12, 0.51, 0.52, 0.61, 0.41)
            ),
            variances = mapOf(
                "coding"     to doubleArrayOf(0.04, 0.03, 0.18, 0.06, 0.04, 0.05, 0.08),
                "non-coding" to doubleArrayOf(0.05, 0.04, 0.22, 0.07, 0.05, 0.06, 0.09)
            ),
            featureNames = listOf("gc_content", "at_skew", "linguistic_complexity",
                "codon_bias", "hexamer_score", "repeat_fraction", "orf_density")
        )
    }

    /** Promoter region detector (TATA-box, GC box features) */
    val promoterClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Promoter Classifier",
            classPriors  = mapOf("promoter" to 0.15, "non-promoter" to 0.85),
            means = mapOf(
                "promoter"     to doubleArrayOf(0.58, 0.72, 0.18, 0.65, 0.44),
                "non-promoter" to doubleArrayOf(0.49, 0.41, 0.08, 0.38, 0.22)
            ),
            variances = mapOf(
                "promoter"     to doubleArrayOf(0.04, 0.05, 0.02, 0.06, 0.04),
                "non-promoter" to doubleArrayOf(0.03, 0.04, 0.01, 0.05, 0.03)
            ),
            featureNames = listOf("gc_content", "cpg_density", "tata_box_score",
                "gc_box_score", "caat_box_score")
        )
    }

    /** Splice site detector (donor/acceptor/non-splice) */
    val spliceSiteClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Splice Site Classifier",
            classPriors  = mapOf("donor" to 0.12, "acceptor" to 0.12, "non-splice" to 0.76),
            means = mapOf(
                "donor"      to doubleArrayOf(0.82, 0.91, 0.28, 0.73, 0.65),
                "acceptor"   to doubleArrayOf(0.71, 0.38, 0.84, 0.29, 0.71),
                "non-splice" to doubleArrayOf(0.49, 0.44, 0.45, 0.42, 0.38)
            ),
            variances = mapOf(
                "donor"      to doubleArrayOf(0.03, 0.02, 0.04, 0.05, 0.04),
                "acceptor"   to doubleArrayOf(0.04, 0.03, 0.03, 0.04, 0.03),
                "non-splice" to doubleArrayOf(0.06, 0.05, 0.06, 0.06, 0.05)
            ),
            featureNames = listOf("gc_local", "gt_ag_score", "py_tract_score",
                "branch_point_score", "conservation_score")
        )
    }

    // ── Protein Models ────────────────────────────────────────

    /** Enzyme vs non-enzyme (uses physicochemical descriptors) */
    val enzymeClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Enzyme Classifier",
            classPriors  = mapOf("enzyme" to 0.42, "non-enzyme" to 0.58),
            means = mapOf(
                "enzyme"     to doubleArrayOf(-0.12, 0.38, 94.2, 5.8, 0.44, 7.12, 1480.0, 0.22),
                "non-enzyme" to doubleArrayOf( 0.08, 0.29, 78.1, 6.4, 0.51, 6.82, 1820.0, 0.31)
            ),
            variances = mapOf(
                "enzyme"     to doubleArrayOf(0.18, 0.06, 82.0, 0.9, 0.05, 1.2, 280000.0, 0.04),
                "non-enzyme" to doubleArrayOf(0.22, 0.07, 91.0, 1.1, 0.06, 1.4, 310000.0, 0.05)
            ),
            featureNames = listOf("gravy", "hydrophobic_ratio", "aliphatic_index",
                "pi", "aromatic_fraction", "charge_ph7", "molecular_weight", "instability_idx_norm")
        )
    }

    /** Membrane vs soluble protein */
    val membraneClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Membrane Protein Classifier",
            classPriors  = mapOf("membrane" to 0.27, "soluble" to 0.73),
            means = mapOf(
                "membrane" to doubleArrayOf(1.12, 0.62, 0.14, 102.1, 0.08, 3.9),
                "soluble"  to doubleArrayOf(-0.32, 0.38, 0.38, 79.4,  0.31, 7.1)
            ),
            variances = mapOf(
                "membrane" to doubleArrayOf(0.22, 0.05, 0.03, 68.0, 0.02, 1.8),
                "soluble"  to doubleArrayOf(0.28, 0.06, 0.05, 72.0, 0.04, 2.1)
            ),
            featureNames = listOf("gravy", "hydrophobic_ratio", "charged_fraction",
                "aliphatic_index", "aromatic_fraction", "pi")
        )
    }

    // ── Drug / Molecule Models ────────────────────────────────

    /** Drug active vs inactive (binary activity classifier) */
    val drugActivityClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Drug Activity Classifier",
            classPriors  = mapOf("active" to 0.35, "inactive" to 0.65),
            means = mapOf(
                "active"   to doubleArrayOf(320.0, 2.8, 2.0, 5.5, 62.0, 4.0, 1.8, 0.52),
                "inactive" to doubleArrayOf(420.0, 3.9, 3.1, 7.2, 88.0, 5.8, 2.9, 0.38)
            ),
            variances = mapOf(
                "active"   to doubleArrayOf(8100.0, 1.8, 1.2, 4.1, 680.0, 2.4, 0.9, 0.04),
                "inactive" to doubleArrayOf(9200.0, 2.1, 1.4, 5.2, 820.0, 3.1, 1.2, 0.05)
            ),
            featureNames = listOf("molecular_weight", "logP", "hbd", "hba",
                "tpsa", "rot_bonds", "rings", "fsp3")
        )
    }

    /** Toxicity class I/II/III (Cramer rules enhanced with Lipinski) */
    val toxicityClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "Toxicity Classifier (Cramer-enhanced)",
            classPriors  = mapOf("class_I_low" to 0.48, "class_II_moderate" to 0.32, "class_III_high" to 0.20),
            means = mapOf(
                "class_I_low"      to doubleArrayOf(180.0, 0.8, 1.1, 3.2, 42.0, 2.1, 0.8),
                "class_II_moderate" to doubleArrayOf(310.0, 2.4, 2.3, 5.1, 72.0, 4.2, 1.6),
                "class_III_high"   to doubleArrayOf(450.0, 4.1, 3.8, 7.8, 98.0, 6.1, 2.4)
            ),
            variances = mapOf(
                "class_I_low"      to doubleArrayOf(3200.0, 0.6, 0.8, 2.1, 320.0, 1.4, 0.4),
                "class_II_moderate" to doubleArrayOf(6800.0, 1.2, 1.1, 3.4, 580.0, 2.2, 0.8),
                "class_III_high"   to doubleArrayOf(9400.0, 2.2, 1.8, 4.8, 920.0, 3.6, 1.4)
            ),
            featureNames = listOf("molecular_weight", "logP", "hbd", "hba",
                "tpsa", "rot_bonds", "aromatic_rings")
        )
    }

    /** hERG cardiotoxicity predictor */
    val hergClassifier: NaiveBayesClassifier by lazy {
        NaiveBayesClassifier(
            name = "hERG Inhibition Classifier",
            classPriors  = mapOf("inhibitor" to 0.23, "non-inhibitor" to 0.77),
            means = mapOf(
                "inhibitor"     to doubleArrayOf(3.8, 350.0, 1.2, 3.8, 2.1, 42.0),
                "non-inhibitor" to doubleArrayOf(1.4, 280.0, 1.9, 5.8, 0.8, 78.0)
            ),
            variances = mapOf(
                "inhibitor"     to doubleArrayOf(1.4, 7200.0, 0.9, 2.8, 1.6, 420.0),
                "non-inhibitor" to doubleArrayOf(1.8, 6100.0, 1.2, 3.4, 0.8, 580.0)
            ),
            featureNames = listOf("logP", "molecular_weight", "hbd",
                "hba", "aromatic_rings", "tpsa")
        )
    }
}

// ─────────────────────────────────────────────────────────────
// FEATURE EXTRACTORS
// Convert BioKt objects to feature vectors for each model
// ─────────────────────────────────────────────────────────────

object FeatureExtractor {

    // ── DNA / Sequence features ────────────────────────────────

    fun codingSequenceFeatures(seq: DNASequence): DoubleArray {
        val gc       = seq.gcContent() / 100.0
        val freqs    = SeqStats.nucleotideFrequency(seq)
        val a        = freqs['A'] ?: 0.25
        val t        = freqs['T'] ?: 0.25
        val atSkew   = if (a + t > 0) (a - t) / (a + t) else 0.0
        val complexity = SeqStats.linguisticComplexity(seq)
        val orfs     = seq.findOrfs(minLength = 30)
        val orfDensity = orfs.size.toDouble() / (seq.length / 100.0).coerceAtLeast(1.0)

        // Codon bias proxy — count in-frame stop codons (lower = more likely coding)
        val triplets = (0 until seq.length - 2 step 3).count { i ->
            val codon = seq.sequence.substring(i, i + 3)
            CodonTables.STANDARD.translate(codon) == '*'
        }
        val codonBias = 1.0 - (triplets.toDouble() / ((seq.length / 3).coerceAtLeast(1)))

        // Hexamer score proxy — GC in first/third codon positions
        val hex = (0 until seq.length - 2 step 3).map { i ->
            val c = seq.sequence[i + 2]
            if (c == 'G' || c == 'C') 1.0 else 0.0
        }.average().let { if (it.isNaN()) 0.5 else it }

        // Repeat fraction — dinucleotide entropy proxy
        val repeatFrac = 1.0 - complexity

        return doubleArrayOf(gc, atSkew, complexity, codonBias, hex, repeatFrac, orfDensity)
    }

    fun promoterFeatures(seq: DNASequence): DoubleArray {
        val s = seq.sequence
        val gc = seq.gcContent() / 100.0

        // CpG density (C followed by G)
        val cpgCount = (0 until s.length - 1).count { i -> s[i] == 'C' && s[i+1] == 'G' }
        val cpgDensity = cpgCount.toDouble() / (s.length - 1).coerceAtLeast(1)

        // TATA box score — look for TATAAA or similar
        val tataBoxes = listOf("TATAAAA","TATAAAT","TATAAAG","TATAATG")
        val tataScore = if (tataBoxes.any { s.contains(it) }) 1.0
                        else if (s.contains("TATAA")) 0.5
                        else if (s.contains("TATA")) 0.2 else 0.0

        // GC box: GGGCGG
        val gcBoxScore = when {
            s.contains("GGGCGG") -> 1.0
            s.contains("GGCGG")  -> 0.6
            gc > 0.65            -> 0.3
            else                 -> 0.0
        }

        // CAAT box: CCAAT
        val caatScore = if (s.contains("CCAAT") || s.contains("CAAT")) 0.8 else 0.0

        return doubleArrayOf(gc, cpgDensity, tataScore, gcBoxScore, caatScore)
    }

    fun spliceSiteFeatures(seq: DNASequence, position: Int = 0): DoubleArray {
        val s = seq.sequence
        val window = 20
        val start  = maxOf(0, position - window / 2)
        val end    = minOf(s.length, position + window / 2)
        val local  = if (start < end) s.substring(start, end) else s

        val gcLocal  = local.count { it == 'G' || it == 'C' }.toDouble() / local.length.coerceAtLeast(1)

        // GT-AG rule score
        val gtScore = if (position + 1 < s.length && s[position] == 'G' && s[position+1] == 'T') 1.0 else 0.0
        val agScore = if (position + 1 < s.length && s[position] == 'A' && s[position+1] == 'G') 1.0 else 0.0
        val gtAgScore = maxOf(gtScore, agScore)

        // Polypyrimidine tract score (CT-rich region before acceptor)
        val pyScore = local.count { it == 'C' || it == 'T' }.toDouble() / local.length.coerceAtLeast(1)

        // Branch point proxy (YNYURAY pattern — simplified)
        val branchScore = if (local.contains("TACTAAC") || local.contains("TACTAAG")) 0.9
                          else if (local.count { it == 'A' } / local.length.coerceAtLeast(1).toDouble() > 0.35) 0.4
                          else 0.1

        val conservScore = gtAgScore * 0.5 + gcLocal * 0.3 + pyScore * 0.2

        return doubleArrayOf(gcLocal, gtAgScore, pyScore, branchScore, conservScore)
    }

    // ── Protein features ──────────────────────────────────────

    fun proteinFeatures(seq: ProteinSequence): DoubleArray = doubleArrayOf(
        ProteinDescriptors.gravyIndex(seq),
        ProteinDescriptors.hydrophobicRatio(seq),
        ProteinDescriptors.aliphaticIndex(seq),
        seq.isoelectricPoint(),
        seq.aromaticity(),
        ProteinDescriptors.chargeAtPH(seq, 7.0),
        seq.molecularWeight(),
        (seq.instabilityIndex() / 100.0).coerceIn(0.0, 1.0)
    )

    fun membraneFeatures(seq: ProteinSequence): DoubleArray {
        val charged = listOf('R','K','H','D','E').map { seq.count(it) }.sum().toDouble() / seq.length
        return doubleArrayOf(
            ProteinDescriptors.gravyIndex(seq),
            ProteinDescriptors.hydrophobicRatio(seq),
            charged,
            ProteinDescriptors.aliphaticIndex(seq),
            seq.aromaticity(),
            seq.isoelectricPoint()
        )
    }

    // ── Drug / Molecule features ──────────────────────────────

    fun drugFeatures(mol: Molecule): DoubleArray {
        val d = MolDescriptors.calculate(mol)
        return doubleArrayOf(
            d.molecularWeight,
            d.logP,
            d.hBondDonors.toDouble(),
            d.hBondAcceptors.toDouble(),
            d.topologicalPolarSurfaceArea,
            d.numRotatableBonds.toDouble(),
            d.numRings.toDouble(),
            d.fractionCSP3
        )
    }

    fun toxicityFeatures(mol: Molecule): DoubleArray {
        val d = MolDescriptors.calculate(mol)
        return doubleArrayOf(
            d.molecularWeight,
            d.logP,
            d.hBondDonors.toDouble(),
            d.hBondAcceptors.toDouble(),
            d.topologicalPolarSurfaceArea,
            d.numRotatableBonds.toDouble(),
            d.numAromaticRings.toDouble()
        )
    }

    fun hergFeatures(mol: Molecule): DoubleArray {
        val d = MolDescriptors.calculate(mol)
        return doubleArrayOf(
            d.logP,
            d.molecularWeight,
            d.hBondDonors.toDouble(),
            d.hBondAcceptors.toDouble(),
            d.numAromaticRings.toDouble(),
            d.topologicalPolarSurfaceArea
        )
    }
}

// ─────────────────────────────────────────────────────────────
// HIGH-LEVEL PREDICTOR APIs
// Simple one-call interfaces for each task domain
// ─────────────────────────────────────────────────────────────

object SequenceMLPredictor {

    /** Is this sequence likely to be a coding region (CDS)? */
    fun predictCoding(seq: DNASequence): MLPrediction =
        PretrainedModels.codingSequenceClassifier
            .predict(FeatureExtractor.codingSequenceFeatures(seq))

    /** Is this sequence likely to be a promoter? */
    fun predictPromoter(seq: DNASequence): MLPrediction =
        PretrainedModels.promoterClassifier
            .predict(FeatureExtractor.promoterFeatures(seq))

    /** Is this position likely to be a splice site? */
    fun predictSpliceSite(seq: DNASequence, position: Int = 0): MLPrediction =
        PretrainedModels.spliceSiteClassifier
            .predict(FeatureExtractor.spliceSiteFeatures(seq, position))

    /** Run all three classifiers and return an ensemble */
    fun fullAnalysis(seq: DNASequence): Map<String, MLPrediction> = mapOf(
        "coding"      to predictCoding(seq),
        "promoter"    to predictPromoter(seq),
        "splice_site" to predictSpliceSite(seq)
    )
}

object ProteinMLPredictor {

    /** Is this protein likely to be an enzyme? */
    fun predictEnzyme(seq: ProteinSequence): MLPrediction =
        PretrainedModels.enzymeClassifier
            .predict(FeatureExtractor.proteinFeatures(seq))

    /** Is this protein membrane-associated or soluble? */
    fun predictMembrane(seq: ProteinSequence): MLPrediction =
        PretrainedModels.membraneClassifier
            .predict(FeatureExtractor.membraneFeatures(seq))

    /** Run both classifiers */
    fun fullAnalysis(seq: ProteinSequence): Map<String, MLPrediction> = mapOf(
        "enzyme_prediction"   to predictEnzyme(seq),
        "membrane_prediction" to predictMembrane(seq)
    )
}

object DrugMLPredictor {

    /** Predict drug activity class (active / inactive) */
    fun predictActivity(mol: Molecule): MLPrediction =
        PretrainedModels.drugActivityClassifier
            .predict(FeatureExtractor.drugFeatures(mol))

    /** Predict toxicity class (I / II / III) */
    fun predictToxicity(mol: Molecule): MLPrediction =
        PretrainedModels.toxicityClassifier
            .predict(FeatureExtractor.toxicityFeatures(mol))

    /** Predict hERG cardiac risk */
    fun predictHerg(mol: Molecule): MLPrediction =
        PretrainedModels.hergClassifier
            .predict(FeatureExtractor.hergFeatures(mol))

    /** Run all three classifiers */
    fun fullAnalysis(mol: Molecule): Map<String, MLPrediction> = mapOf(
        "activity"  to predictActivity(mol),
        "toxicity"  to predictToxicity(mol),
        "herg_risk" to predictHerg(mol)
    )
}

// ─────────────────────────────────────────────────────────────
// CUSTOM MODEL LOADER
// Load any pretrained model from a JSON weight file
// ─────────────────────────────────────────────────────────────

object CustomModelLoader {

    /**
     * Load a pretrained model from a JSON file.
     *
     * The JSON must contain a "algorithm" field:
     *   "naive_bayes", "knn", "decision_tree",
     *   "random_forest", "linear"
     *
     * Example JSON structure:
     * {
     *   "name": "My Classifier",
     *   "algorithm": "naive_bayes",
     *   "features": ["feat1", "feat2"],
     *   "priors": {"classA": 0.5, "classB": 0.5},
     *   "means": {"classA": [0.1, 0.2], "classB": [0.3, 0.4]},
     *   "variances": {"classA": [0.01, 0.01], "classB": [0.01, 0.01]}
     * }
     */
    fun load(jsonPath: String): Any {
        val weights = ModelStore.loadFromJson(jsonPath)
        val algo    = weights["algorithm"]?.toString() ?: "naive_bayes"
        val name    = weights["name"]?.toString() ?: jsonPath
        return when (algo.toLowerCase()) {
            "naive_bayes", "nb"  -> NaiveBayesClassifier.fromWeights(name, weights)
            "knn"                -> KNNClassifier.fromWeights(name, weights)
            "decision_tree", "dt"-> DecisionTreeClassifier.fromWeights(name, weights)
            "random_forest", "rf"-> RandomForestClassifier.fromWeights(name, weights)
            "linear", "lr"       -> LinearClassifier.fromWeights(name, weights)
            else -> throw IllegalArgumentException("Unknown algorithm: $algo")
        }
    }

    /** Predict using a loaded model (pass the result of load()) */
    fun predict(model: Any, features: DoubleArray): MLPrediction = when (model) {
        is NaiveBayesClassifier  -> model.predict(features)
        is KNNClassifier         -> model.predict(features)
        is DecisionTreeClassifier-> model.predict(features)
        is RandomForestClassifier-> model.predict(features)
        is LinearClassifier      -> model.predict(features)
        else -> throw IllegalArgumentException("Unknown model type")
    }
}

// ─────────────────────────────────────────────────────────────
// ENSEMBLE
// Combine multiple models with majority vote
// ─────────────────────────────────────────────────────────────

object Ensemble {

    fun vote(predictions: List<MLPrediction>): EnsemblePrediction {
        val votes   = mutableMapOf<String, Int>()
        val confAcc = mutableMapOf<String, Double>()
        predictions.forEach { p ->
            votes[p.label]   = (votes[p.label] ?: 0) + 1
            confAcc[p.label] = (confAcc[p.label] ?: 0.0) + p.confidence
        }
        val best  = votes.maxBy { it.value }!!.key
        val avgConf = confAcc[best]!! / votes[best]!!
        return EnsemblePrediction(best, avgConf, votes, predictions)
    }

    fun weightedVote(predictions: List<Pair<MLPrediction, Double>>): EnsemblePrediction {
        val scoreMap = mutableMapOf<String, Double>()
        predictions.forEach { (pred, weight) ->
            pred.probabilities.forEach { (lbl, prob) ->
                scoreMap[lbl] = (scoreMap[lbl] ?: 0.0) + prob * weight
            }
        }
        val totalWeight = predictions.map { it.second }.sum()
        val probs = scoreMap.mapValues { it.value / totalWeight }
        val best  = probs.maxBy { it.value }!!
        val votes = predictions.groupBy { it.first.label }.mapValues { it.value.size }
        return EnsemblePrediction(best.key, best.value, votes, predictions.map { it.first })
    }
}
