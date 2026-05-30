package biokt

import kotlin.math.ln

// ─── Phylogenetic Tree ──────────────────────────────────────

sealed class PhyloNode {
    abstract val name: String
    abstract val branchLength: Double

    data class Leaf(
        override val name: String,
        override val branchLength: Double = 0.0
    ) : PhyloNode()

    data class Internal(
        override val name: String = "",
        override val branchLength: Double = 0.0,
        val children: List<PhyloNode> = emptyList()
    ) : PhyloNode()
}

class PhyloTree(val root: PhyloNode) {

    fun toNewick(): String = nodeToNewick(root) + ";"

    private fun nodeToNewick(node: PhyloNode): String {
        return when (node) {
            is PhyloNode.Leaf -> "${node.name}:${formatBranch(node.branchLength)}"
            is PhyloNode.Internal -> {
                val children = node.children.joinToString(",") { nodeToNewick(it) }
                val name = node.name
                val bl = formatBranch(node.branchLength)
                "($children)$name:$bl"
            }
        }
    }

    private fun formatBranch(len: Double) = "%.6f".format(len)

    fun leaves(): List<PhyloNode.Leaf> {
        val result = mutableListOf<PhyloNode.Leaf>()
        collectLeaves(root, result)
        return result
    }

    private fun collectLeaves(node: PhyloNode, acc: MutableList<PhyloNode.Leaf>) {
        when (node) {
            is PhyloNode.Leaf -> acc.add(node)
            is PhyloNode.Internal -> node.children.forEach { collectLeaves(it, acc) }
        }
    }

    fun depth(): Int {
        fun nodeDepth(n: PhyloNode): Int = when (n) {
            is PhyloNode.Leaf -> 1
            is PhyloNode.Internal -> 1 + (n.children.map { nodeDepth(it) }.max() ?: 0)
        }
        return nodeDepth(root)
    }

    fun printAscii(indent: Int = 0, node: PhyloNode = root): String {
        val sb = StringBuilder()
        val prefix = "  ".repeat(indent)
        when (node) {
            is PhyloNode.Leaf     -> sb.append("$prefix+- ${node.name} (${formatBranch(node.branchLength)})\n")
            is PhyloNode.Internal -> {
                val nm = if (node.name.isEmpty()) "internal" else node.name
                sb.append("$prefix+- [$nm] (${formatBranch(node.branchLength)})\n")
                node.children.forEach { sb.append(printAscii(indent + 1, it)) }
            }
        }
        return sb.toString()
    }

    companion object {
        fun fromNewick(newick: String): PhyloTree {
            val clean = newick.trim().trimEnd(';')
            return PhyloTree(parseNewick(clean))
        }

        private fun parseNewick(s: String): PhyloNode {
            if (!s.startsWith("(")) {
                // Leaf
                val parts = s.split(":")
                return PhyloNode.Leaf(parts[0], parts.getOrNull(1)?.toDoubleOrNull() ?: 0.0)
            }
            // Find matching closing paren
            var depth = 0; var closeIdx = -1
            for (i in s.indices) {
                if (s[i] == '(') depth++
            else if (s[i] == ')') { depth--; if (depth == 0) { closeIdx = i } }
            if (closeIdx >= 0 && s[i] == ')' && depth == 0) break
            }
            val inner = s.substring(1, closeIdx)
            val afterParen = s.substring(closeIdx + 1)
            val afterParts = afterParen.split(":")
            val nodeName = afterParts[0]
            val branchLen = afterParts.getOrNull(1)?.toDoubleOrNull() ?: 0.0

            // Split children
            val children = splitTopLevel(inner).map { parseNewick(it.trim()) }
            return PhyloNode.Internal(nodeName, branchLen, children)
        }

        private fun splitTopLevel(s: String): List<String> {
            val parts = mutableListOf<String>()
            var depth = 0; var start = 0
            for (i in s.indices) {
                when (s[i]) {
                    '(' -> depth++; ')' -> depth--
                    ',' -> if (depth == 0) { parts.add(s.substring(start, i)); start = i + 1 }
                }
            }
            parts.add(s.substring(start))
            return parts
        }
    }
}

// ─── Distance-based phylogeny ────────────────────────────────

class DistanceMatrix(
    val taxa: List<String>,
    val distances: Array<DoubleArray>
) {
    val size: Int get() = taxa.size

    operator fun get(i: Int, j: Int): Double = distances[i][j]
    operator fun get(a: String, b: String): Double =
        distances[taxa.indexOf(a)][taxa.indexOf(b)]

    companion object {
        fun fromAlignment(msa: MultipleAlignment): DistanceMatrix {
            val n = msa.numSeqs
            val dist = Array(n) { DoubleArray(n) }
            for (i in 0 until n) for (j in i + 1 until n) {
                val a = msa.alignedSequences[i]; val b = msa.alignedSequences[j]
                val pairs = a.zip(b).filter { (x, y) -> x != '-' && y != '-' }
                val diff = pairs.count { (x, y) -> x != y }
                val d = if (pairs.isEmpty()) 1.0 else diff.toDouble() / pairs.size
                dist[i][j] = d; dist[j][i] = d
            }
            return DistanceMatrix(msa.ids, dist)
        }

        fun jukesCantor(rawDistance: Double): Double {
            if (rawDistance >= 0.75) return Double.MAX_VALUE / 2
            return -0.75 * ln(1.0 - rawDistance * 4.0 / 3.0)
        }
    }

    override fun toString(): String {
        val sb = StringBuilder()
        val nameLen = taxa.map { it.length }.max()!! + 2
        sb.append(" ".repeat(nameLen))
        taxa.forEach { sb.append(it.take(8).padEnd(9)) }
        sb.append("\n")
        for (i in 0 until size) {
            sb.append(taxa[i].padEnd(nameLen))
            for (j in 0 until size) sb.append("%.4f  ".format(distances[i][j]))
            sb.append("\n")
        }
        return sb.toString()
    }
}

// ─── UPGMA Tree Builder ──────────────────────────────────────

object TreeBuilder {

    fun upgma(matrix: DistanceMatrix): PhyloTree {
        data class Cluster(val indices: List<Int>, val node: PhyloNode, val height: Double)

        val dist = Array(matrix.size) { i -> matrix.distances[i].clone() }
        val clusters = MutableList(matrix.size) { i ->
            Cluster(listOf(i), PhyloNode.Leaf(matrix.taxa[i]), 0.0)
        }

        while (clusters.size > 1) {
            // Find minimum distance
            var minD = Double.MAX_VALUE; var minI = 0; var minJ = 1
            for (i in 0 until clusters.size) {
                for (j in i + 1 until clusters.size) {
                    val si = clusters[i].indices; val sj = clusters[j].indices
                    val d = si.flatMap { a -> sj.map { b -> dist[a][b] } }.average()
                    if (d < minD) { minD = d; minI = i; minJ = j }
                }
            }

            val ci = clusters[minI]; val cj = clusters[minJ]
            val newHeight = minD / 2.0
            val childI = adjustBranch(ci.node, ci.height, newHeight)
            val childJ = adjustBranch(cj.node, cj.height, newHeight)
            val newCluster = Cluster(
                ci.indices + cj.indices,
                PhyloNode.Internal(children = listOf(childI, childJ), branchLength = 0.0),
                newHeight
            )

            clusters.removeAt(minJ); clusters.removeAt(minI)
            clusters.add(newCluster)
        }

        return PhyloTree(clusters.first().node)
    }

    private fun adjustBranch(node: PhyloNode, currentHeight: Double, parentHeight: Double): PhyloNode {
        val bl = parentHeight - currentHeight
        return when (node) {
            is PhyloNode.Leaf     -> node.copy(branchLength = bl)
            is PhyloNode.Internal -> node.copy(branchLength = bl)
        }
    }

    // Neighbor Joining
    fun neighborJoining(matrix: DistanceMatrix): PhyloTree {
        val n = matrix.size
        val dist = Array(n) { i -> matrix.distances[i].clone() }
        val labels = matrix.taxa.toMutableList()
        val nodes: MutableList<PhyloNode> = MutableList(n) { PhyloNode.Leaf(matrix.taxa[it]) }

        var activeN = n
        while (activeN > 2) {
            // Compute Q matrix
            val rowSums = DoubleArray(activeN) { i -> (0 until activeN).map { j -> dist[i][j] }.sum() }
            var minQ = Double.MAX_VALUE; var minI = 0; var minJ = 1
            for (i in 0 until activeN) for (j in i + 1 until activeN) {
                val q = (activeN - 2) * dist[i][j] - rowSums[i] - rowSums[j]
                if (q < minQ) { minQ = q; minI = i; minJ = j }
            }

            // Branch lengths
            val dij = dist[minI][minJ]
            val biLen = dij / 2 + (rowSums[minI] - rowSums[minJ]) / (2 * (activeN - 2))
            val bjLen = dij - biLen

            // New node
            val newNode = PhyloNode.Internal(children = listOf(
                adjustBranch(nodes[minI], 0.0, biLen),
                adjustBranch(nodes[minJ], 0.0, bjLen)
            ))

            // Update distances
            val newDist = DoubleArray(activeN) { k ->
                if (k == minI || k == minJ) 0.0
                else (dist[minI][k] + dist[minJ][k] - dij) / 2.0
            }

            // Collapse
            val newDistMatrix = Array(activeN - 1) { DoubleArray(activeN - 1) }
            val keepIdx = (0 until activeN).filter { it != minI && it != minJ }
            for (ai in keepIdx.indices) for (aj in keepIdx.indices)
                newDistMatrix[ai][aj] = dist[keepIdx[ai]][keepIdx[aj]]
            for (ai in keepIdx.indices) {
                newDistMatrix[ai][keepIdx.size] = newDist[keepIdx[ai]]
                newDistMatrix[keepIdx.size][ai] = newDist[keepIdx[ai]]
            }
            for (i in 0 until newDistMatrix.size) dist[i] = newDistMatrix[i]
            nodes.removeAt(maxOf(minI, minJ)); nodes.removeAt(minOf(minI, minJ))
            nodes.add(newNode)
            labels.removeAt(maxOf(minI, minJ)); labels.removeAt(minOf(minI, minJ))
            labels.add("_node_")
            activeN--
        }

        val root = PhyloNode.Internal(children = nodes, branchLength = 0.0)
        return PhyloTree(root)
    }
}
