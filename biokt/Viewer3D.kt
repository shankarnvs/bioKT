package biokt

object Viewer3D {
    fun buildDNAViewer(): String = load("dna_viewer.html")
    fun buildProteinViewer(): String = load("protein_viewer.html")
    fun buildMoleculeViewer(): String = load("molecule_viewer.html")

    private fun load(filename: String): String =
        checkNotNull(Viewer3D::class.java.getResourceAsStream("/$filename")) {
            "Missing bundled viewer: $filename"
        }.bufferedReader(Charsets.UTF_8).use { it.readText() }
}
