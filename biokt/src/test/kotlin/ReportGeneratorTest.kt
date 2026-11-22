package biokt

import java.nio.file.Path
import kotlin.test.Test
import kotlin.test.assertContains
import org.junit.jupiter.api.io.TempDir

class ReportGeneratorTest {
    @Test
    fun writesReportsAndTheirViewers(@TempDir directory: Path) {
        val dnaReport = directory.resolve("dna.html").toFile()
        val proteinReport = directory.resolve("protein.html").toFile()
        val drugReport = directory.resolve("drug.html").toFile()

        ReportGenerator.dnaReport(DNASequence("ATGGCTGAATGA"), dnaReport.path)
        ReportGenerator.proteinReport(ProteinSequence("ACDEFGHIKLMNPQRSTVWY"), proteinReport.path)
        ReportGenerator.drugReport(Molecule("CC(=O)Oc1ccccc1C(=O)O", name = "Aspirin"), drugReport.path)

        assertContains(dnaReport.readText(), "dna_viewer.html?seq=")
        assertContains(proteinReport.readText(), "protein_viewer.html?seq=")
        assertContains(drugReport.readText(), "molecule_viewer.html?smiles=")
        for (viewer in listOf("dna_viewer.html", "protein_viewer.html", "molecule_viewer.html")) {
            assertContains(directory.resolve(viewer).toFile().readText(), "<html", ignoreCase = true)
        }
    }
}
