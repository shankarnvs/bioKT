package biokt.benchmark;

import biokt.DNAInteraction;
import biokt.DNASequence;
import biokt.DrugDNABinding;
import biokt.DrugRNABinding;
import biokt.Fingerprints;
import biokt.InteractionScore;
import biokt.Molecule;
import biokt.ProteinInteraction;
import biokt.ProteinSequence;
import biokt.RNASequence;
import java.util.concurrent.TimeUnit;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Fork(value = 2, jvmArgsAppend = {"-Xms256m", "-Xmx512m"})
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
public class InteractionBenchmarks {
    @State(Scope.Benchmark)
    public static class Fixture {
        private ProteinSequence proteinA;
        private ProteinSequence proteinB;
        private DNASequence dna;
        private DNASequence complement;
        private RNASequence rna;
        private Molecule aspirin;
        private Molecule acetaminophen;
        private boolean[] aspirinFingerprint;
        private boolean[] acetaminophenFingerprint;

        @Setup
        public void prepare() {
            proteinA = new ProteinSequence("ACDEFGHIKLMNPQRSTVWY".repeat(6), "synthetic-a", "");
            proteinB = new ProteinSequence("AVILMFYWSTNQDEKRHGPC".repeat(6), "synthetic-b", "");
            dna = new DNASequence("ACTGATGCTAGTCCAGTACG", "synthetic-dna", "");
            complement = dna.reverseComplement();
            rna = new RNASequence("GCAUCUAGUCGAUAGGCUAC", "synthetic-rna", "");
            aspirin = new Molecule("CC(=O)Oc1ccccc1C(=O)O", "aspirin", "");
            acetaminophen = new Molecule("CC(=O)Nc1ccc(O)cc1", "acetaminophen", "");
            aspirinFingerprint = fingerprint(aspirin);
            acetaminophenFingerprint = fingerprint(acetaminophen);
        }
    }

    @Benchmark
    public InteractionScore proteinProtein(Fixture fixture) {
        return ProteinInteraction.INSTANCE.score(fixture.proteinA, fixture.proteinB);
    }

    @Benchmark
    public DNAInteraction.HybridizationResult dnaDuplex(Fixture fixture) {
        return DNAInteraction.INSTANCE.hybridize(fixture.dna, fixture.complement, 0.05, 250e-9);
    }

    @Benchmark
    public DrugDNABinding.DrugDNAProfile moleculeDna(Fixture fixture) {
        return DrugDNABinding.INSTANCE.profile(fixture.aspirin, fixture.dna);
    }

    @Benchmark
    public double moleculeRna(Fixture fixture) {
        return DrugRNABinding.INSTANCE.score(fixture.aspirin, fixture.rna);
    }

    @Benchmark
    public double moleculeSimilarity(Fixture fixture) {
        return Fingerprints.INSTANCE.tanimoto(fixture.aspirinFingerprint, fixture.acetaminophenFingerprint);
    }

    @Benchmark
    public double fingerprintAndCompare(Fixture fixture) {
        return Fingerprints.INSTANCE.tanimoto(fingerprint(fixture.aspirin), fingerprint(fixture.acetaminophen));
    }

    private static boolean[] fingerprint(Molecule molecule) {
        return Fingerprints.INSTANCE.morgan(molecule, 2, 1024);
    }

    public static void main(String[] args) {
        var fixture = new Fixture();
        fixture.prepare();
        var benchmark = new InteractionBenchmarks();
        double forward = benchmark.proteinProtein(fixture).getScore();
        double reverse = ProteinInteraction.INSTANCE.score(fixture.proteinB, fixture.proteinA).getScore();
        require(Double.isFinite(forward) && Math.abs(forward - reverse) < 1e-12, "Protein score symmetry");
        var exact = benchmark.dnaDuplex(fixture);
        require(exact.getMismatches() == 0 && Double.isFinite(exact.getTm()), "Perfect duplex");
        var moleculeDna = benchmark.moleculeDna(fixture);
        double moleculeRna = benchmark.moleculeRna(fixture);
        require(Double.isFinite(moleculeDna.getDeltaG()) && Double.isFinite(moleculeRna), "Finite ligand scores");

        boolean[] left = fixture.aspirinFingerprint;
        boolean[] right = fixture.acetaminophenFingerprint;
        double actual = benchmark.moleculeSimilarity(fixture);
        require(Fingerprints.INSTANCE.tanimoto(left, left) == 1.0, "Fingerprint identity");
        require(Fingerprints.INSTANCE.tanimoto(right, left) == actual, "Fingerprint symmetry");
        int intersection = 0;
        int union = 0;
        for (int i = 0; i < left.length; i++) {
            if (left[i] && right[i]) intersection++;
            if (left[i] || right[i]) union++;
        }
        require(Math.abs(actual - (double) intersection / union) < 1e-12, "Fingerprint Jaccard calculation");

        String sequence = fixture.complement.getSequence();
        char replacement = sequence.charAt(0) == 'A' ? 'C' : 'A';
        var mismatch = new DNASequence(replacement + sequence.substring(1), "terminal-mismatch", "");
        var terminal = DNAInteraction.INSTANCE.hybridize(fixture.dna, mismatch, 0.05, 250e-9);
        System.out.printf("protein-protein symmetry PASS (score %.6f)%n", forward);
        System.out.printf("perfect duplex PASS (Tm %.3f, mismatches %d)%n", exact.getTm(), exact.getMismatches());
        require(terminal.getMismatches() == 1 && terminal.getTm() < exact.getTm(), "Terminal duplex mismatch");
        System.out.printf("terminal mismatch PASS (expected=1 observed=%d); perfect Tm=%.3f mismatch Tm=%.3f%n",
            terminal.getMismatches(), exact.getTm(), terminal.getTm());
        System.out.printf("finite ligand-DNA and ligand-RNA scores PASS (%.6f, %.6f)%n",
            moleculeDna.getDeltaG(), moleculeRna);
        System.out.printf("fingerprint identity, symmetry and Jaccard PASS (cross-pair %.6f)%n", actual);
    }

    private static void require(boolean condition, String check) {
        if (!condition) throw new AssertionError(check);
    }
}
