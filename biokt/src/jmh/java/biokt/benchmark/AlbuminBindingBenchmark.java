package biokt.benchmark;

import biokt.BindingResult;
import biokt.DrugProteinBinding;
import biokt.Molecule;
import biokt.ProteinSequence;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Fork(value = 2, jvmArgsAppend = {"-Xms256m", "-Xmx512m"})
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
public class AlbuminBindingBenchmark {
    private record Ligand(String name, String smiles, double associationConstant) {
        Molecule parse() {
            return new Molecule(smiles, name, name);
        }
    }

    // Association constants from the same human serum albumin study (PMID 7464807).
    private static final List<Ligand> LIGANDS = List.of(
        new Ligand("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 1.3e6),
        new Ligand("naproxen", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1", 1.8e6),
        new Ligand("flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1", 5.0e6)
    );

    @State(Scope.Benchmark)
    public static class Fixture {
        private ProteinSequence protein;
        private List<Molecule> molecules;

        @Setup
        public void prepare() throws IOException {
            try (InputStream stream = AlbuminBindingBenchmark.class.getResourceAsStream("/P02768.fasta")) {
                if (stream == null) {
                    throw new IOException("Missing UniProt albumin fixture");
                }
                String precursor = new String(stream.readAllBytes(), StandardCharsets.US_ASCII)
                    .lines().filter(line -> !line.startsWith(">")).reduce("", String::concat);
                if (precursor.length() != 609) {
                    throw new IOException("Expected the canonical albumin precursor");
                }
                protein = new ProteinSequence(precursor.substring(24), "P02768", "Mature human serum albumin");
            }
            molecules = LIGANDS.stream().map(Ligand::parse).toList();
            for (Molecule molecule : molecules) {
                if (molecule.getAtoms().isEmpty()) {
                    throw new IllegalStateException("Empty molecule: " + molecule.getName());
                }
                checkScore(DrugProteinBinding.INSTANCE.score(molecule, protein));
            }
        }
    }

    @Benchmark
    @OperationsPerInvocation(3)
    public void scoreParsedLigands(Fixture fixture, Blackhole sink) {
        for (Molecule molecule : fixture.molecules) {
            sink.consume(DrugProteinBinding.INSTANCE.score(molecule, fixture.protein));
        }
    }

    @Benchmark
    @OperationsPerInvocation(3)
    public void parseAndScoreLigands(Fixture fixture, Blackhole sink) {
        for (Ligand ligand : LIGANDS) {
            sink.consume(DrugProteinBinding.INSTANCE.score(ligand.parse(), fixture.protein));
        }
    }

    public static void main(String[] args) throws IOException {
        Fixture fixture = new Fixture();
        fixture.prepare();
        System.out.println("ligand\tapparent_kd_um\theuristic_score\tlegacy_kd_proxy_um\t" +
            "held_out_kd_um\tfold_error\tgeometric_mean_baseline_fold_error\textrapolation");
        for (int i = 0; i < LIGANDS.size(); i++) {
            Ligand ligand = LIGANDS.get(i);
            var training = new ArrayList<DrugProteinBinding.MeasuredAffinity>();
            for (int j = 0; j < LIGANDS.size(); j++) {
                if (j != i) {
                    training.add(new DrugProteinBinding.MeasuredAffinity(
                        fixture.molecules.get(j), 1e6 / LIGANDS.get(j).associationConstant()
                    ));
                }
            }
            var scoring = DrugProteinBinding.INSTANCE;
            var calibration = scoring.fitAffinityCalibration(fixture.protein, training);
            var result = scoring.score(fixture.molecules.get(i), fixture.protein, calibration);
            double measured = 1e6 / ligand.associationConstant();
            double predicted = result.getKdMicromolar();
            double baseline = Math.exp(training.stream()
                .mapToDouble(row -> Math.log(row.getKdMicromolar())).average().orElseThrow());
            System.out.printf(Locale.ROOT, "%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f\t%b%n",
                ligand.name(), measured, result.getHeuristic().getBindingScore(),
                result.getHeuristic().getEstimatedKd(), predicted,
                foldError(predicted, measured), foldError(baseline, measured), result.isExtrapolation());
        }
    }

    private static double foldError(double predicted, double observed) {
        return Math.max(predicted / observed, observed / predicted);
    }

    private static void checkScore(BindingResult result) {
        if (!Double.isFinite(result.getBindingScore()) ||
            !Double.isFinite(result.getEstimatedKd()) || result.getEstimatedKd() <= 0) {
            throw new IllegalStateException("Nonfinite score or invalid estimated Kd");
        }
    }
}
