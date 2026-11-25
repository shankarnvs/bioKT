import biokt.DNAInteraction;
import biokt.DNASequence;
import java.util.Locale;

class DuplexMeltingCheck {
    private record Duplex(String id, String sequence, double measuredTm) {}

    public static void main(String[] args) {
        // Rozners and Moulder (2004), Table 1, unmodified DNA only.
        Duplex[] duplexes = {
            new Duplex("1d", "GCGAATTCGC", 47.9),
            new Duplex("2d", "TA".repeat(6), 22.9),
            new Duplex("3a", "CG".repeat(3), 39.2),
            new Duplex("4d", "CGCAAATTTGCG", 56.6),
        };
        double absoluteError = 0;
        double signedError = 0;
        System.out.println("case_id\tdna\tmeasured_tm_c\tbiokt_tm_c\terror_c");
        for (Duplex duplex : duplexes) {
            var strand = new DNASequence(duplex.sequence(), duplex.id(), "");
            if (!strand.getSequence().equals(strand.reverseComplement().getSequence())) {
                throw new AssertionError("Fixture must be self-complementary");
            }
            var prediction = DNAInteraction.INSTANCE.hybridize(strand, strand, 0.300, 2e-6);
            double error = prediction.getTm() - duplex.measuredTm();
            if (!Double.isFinite(error) || prediction.getMismatches() != 0) {
                throw new AssertionError("Unexpected duplex result: " + duplex.id());
            }
            signedError += error;
            absoluteError += Math.abs(error);
            System.out.printf(Locale.ROOT, "%s\t%s\t%.1f\t%.4f\t%.4f%n",
                duplex.id(), duplex.sequence(), duplex.measuredTm(), prediction.getTm(), error);
        }
        System.err.printf(Locale.ROOT, "BioKt MAE %.4f °C, mean signed error %.4f °C%n",
            absoluteError / duplexes.length, signedError / duplexes.length);
    }
}
