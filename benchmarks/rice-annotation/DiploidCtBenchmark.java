import biokt.rice.VcfParser;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

class DiploidCtBenchmark {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException("Usage: DiploidCtBenchmark.java <fixture-directory>");
        }
        Path fixtures = Path.of(args[0]);
        var vcf = VcfParser.INSTANCE.parseString(Files.readString(fixtures.resolve("diploid-ct.vcf")),
            null, null, null, false, false, 0, 1);
        if (vcf.getRecords().size() != 1) {
            throw new AssertionError("Expected one VCF record");
        }
        var site = vcf.getRecords().get(0);
        if (!site.getRef().equals("C") || !site.getAlt().equals(List.of("T"))
            || !site.isSNP() || !site.isTransition() || Math.abs(site.hetFrequency() - 0.6) > 1e-12) {
            throw new AssertionError("Unexpected C-to-T site classification or heterozygosity");
        }

        var expected = Files.readAllLines(fixtures.resolve("diploid-ct-expected.tsv"));
        var samples = vcf.getHeader().getSamples();
        if (expected.size() != samples.size() + 1 || samples.size() != site.getGenotypes().size()) {
            throw new AssertionError("Different sample counts in VCF, parsed data, and expected results");
        }
        var bases = List.of(site.getRef(), site.getAlt().get(0));
        int matches = 0;
        System.out.println(expected.get(0) + "\tstatus");
        for (int i = 0; i < samples.size(); i++) {
            String sample = samples.get(i);
            var genotype = site.getGenotypes().get(sample);
            String decoded = genotype.isMissing() ? genotype.getGt() : genotype.alleles().stream()
                .map(index -> bases.get(Integer.parseInt(index)))
                .collect(Collectors.joining(genotype.isPhased() ? "|" : "/"));
            String actual = String.join("\t", sample, genotype.getGt(), decoded,
                Objects.toString(genotype.dosage(), "NA"), Boolean.toString(genotype.isPhased()),
                Boolean.toString(genotype.isHet()), Boolean.toString(genotype.isMissing()));
            boolean match = actual.equals(expected.get(i + 1));
            System.out.println(actual + "\t" + (match ? "MATCH" : "MISMATCH"));
            if (match) {
                matches++;
            } else {
                System.err.println("Expected: " + expected.get(i + 1));
            }
        }
        System.err.printf("Diploid C/T sample checks: %d/%d matched%n", matches, samples.size());
        if (matches != samples.size()) {
            System.exit(1);
        }
    }
}
