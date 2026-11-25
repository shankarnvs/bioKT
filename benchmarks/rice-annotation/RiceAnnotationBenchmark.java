import biokt.rice.AnnotationOverlap;
import biokt.rice.AnnotationParser;
import biokt.rice.VcfParser;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

class RiceAnnotationBenchmark {
    private record CaseKey(String variantId, String geneId) {}

    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            System.err.println("Usage: RiceAnnotationBenchmark.java <fixture-directory>");
            System.exit(2);
        }

        Path fixtures = Path.of(args[0]);
        var expected = new LinkedHashMap<CaseKey, String>();
        var lines = Files.readAllLines(fixtures.resolve("expected.tsv"));
        if (lines.size() < 2 || !lines.get(0).equals("case_id\tgene_id\texpected")) {
            throw new IllegalArgumentException("Expected a TSV header and at least one benchmark case");
        }
        for (String line : lines.subList(1, lines.size())) {
            String[] columns = line.split("\t", -1);
            if (columns.length != 3 || List.of(columns).stream().anyMatch(String::isBlank)) {
                throw new IllegalArgumentException("Invalid expected result: " + line);
            }
            var key = new CaseKey(columns[0], columns[1]);
            if (expected.putIfAbsent(key, columns[2]) != null) {
                throw new IllegalArgumentException("Duplicate expected case: " + key);
            }
        }

        var annotation = AnnotationParser.INSTANCE.parseGff3(fixtures.resolve("fixture.gff3").toString());
        var variants = VcfParser.INSTANCE.parseString(
            Files.readString(fixtures.resolve("variants.vcf")), null, null, null, false, false, 0, 1
        );
        var actual = new LinkedHashMap<CaseKey, List<String>>();
        for (var result : AnnotationOverlap.INSTANCE.annotateVariants(variants, annotation)) {
            var key = new CaseKey(result.getRecord().getId(), result.getGeneId());
            actual.computeIfAbsent(key, ignored -> new ArrayList<>()).add(result.getConsequence());
        }

        int matches = 0;
        System.out.println("case_id\tgene_id\texpected\tactual\tstatus");
        for (var entry : expected.entrySet()) {
            var key = entry.getKey();
            var labels = actual.remove(key);
            var wanted = entry.getValue();
            boolean match = List.of(wanted).equals(labels);
            if (match) {
                matches++;
            }
            printRow(key, wanted, labels == null ? "." : String.join(",", labels), match ? "PASS" : "FAIL");
        }
        for (var entry : actual.entrySet()) {
            printRow(entry.getKey(), ".", String.join(",", entry.getValue()), "UNEXPECTED");
        }

        int failures = expected.size() - matches + actual.size();
        System.err.printf("%d matched, %d failed or unexpected%n", matches, failures);
        if (failures != 0) {
            System.exit(1);
        }
    }

    private static void printRow(CaseKey key, String expected, String actual, String status) {
        System.out.printf("%s\t%s\t%s\t%s\t%s%n", key.variantId(), key.geneId(), expected, actual, status);
    }
}
