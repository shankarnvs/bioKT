package biokt.benchmark;

import biokt.rice.AnnotationDb;
import biokt.rice.AnnotationOverlap;
import biokt.rice.AnnotationParser;
import biokt.rice.VcfFile;
import biokt.rice.VcfParser;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Fork(value = 2, jvmArgsAppend = {"-Xms256m", "-Xmx512m"})
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
public class RiceAnnotationBenchmark {
    @State(Scope.Benchmark)
    public static class Fixture {
        private static final int GENE_COUNT = 128;
        private static final int BASES_PER_GENE = 24;

        @Param({"128", "2048"})
        public int variants;

        private String vcfText;
        private VcfFile vcf;
        private AnnotationDb annotation;
        private Map<String, String> reference;

        @Setup
        public void prepare() {
            StringBuilder gff = new StringBuilder("##gff-version 3\n");
            for (int gene = 0; gene < GENE_COUNT; gene++) {
                int start = gene * BASES_PER_GENE + 1;
                int end = start + BASES_PER_GENE - 1;
                addFeature(gff, "gene", start, end, ".", "ID=g" + gene);
                addFeature(gff, "mRNA", start, end, ".", "ID=t" + gene + ";Parent=g" + gene);
                addFeature(gff, "CDS", start, end, "0", "Parent=t" + gene);
            }
            annotation = AnnotationParser.INSTANCE.parseGff3String(gff.toString());
            String bases = "GCT".repeat(GENE_COUNT * BASES_PER_GENE / 3);
            reference = Map.of("chr1", bases);

            StringBuilder input = new StringBuilder(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            );
            for (int i = 0; i < variants; i++) {
                int position = i * bases.length() / variants;
                input.append("chr1\t").append(position + 1).append("\tv").append(i).append('\t')
                    .append(bases.charAt(position)).append("\tA\t60\tPASS\t.\n");
            }
            vcfText = input.toString();
            vcf = parse(vcfText);

            List<AnnotationOverlap.VariantAnnotation> resolved =
                AnnotationOverlap.INSTANCE.annotateVariants(vcf, annotation, reference);
            if (vcf.getRecords().size() != variants || resolved.size() != variants ||
                resolved.stream().anyMatch(row -> row.getCodingSnpEffects().size() != 1)) {
                throw new IllegalStateException("The benchmark fixture must resolve every coding SNP");
            }
        }

        private static void addFeature(
            StringBuilder gff, String type, int start, int end, String phase, String attributes
        ) {
            gff.append("chr1\tfixture\t").append(type).append('\t').append(start).append('\t')
                .append(end).append("\t.\t+\t").append(phase).append('\t').append(attributes).append('\n');
        }
    }

    @Benchmark
    public VcfFile parseVcf(Fixture fixture) {
        return parse(fixture.vcfText);
    }

    @Benchmark
    public List<AnnotationOverlap.VariantAnnotation> annotateWithoutReference(Fixture fixture) {
        return AnnotationOverlap.INSTANCE.annotateVariants(fixture.vcf, fixture.annotation);
    }

    @Benchmark
    public List<AnnotationOverlap.VariantAnnotation> annotateWithReference(Fixture fixture) {
        return AnnotationOverlap.INSTANCE.annotateVariants(fixture.vcf, fixture.annotation, fixture.reference);
    }

    private static VcfFile parse(String input) {
        return VcfParser.INSTANCE.parseString(input, null, null, null, false, false, 0, 1);
    }
}
