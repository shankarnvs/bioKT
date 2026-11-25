import biokt.rice.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
class TimingBioKt {
    static BufferedReader gzip(Path path) throws IOException {
        return new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))));
    }
    public static void main(String[] args) throws Exception {
        Path root = Path.of(args[0]);
        int batches = Integer.parseInt(args[1]);
        String ref;
        AnnotationDb db;
        try (var fasta = gzip(root.resolve("rice-chr9.fa.gz"));
             var gff = gzip(root.resolve("rice-exercise-region.gff3.gz"))) {
            if (!fasta.readLine().startsWith(">9 ")) throw new IllegalArgumentException("Wrong FASTA contig");
            ref = fasta.lines().collect(Collectors.joining());
            db = AnnotationParser.INSTANCE.parseGff3String(gff.lines().collect(Collectors.joining("\n")));
        }
        var vcf = new StringBuilder("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (String row : Files.readAllLines(root.resolve("rice_vep.vcf"))) {
            if (!row.isBlank()) vcf.append(row).append("\t.\t.\t.\n");
        }
        var variants = VcfParser.INSTANCE.parseString(vcf.toString(),null,null,null,false,false,0,1);
        if (variants.getRecords().size()!=50) throw new AssertionError("Wrong fixture size");
        var refs = Map.<String,CharSequence>of("9",ref);
        System.out.println(effect(variants,db,refs));
        if (batches > 0) {
            int count = 0;
            for (int i=0; i<1000; i++) count += effect(variants,db,refs).size();
            long start = System.nanoTime();
            for (int i=0; i<batches; i++) count += effect(variants,db,refs).size();
            double millis = (System.nanoTime()-start)/1e6/batches;
            if (count != batches + 1000) throw new AssertionError("Unexpected batch count");
            System.out.println(millis);
        }
    }
    static List<String> effect(VcfFile vcf, AnnotationDb db, Map<String,CharSequence> ref) {
        var rows = new ArrayList<String>();
        for (var hit : AnnotationOverlap.INSTANCE.annotateVariants(vcf,db,ref)) {
            for (var effect : hit.getCodingSnpEffects()) {
                rows.add(hit.getRecord().getPos()+" "+effect.getTranscriptId()+" "+
                    effect.getReferenceCodon()+">"+effect.getAlternateCodon()+" "+
                    effect.getReferenceAminoAcid()+effect.getAminoAcidPosition().toString()+effect.getAlternateAminoAcid());
            }
        }
        return rows;
    }
}