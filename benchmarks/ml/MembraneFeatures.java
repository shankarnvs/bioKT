import biokt.FeatureExtractor;
import biokt.ProteinMLPredictor;
import biokt.ProteinSequence;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.DoubleStream;

class MembraneFeatures {
    public static void main(String[] args) throws Exception {
        try (var regions = Files.newBufferedWriter(Path.of(args[1]))) {
        regions.write("accession\tstart_1based\tend_1based\n");
        for (String line : Files.readAllLines(Path.of(args[0]))) {
            String[] fields = line.split("\t");
            var protein = new ProteinSequence(fields[2], fields[0], "");
            var extractor = FeatureExtractor.INSTANCE;
            for (var range : extractor.membraneCandidateRegions(protein)) {
                regions.write(fields[0] + "\t" + (range.getFirst() + 1) + "\t" +
                    (range.getLast() + 1) + "\n");
            }
            var features = DoubleStream.concat(Arrays.stream(extractor.membraneFeatures(protein)),
                Arrays.stream(extractor.localMembraneFeatures(protein))).toArray();
            System.out.print(fields[0] + "\t" +
                ProteinMLPredictor.INSTANCE.predictMembrane(protein).getLabel());
            for (double value : features) System.out.print("\t" + value);
            System.out.println();
        }
        }
    }
}
