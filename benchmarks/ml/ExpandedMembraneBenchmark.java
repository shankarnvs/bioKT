import biokt.*;
import java.nio.file.Files;
import java.nio.file.Path;

class ExpandedMembraneBenchmark {
    public static void main(String[] args) throws Exception {
        boolean header = true;
        for (String row : Files.readAllLines(Path.of(args[0]))) {
            String[] fields = row.split("\t");
            var sequence = new ProteinSequence(fields[2], fields[0], "");
            if (args.length > 1 && args[1].equals("predict")) {
                var prediction = ProteinMLPredictor.INSTANCE.predictMembrane(sequence);
                System.out.println(fields[0] + "\t" + prediction.getLabel() + "\t" +
                    prediction.getProbabilities().get("membrane"));
            } else {
                var features = MembraneInputs.INSTANCE.extract(sequence);
                if (header) System.out.println("accession\t" + String.join("\t", features.keySet()));
                header = false;
                System.out.print(fields[0]);
                for (double value : features.values()) System.out.print("\t" + value);
                System.out.println();
            }
        }
    }
}
