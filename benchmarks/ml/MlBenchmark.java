import biokt.CustomModelLoader;
import biokt.MLPrediction;
import biokt.ProteinMLPredictor;
import biokt.ProteinSequence;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

class MlBenchmark {
    public static void main(String[] args) throws IOException {
        if (args.length != 1) throw new IllegalArgumentException("Pass the prepared input directory");
        Path folder = Path.of(args[0]);
        benchmarkIris(folder);
        benchmarkMembraneClassifier(folder);
    }

    private static void benchmarkIris(Path folder) throws IOException {
        List<String> reference = Files.readAllLines(folder.resolve("iris-reference.tsv"));
        String[] columns = reference.get(0).split("\t");
        List<String> labels = Arrays.asList(columns).subList(9, columns.length);
        var models = new HashMap<String, Object>();
        try (BufferedWriter output = Files.newBufferedWriter(folder.resolve("iris-results.tsv"))) {
            output.write("model\tseed\trow\texpected\treference\tbiokt\tmax_probability_error\t");
            output.write(String.join("\t", labels));
            output.newLine();
            for (String line : reference.subList(1, reference.size())) {
                String[] row = line.split("\t");
                Object model = models.computeIfAbsent(row[0],
                    name -> CustomModelLoader.INSTANCE.load(folder.resolve(name).toString()));
                double[] features = Arrays.stream(row, 5, 9).mapToDouble(Double::parseDouble).toArray();
                MLPrediction prediction = CustomModelLoader.INSTANCE.predict(model, features);
                double maxError = 0;
                for (int i = 0; i < labels.size(); i++) {
                    double expected = Double.parseDouble(row[i + 9]);
                    double observed = prediction.getProbabilities().get(labels.get(i));
                    maxError = Math.max(maxError, Math.abs(expected - observed));
                }
                output.write(String.join("\t", Arrays.copyOf(row, 5)));
                output.write(String.format(Locale.ROOT, "\t%s\t%.16g", prediction.getLabel(), maxError));
                for (String label : labels) {
                    output.write(String.format(Locale.ROOT, "\t%.16g", prediction.getProbabilities().get(label)));
                }
                output.newLine();
            }
        }
    }

    private static void benchmarkMembraneClassifier(Path folder) throws IOException {
        List<String> panel = Files.readAllLines(folder.resolve("protein-panel.tsv"));
        try (BufferedWriter output = Files.newBufferedWriter(folder.resolve("protein-results.tsv"))) {
            output.write("accession\texpected\tbiokt\tprobability_membrane");
            output.newLine();
            for (String line : panel.subList(1, panel.size())) {
                String[] row = line.split("\t");
                var protein = new ProteinSequence(row[2], row[0], "");
                var prediction = ProteinMLPredictor.INSTANCE.predictMembrane(protein);
                output.write(String.format(Locale.ROOT, "%s\t%s\t%s\t%.16g%n",
                    row[0], row[1], prediction.getLabel(), prediction.getProbabilities().get("membrane")));
            }
        }
    }
}
