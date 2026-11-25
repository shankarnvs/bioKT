import biokt.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

class ProteinMetricValues {
    static void emit(String id, String metric, Collection<? extends Number> values) {
        System.out.println(id + "\t" + metric + "\t" + values.stream().map(Object::toString).collect(Collectors.joining(",")));
    }
    public static void main(String[] args) throws Exception {
        for (var line : Files.readAllLines(Path.of(args[0]))) {
            var row = line.split("\t");
            var id = row[0];
            var seq = new ProteinSequence(row[2], id, "");
            var p = new ProteinAnalysis(seq);
            emit(id,"counts",p.getCounts().values());
            emit(id,"percentages",p.getPercentages().values());
            emit(id,"mass",List.of(p.molecularWeight(false)));
            emit(id,"mono_mass",List.of(p.molecularWeight(true)));
            emit(id,"aromaticity",List.of(p.aromaticity()));
            emit(id,"instability",List.of(p.instabilityIndex()));
            emit(id,"flexibility",p.flexibility());
            emit(id,"charge",List.of(p.chargeAtPH(4),p.chargeAtPH(7),p.chargeAtPH(10)));
            emit(id,"pI",List.of(p.isoelectricPoint()));
            emit(id,"secondary",p.secondaryStructureFractions());
            emit(id,"extinction",p.extinctionCoefficients());
            for (var scale : args[1].split(",")) emit(id,"gravy_"+scale,List.of(p.gravy(scale)));
            for (int w : new int[]{9,19,25}) for(double e : new double[]{0.4,1.0})
                emit(id,"profile_"+w+"_"+e,p.proteinScale("KyteDoolitle",w,e));
            emit(id,"legacy_six",Arrays.stream(FeatureExtractor.INSTANCE.membraneFeatures(seq)).boxed().toList());
        }
    }
}
