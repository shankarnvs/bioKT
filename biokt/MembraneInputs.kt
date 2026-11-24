package biokt

/** Fixed-width inputs built from the validated ProtParam metrics and local membrane candidates. */
object MembraneInputs {
    fun extract(seq: ProteinSequence): Map<String, Double> = buildMap {
        val p = ProteinAnalysis(seq)
        put("length", seq.length.toDouble())
        p.counts.forEach { (aa, v) -> put("count_$aa", v.toDouble()) }
        p.percentages.forEach { (aa, v) -> put("percent_$aa", v) }
        put("mass", p.molecularWeight())
        put("monoisotopic_mass", p.molecularWeight(true))
        put("instability", p.instabilityIndex())
        put("aromatic_fraction", p.aromaticity())
        put("pi", p.isoelectricPoint())
        put("aliphatic_index", ProteinDescriptors.aliphaticIndex(seq))
        put("hydrophobic_ratio", ProteinDescriptors.hydrophobicRatio(seq))
        put("charged_fraction", "RKHDE".sumOf { p.counts.getValue(it) }.toDouble() / seq.length)
        for (ph in listOf(4, 7, 10)) put("charge_$ph", p.chargeAtPH(ph.toDouble()))
        for ((name, value) in listOf("helix","turn","sheet").zip(p.secondaryStructureFractions())) put(name, value)
        val extinction = p.extinctionCoefficients()
        put("extinction_reduced", extinction[0].toDouble())
        put("extinction_cystines", extinction[1].toDouble())
        for (scale in ProteinAnalysisData.scales.keys.sorted()) put("gravy_$scale", p.gravy(scale))
        fun profile(name: String, values: List<Double>) {
            put("${name}_min", values.minOrNull() ?: 0.0)
            put("${name}_mean", if (values.isEmpty()) 0.0 else values.average())
            put("${name}_max", values.maxOrNull() ?: 0.0)
        }
        profile("flexibility", p.flexibility())
        for (width in listOf(9,19,25)) for (edge in listOf(0.4,1.0)) {
            profile("hydropathy_${width}_$edge", p.proteinScale("KyteDoolitle",width,edge))
        }
        val localNames = listOf("region_mean","region_count","region_coverage","region_charged",
            "flank_charged","region_start","later_region","region_length","region_peak","region_support")
        for ((name,value) in localNames.zip(FeatureExtractor.localMembraneFeatures(seq).asList())) put(name,value)
    }
}
