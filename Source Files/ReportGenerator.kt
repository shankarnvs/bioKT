package biokt

import kotlin.math.*

// ============================================================
// REPORT GENERATOR
// Produces self-contained HTML reports for:
//   - DNA Analysis
//   - Protein Analysis
//   - Drug / Molecule Analysis
// Each report has a 1-page summary + detailed appendix.
// ============================================================

object ReportGenerator {

    // ── Entry points ─────────────────────────────────────────

    fun dnaReport(seq: DNASequence, outputPath: String = "dna_report.html") {
        val dir = java.io.File(outputPath).parent ?: "."
        val viewerPath = "$dir/dna_viewer.html"
        java.io.File(viewerPath).writeText(Viewer3D.buildDNAViewer())
        val html = buildDNAReport(seq, viewerPath = viewerPath)
        java.io.File(outputPath).writeText(html)
        println("DNA report written: $outputPath")
        println("DNA 3D viewer written: $viewerPath")
    }

    fun proteinReport(seq: ProteinSequence, outputPath: String = "protein_report.html") {
        val dir = java.io.File(outputPath).parent ?: "."
        val viewerPath = "$dir/protein_viewer.html"
        java.io.File(viewerPath).writeText(Viewer3D.buildProteinViewer())
        val html = buildProteinReport(seq, viewerPath = viewerPath)
        java.io.File(outputPath).writeText(html)
        println("Protein report written: $outputPath")
        println("Protein 3D viewer written: $viewerPath")
    }

    fun drugReport(mol: Molecule, outputPath: String = "drug_report.html") {
        val dir = java.io.File(outputPath).parent ?: "."
        val viewerPath = "$dir/molecule_viewer.html"
        java.io.File(viewerPath).writeText(Viewer3D.buildMoleculeViewer())
        val html = buildDrugReport(mol, viewerPath = viewerPath)
        java.io.File(outputPath).writeText(html)
        println("Drug report written: $outputPath")
        println("Molecule 3D viewer written: $viewerPath")
    }

    // ── Shared CSS + JS ──────────────────────────────────────

    private fun sharedStyles() = """
<style>
:root {
  --bg: #f8fafc;
  --surface: #ffffff;
  --border: #e2e8f0;
  --blue: #1e40af;
  --teal: #0f766e;
  --green: #15803d;
  --orange: #c2410c;
  --purple: #6d28d9;
  --red: #dc2626;
  --gray: #374151;
  --lgray: #6b7280;
  --xlgray: #f1f5f9;
  --text: #111827;
  --pass: #dcfce7;
  --pass-text: #15803d;
  --fail: #fee2e2;
  --fail-text: #dc2626;
  --warn: #fef3c7;
  --warn-text: #92400e;
}
* { margin:0; padding:0; box-sizing:border-box; }
body { font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  background:var(--bg); color:var(--text); font-size:14px; line-height:1.6; }
a { color:var(--blue); }

/* ── Layout ── */
.page { max-width:1100px; margin:0 auto; padding:32px 24px 60px; }

/* ── Report header ── */
.report-header {
  background:linear-gradient(135deg,#1e3a8a 0%,#1e40af 50%,#0f766e 100%);
  color:#fff; border-radius:14px; padding:32px 36px; margin-bottom:32px;
  display:flex; justify-content:space-between; align-items:flex-start;
  flex-wrap:wrap; gap:16px;
}
.report-header h1 { font-size:28px; font-weight:800; margin-bottom:6px; }
.report-header .subtitle { font-size:14px; opacity:.85; }
.report-header .meta { text-align:right; font-size:12px; opacity:.8; line-height:1.8; }
.report-type-badge {
  display:inline-block; background:rgba(255,255,255,.2); border:1px solid rgba(255,255,255,.3);
  border-radius:20px; padding:3px 12px; font-size:12px; font-weight:700;
  text-transform:uppercase; letter-spacing:.5px; margin-bottom:8px;
}

/* ── Section ── */
.section { margin-bottom:32px; }
.section-title {
  font-size:16px; font-weight:700; color:var(--blue);
  border-bottom:2px solid var(--border); padding-bottom:8px; margin-bottom:16px;
  display:flex; align-items:center; gap:8px;
}
.section-title .icon { font-size:18px; }

/* ── Summary grid ── */
.summary-grid {
  display:grid; grid-template-columns:repeat(auto-fill,minmax(160px,1fr));
  gap:12px; margin-bottom:24px;
}
.stat-card {
  background:var(--surface); border:1px solid var(--border); border-radius:10px;
  padding:14px 16px; text-align:center;
}
.stat-card .val { font-size:24px; font-weight:800; color:var(--blue); line-height:1.2; }
.stat-card .lbl { font-size:11px; color:var(--lgray); margin-top:4px; text-transform:uppercase; letter-spacing:.4px; }
.stat-card.green .val  { color:var(--green); }
.stat-card.orange .val { color:var(--orange); }
.stat-card.purple .val { color:var(--purple); }
.stat-card.teal .val   { color:var(--teal); }
.stat-card.red .val    { color:var(--red); }

/* ── Tables ── */
.data-table { width:100%; border-collapse:collapse; border:1px solid var(--border);
  border-radius:8px; overflow:hidden; }
.data-table thead tr { background:var(--blue); }
.data-table thead th { padding:9px 14px; text-align:left; font-size:12px;
  font-weight:700; color:#fff; text-transform:uppercase; letter-spacing:.4px; }
.data-table tbody tr:nth-child(even) { background:var(--xlgray); }
.data-table tbody tr:hover { background:#e0f2fe; }
.data-table td { padding:8px 14px; font-size:13px; border-bottom:1px solid var(--border); }
.data-table td.mono { font-family:'Cascadia Code','Fira Code',monospace; font-size:12px; }
.data-table td.num { text-align:right; font-variant-numeric:tabular-nums; }

/* ── Badge chips ── */
.badge {
  display:inline-block; padding:2px 9px; border-radius:12px; font-size:11px;
  font-weight:700; text-transform:uppercase; letter-spacing:.3px;
}
.badge.pass  { background:var(--pass);  color:var(--pass-text); }
.badge.fail  { background:var(--fail);  color:var(--fail-text); }
.badge.warn  { background:var(--warn);  color:var(--warn-text); }
.badge.info  { background:#dbeafe; color:#1e40af; }
.badge.purple{ background:#ede9fe; color:#6d28d9; }

/* ── Sequence display ── */
.seq-display {
  font-family:'Cascadia Code','Fira Code',monospace; font-size:12px;
  background:#0f172a; color:#e2e8f0; border-radius:8px; padding:16px 18px;
  line-height:2; word-break:break-all; overflow-x:auto;
}
.seq-display .pos { color:#64748b; font-size:10px; margin-right:8px; user-select:none; }
.seq-A { color:#f87171; }
.seq-T { color:#60a5fa; }
.seq-G { color:#4ade80; }
.seq-C { color:#facc15; }
.seq-U { color:#c084fc; }

/* ── GC skew bar ── */
.skew-bar { height:28px; display:flex; border-radius:6px; overflow:hidden; margin:8px 0; }
.skew-pos { background:#22c55e; display:flex; align-items:center; justify-content:center;
  color:#fff; font-size:11px; font-weight:700; }
.skew-neg { background:#ef4444; display:flex; align-items:center; justify-content:center;
  color:#fff; font-size:11px; font-weight:700; }

/* ── Horizontal bar chart ── */
.bar-chart { display:flex; flex-direction:column; gap:5px; }
.bar-row { display:flex; align-items:center; gap:10px; }
.bar-row .bar-label { width:40px; text-align:right; font-size:12px; font-weight:700;
  color:var(--gray); font-family:monospace; }
.bar-row .bar-track { flex:1; height:18px; background:var(--xlgray);
  border-radius:4px; overflow:hidden; }
.bar-row .bar-fill { height:100%; border-radius:4px; display:flex;
  align-items:center; padding-left:6px; font-size:11px; color:#fff; font-weight:600; }
.bar-row .bar-val { width:60px; font-size:12px; color:var(--lgray); }

/* ── Dot-bracket structure ── */
.dot-bracket {
  font-family:'Cascadia Code','Fira Code',monospace; font-size:12px;
  background:#1e293b; color:#a5f3fc; border-radius:6px; padding:12px 14px;
  word-break:break-all; line-height:1.8;
}
.dot-bracket .seq-line { color:#94a3b8; }
.dot-bracket .str-open { color:#f59e0b; }
.dot-bracket .str-close { color:#f59e0b; }
.dot-bracket .str-dot  { color:#475569; }

/* ── Radar-like property circle ── */
.prop-grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(200px,1fr)); gap:10px; }
.prop-item {
  background:var(--surface); border:1px solid var(--border); border-radius:8px;
  padding:12px 14px; display:flex; flex-direction:column; gap:4px;
}
.prop-item .prop-name { font-size:11px; color:var(--lgray); text-transform:uppercase; letter-spacing:.4px; }
.prop-item .prop-val  { font-size:18px; font-weight:700; color:var(--blue); }
.prop-item .prop-unit { font-size:11px; color:var(--lgray); }
.prop-item .prop-bar  { height:5px; background:var(--xlgray); border-radius:3px; margin-top:4px; }
.prop-item .prop-bar-fill { height:100%; border-radius:3px; background:var(--blue); }

/* ── ORF table ── */
.orf-row td:first-child { font-family:monospace; }

/* ── Alert box ── */
.alert { border-radius:8px; padding:12px 16px; margin-bottom:12px;
  display:flex; gap:12px; align-items:flex-start; }
.alert.info  { background:#dbeafe; border-left:4px solid var(--blue); }
.alert.warn  { background:var(--warn); border-left:4px solid #d97706; }
.alert.danger{ background:var(--fail); border-left:4px solid var(--red); }
.alert.good  { background:var(--pass); border-left:4px solid var(--green); }
.alert .alert-icon { font-size:18px; flex-shrink:0; margin-top:1px; }
.alert .alert-body { font-size:13px; color:var(--gray); }
.alert .alert-title { font-weight:700; margin-bottom:2px; }

/* ── Divider ── */
.divider { border:none; border-top:2px solid var(--border); margin:32px 0; }

/* ── Appendix header ── */
.appendix-header {
  background:var(--xlgray); border:1px solid var(--border); border-radius:10px;
  padding:16px 20px; margin-bottom:24px;
}
.appendix-header h2 { font-size:20px; font-weight:800; color:var(--gray); }
.appendix-header p  { font-size:13px; color:var(--lgray); margin-top:4px; }

/* ── Print ── */
@media print {
  body { background:#fff; }
  .page { padding:0; }
  .report-header { border-radius:0; }
  .no-print { display:none; }
}

/* ── Responsive ── */
@media(max-width:600px) {
  .report-header { flex-direction:column; }
  .report-header .meta { text-align:left; }
  .summary-grid { grid-template-columns:repeat(2,1fr); }
}

/* ── 3D Viewer ── */
.viewer3d-container {
  background:#0d1117; border-radius:12px; overflow:hidden;
  border:1px solid #1e3a5f; position:relative; margin-top:8px;
}
.viewer3d-container iframe {
  width:100%; height:400px; border:none; display:block;
}
.viewer3d-overlay {
  position:absolute; top:10px; right:10px; display:flex; gap:8px; z-index:5;
}
.viewer3d-btn {
  background:rgba(31,97,235,.9); color:#fff; border:none; border-radius:16px;
  padding:6px 14px; font-size:11px; font-weight:700; cursor:pointer;
  text-decoration:none; backdrop-filter:blur(6px); transition:background .2s;
}
.viewer3d-btn:hover { background:#1f6feb; }
.viewer3d-btn.full { background:rgba(15,118,110,.9); }
.viewer3d-btn.full:hover { background:#0f766e; }

</style>"""

    // ── HTML page wrapper ─────────────────────────────────────

    private fun page(title: String, body: String) = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>$title — BioKt Report</title>
${sharedStyles()}
</head>
<body>
<div class="page">
$body
</div>
<script>
// Smooth scroll for anchor links
document.querySelectorAll('a[href^="#"]').forEach(a => {
  a.addEventListener('click', e => {
    e.preventDefault();
    document.querySelector(a.getAttribute('href'))?.scrollIntoView({behavior:'smooth'});
  });
});
</script>
</body>
</html>"""

    // ── Helpers ──────────────────────────────────────────────

    private fun Double.fmt(n: Int = 2) = "%.${n}f".format(this)
    private fun badge(text: String, type: String) = """<span class="badge $type">$text</span>"""
    private fun passFailBadge(pass: Boolean, passText: String = "PASS", failText: String = "FAIL") =
        if (pass) badge(passText, "pass") else badge(failText, "fail")

    private fun sectionTitle(icon: String, title: String) =
        """<div class="section-title"><span class="icon">$icon</span> $title</div>"""

    private fun statCard(value: String, label: String, color: String = "") =
        """<div class="stat-card $color"><div class="val">$value</div><div class="lbl">$label</div></div>"""

    private fun coloredSeq(seq: String, maxLen: Int = 120): String {
        val sb = StringBuilder()
        var pos = 0
        seq.take(maxLen).forEachIndexed { i, c ->
            if (i % 60 == 0) {
                if (i > 0) sb.append("<br>")
                sb.append("""<span class="pos">${i + 1}</span>""")
            }
            val cls = when (c) {
                'A' -> "seq-A"; 'T' -> "seq-T"; 'G' -> "seq-G"
                'C' -> "seq-C"; 'U' -> "seq-U"
                else -> ""
            }
            sb.append(if (cls.isNotEmpty()) """<span class="$cls">$c</span>""" else c.toString())
        }
        if (seq.length > maxLen) sb.append("""<span style="color:#64748b"> …+${seq.length - maxLen} more</span>""")
        return sb.toString()
    }

    private fun barChart(data: Map<String, Double>, maxVal: Double, colorFn: (String) -> String = { "#3b82f6" }): String {
        val sb = StringBuilder("""<div class="bar-chart">""")
        data.entries.sortedByDescending { it.value }.forEach { (label, value) ->
            val pct = if (maxVal > 0) (value / maxVal * 100).coerceIn(0.0, 100.0) else 0.0
            sb.append("""
<div class="bar-row">
  <div class="bar-label">$label</div>
  <div class="bar-track">
    <div class="bar-fill" style="width:${pct.fmt(1)}%;background:${colorFn(label)}">
      ${if (pct > 12) "${value.fmt(1)}%" else ""}
    </div>
  </div>
  <div class="bar-val">${value.fmt(2)}%</div>
</div>""")
        }
        sb.append("</div>")
        return sb.toString()
    }

    private fun propGrid(items: List<Triple<String, String, String>>): String {
        val sb = StringBuilder("""<div class="prop-grid">""")
        items.forEach { (name, value, unit) ->
            sb.append("""
<div class="prop-item">
  <div class="prop-name">$name</div>
  <div class="prop-val">$value</div>
  <div class="prop-unit">$unit</div>
</div>""")
        }
        sb.append("</div>")
        return sb.toString()
    }

    private fun alertBox(type: String, icon: String, title: String, body: String) = """
<div class="alert $type">
  <div class="alert-icon">$icon</div>
  <div class="alert-body"><div class="alert-title">$title</div>$body</div>
</div>"""

    private fun dataTable(headers: List<String>, rows: List<List<String>>, colClasses: List<String> = emptyList()): String {
        val sb = StringBuilder("""<table class="data-table"><thead><tr>""")
        headers.forEach { sb.append("<th>$it</th>") }
        sb.append("</tr></thead><tbody>")
        rows.forEach { row ->
            sb.append("<tr>")
            row.forEachIndexed { i, cell ->
                val cls = colClasses.getOrElse(i) { "" }
                sb.append("""<td class="$cls">$cell</td>""")
            }
            sb.append("</tr>")
        }
        sb.append("</tbody></table>")
        return sb.toString()
    }

    private fun tocLink(anchor: String, text: String) =
        """<a href="#$anchor" style="color:var(--blue);text-decoration:none;font-size:13px">▶ $text</a>"""

    private fun nowString(): String {
        val now = java.util.Date()
        return java.text.SimpleDateFormat("dd MMM yyyy HH:mm").format(now)
    }

    // ─────────────────────────────────────────────────────────
    // DNA REPORT
    // ─────────────────────────────────────────────────────────

    private fun buildDNAReport(seq: DNASequence, viewerPath: String = ""): String {
        val gc       = seq.gcContent()
        val tm       = seq.meltingTemperature()
        val rc       = seq.reverseComplement().sequence
        val comp     = seq.complement().sequence
        val rna      = seq.transcribe()
        val protein  = seq.translate()
        val orfs     = seq.findOrfs(minLength = 30)
        val freqs    = SeqStats.nucleotideFrequency(seq)
        val dinucs   = SeqStats.dinucleotideFrequency(seq)
        val entropy  = SeqStats.entropy(seq)
        val complexity = SeqStats.linguisticComplexity(seq)
        val sites    = RestrictionEnzymes.digest(seq)

        val gcClass = when {
            gc < 40 -> "orange"
            gc > 60 -> "purple"
            else    -> "green"
        }
        val tmClass = when {
            tm < 50 -> "orange"
            tm > 75 -> "red"
            else    -> "teal"
        }

        // ── TOC ──────────────────────────────────────────────
        val toc = """
<div style="background:var(--surface);border:1px solid var(--border);border-radius:10px;
  padding:16px 20px;margin-bottom:28px;" class="no-print">
  <div style="font-weight:700;color:var(--gray);margin-bottom:10px;font-size:13px;">
    📋 Contents
  </div>
  <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:6px;">
    ${tocLink("summary", "Executive Summary")}
    ${tocLink("sequence", "Sequence Viewer")}
    ${tocLink("composition", "Base Composition")}
    ${tocLink("structure", "Structural Properties")}
    ${tocLink("orfs", "Open Reading Frames")}
    ${tocLink("restriction", "Restriction Sites")}
    ${tocLink("translation", "Transcription & Translation")}
    ${if (viewerPath.isNotEmpty()) tocLink("viewer3d", "3D Structure Viewer") else ""}
    ${tocLink("appendix", "Detailed Appendix")}
  </div>
</div>"""

        // ── Header ───────────────────────────────────────────
        val header = """
<div class="report-header">
  <div>
    <div class="report-type-badge">DNA Analysis Report</div>
    <h1>${seq.id.ifEmpty { "Sequence" }}</h1>
    <div class="subtitle">${seq.description.ifEmpty { "Deoxyribonucleic acid sequence analysis" }}</div>
  </div>
  <div class="meta">
    Generated by BioKt v2.0<br>
    ${nowString()}<br>
    Length: ${seq.length} bp
  </div>
</div>"""

        // ── Summary ───────────────────────────────────────────
        val summary = """
<div class="section" id="summary">
  ${sectionTitle("📊", "Executive Summary")}
  <div class="summary-grid">
    ${statCard("${seq.length}", "Length (bp)")}
    ${statCard("${gc.fmt(1)}%", "GC Content", gcClass)}
    ${statCard("${tm.fmt(1)}°C", "Melting Temp", tmClass)}
    ${statCard("${orfs.size}", "ORFs Found", if (orfs.isEmpty()) "orange" else "green")}
    ${statCard("${sites.size}", "Restriction Sites", "purple")}
    ${statCard("${protein.length}", "Protein Length (aa)", "teal")}
    ${statCard(entropy.fmt(3), "Shannon Entropy")}
    ${statCard("${(complexity * 100).fmt(1)}%", "Complexity", if (complexity < 0.5) "orange" else "green")}
  </div>
  ${when {
      gc < 30 -> alertBox("warn", "⚠️", "Low GC Content", "GC content of ${gc.fmt(1)}% is below 30%. This sequence may have reduced thermal stability and could be AT-rich repetitive DNA.")
      gc > 70 -> alertBox("warn", "⚠️", "High GC Content", "GC content of ${gc.fmt(1)}% is above 70%. This is characteristic of certain bacterial or GC-rich genomic regions.")
      else    -> alertBox("good", "✅", "GC Content Normal", "GC content of ${gc.fmt(1)}% is within the typical range (40–60%) for most coding sequences.")
  }}
  ${if (complexity < 0.4) alertBox("warn", "⚠️", "Low Complexity", "Linguistic complexity of ${(complexity*100).fmt(1)}% suggests repetitive sequence content (e.g. tandem repeats, low-complexity regions).") else ""}
  ${if (orfs.isEmpty()) alertBox("warn", "⚠️", "No ORFs Found", "No open reading frames of ≥30 bp were detected. This may be a non-coding RNA, regulatory, or intergenic sequence.") else ""}
</div>"""

        // ── Sequence viewer ────────────────────────────────────
        val seqViewer = """
<div class="section" id="sequence">
  ${sectionTitle("🔬", "Sequence Viewer")}
  <div style="margin-bottom:8px;display:flex;gap:8px;flex-wrap:wrap;">
    <span style="font-size:12px;color:var(--lgray)">Legend:</span>
    <span style="font-size:12px;color:#f87171;font-family:monospace;font-weight:700">■ A</span>
    <span style="font-size:12px;color:#60a5fa;font-family:monospace;font-weight:700">■ T</span>
    <span style="font-size:12px;color:#4ade80;font-family:monospace;font-weight:700">■ G</span>
    <span style="font-size:12px;color:#facc15;font-family:monospace;font-weight:700">■ C</span>
  </div>
  <div class="seq-display">${coloredSeq(seq.sequence, 300)}</div>
  <div style="margin-top:10px;font-size:12px;color:var(--lgray)">
    Showing up to 300 bp. Full length: ${seq.length} bp.
  </div>
</div>"""

        // ── Base composition ───────────────────────────────────
        val baseColors = mapOf("A" to "#f87171", "T" to "#60a5fa", "G" to "#4ade80", "C" to "#facc15")
        val compSection = """
<div class="section" id="composition">
  ${sectionTitle("📈", "Base Composition")}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:24px;flex-wrap:wrap;">
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Mononucleotide Frequencies</div>
      ${barChart(
            freqs.map { it.key.toString() to it.value * 100 }.toMap(),
            ((if (freqs.values.isEmpty()) 0.0 else freqs.values.max()!!)) * 100
        ) { base -> baseColors[base] ?: "#3b82f6" }}
    </div>
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Top 8 Dinucleotides</div>
      ${barChart(
            dinucs.entries.sortedByDescending { it.value }.take(8)
                .associate { it.key to it.value * 100 },
            ((if (dinucs.values.isEmpty()) 0.0 else dinucs.values.max()!!)) * 100
        )}
    </div>
  </div>
  <div style="margin-top:16px;">
  ${dataTable(
        listOf("Base", "Count", "Frequency", "Expected (50% GC)"),
        listOf("A", "T", "G", "C").map { base ->
            val f = freqs[base[0]] ?: 0.0
            val count = (f * seq.length).toInt()
            val expected = if (base == "G" || base == "C") gc / 2 else (100 - gc) / 2
            listOf(base, count.toString(), "${(f * 100).fmt(2)}%", "${expected.fmt(1)}%")
        },
        listOf("", "num", "num", "num")
    )}
  </div>
</div>"""

        // ── Structural properties ──────────────────────────────
        val gcSkew = seq.gcSkew(windowSize = minOf(50, seq.length / 4 + 1), stepSize = minOf(10, seq.length / 10 + 1))
        val avgSkew = if (gcSkew.isNotEmpty()) gcSkew.map { it.second }.average() else 0.0
        val structSection = """
<div class="section" id="structure">
  ${sectionTitle("🧬", "Structural Properties")}
  ${propGrid(listOf(
        Triple("Melting Temperature", tm.fmt(1), "°C"),
        Triple("GC Content", gc.fmt(2), "%"),
        Triple("AT Content", (100 - gc).fmt(2), "%"),
        Triple("GC Skew (avg)", avgSkew.fmt(4), "(G−C)/(G+C)"),
        Triple("Shannon Entropy", entropy.fmt(4), "bits"),
        Triple("Complexity", (complexity * 100).fmt(1), "%"),
        Triple("Length", seq.length.toString(), "bp"),
        Triple("Complement", "5'→3'", "sense")
  ))}
  <div style="margin-top:20px;">
    <div style="font-weight:600;color:var(--gray);margin-bottom:8px;font-size:13px;">GC Skew Profile</div>
    <div style="font-size:12px;color:var(--lgray);margin-bottom:8px;">
      Positive skew (green) = more G than C. Negative skew (red) = more C than G.
      Skew shifts mark replication origins/termini.
    </div>
    <div style="display:flex;flex-direction:column;gap:4px;">
    ${gcSkew.take(20).joinToString("") { (pos, skew) ->
        val pct = ((skew + 1) / 2 * 100).coerceIn(5.0, 95.0)
        val color = if (skew >= 0) "#22c55e" else "#ef4444"
        """<div style="display:flex;align-items:center;gap:8px;">
          <span style="font-size:10px;color:var(--lgray);width:40px;text-align:right">$pos</span>
          <div style="flex:1;height:12px;background:#f1f5f9;border-radius:3px;position:relative;">
            <div style="position:absolute;left:50%;top:0;width:1px;height:100%;background:#94a3b8;"></div>
            ${if (skew >= 0)
                """<div style="position:absolute;left:50%;width:${(skew*50).coerceIn(0.0,50.0).fmt(1)}%;height:100%;background:$color;border-radius:0 3px 3px 0;"></div>"""
            else
                """<div style="position:absolute;right:${(50 + skew*50).coerceIn(0.0,50.0).fmt(1)}%;width:${(-skew*50).coerceIn(0.0,50.0).fmt(1)}%;height:100%;background:$color;border-radius:3px 0 0 3px;"></div>"""}
          </div>
          <span style="font-size:10px;color:${if (skew>=0) "#16a34a" else "#dc2626"};width:55px">${if (skew>=0) "+" else ""}${skew.fmt(3)}</span>
        </div>"""
    }}
    </div>
  </div>
  <div style="margin-top:16px;">
    <div style="font-weight:600;color:var(--gray);margin-bottom:8px;font-size:13px;">Reverse Complement</div>
    <div class="seq-display">${coloredSeq(rc, 150)}</div>
  </div>
</div>"""

        // ── ORFs ──────────────────────────────────────────────
        val orfSection = """
<div class="section" id="orfs">
  ${sectionTitle("🔍", "Open Reading Frames")}
  ${if (orfs.isEmpty())
        alertBox("warn", "⚠️", "No ORFs Found", "No ORFs of ≥30 bp detected in the forward strand.")
    else """
  <div style="margin-bottom:12px;font-size:13px;color:var(--lgray)">
    Found ${orfs.size} ORF(s) in forward reading frames (≥30 bp threshold).
  </div>
  ${dataTable(
        listOf("Frame", "Start", "End", "Length (bp)", "Protein (aa)", "Preview"),
        orfs.map { orf ->
            listOf(
                "+${orf.frame}",
                orf.start.toString(),
                orf.end.toString(),
                orf.length.toString(),
                orf.protein.length.toString(),
                orf.protein.sequence.take(20) + if (orf.protein.length > 20) "…" else ""
            )
        },
        listOf("", "num", "num", "num", "num", "mono")
    )}"""}
</div>"""

        // ── Restriction sites ──────────────────────────────────
        val restrictSection = """
<div class="section" id="restriction">
  ${sectionTitle("✂️", "Restriction Enzyme Sites")}
  ${if (sites.isEmpty())
        alertBox("info", "ℹ️", "No Sites Found", "None of the 25 built-in restriction enzyme recognition sequences were found in this sequence.")
    else """
  <div style="margin-bottom:12px;font-size:13px;color:var(--lgray)">
    Found ${sites.size} restriction site(s) across ${sites.map { it.enzyme }.distinct().size} enzyme(s).
  </div>
  ${dataTable(
        listOf("Enzyme", "Position", "Recognition Sequence", "Cut Offset"),
        sites.sortedBy { it.position }.map { site ->
            listOf(site.enzyme, site.position.toString(), site.recognition, "—")
        },
        listOf("", "num", "mono", "num")
    )}"""}
</div>"""

        // ── Translation ───────────────────────────────────────
        val transSection = """
<div class="section" id="translation">
  ${sectionTitle("🔄", "Transcription & Translation")}
  <div style="margin-bottom:16px;">
    <div style="font-weight:600;color:var(--gray);margin-bottom:6px;font-size:13px;">
      mRNA Transcript (5'→3')
    </div>
    <div class="seq-display" style="color:#c084fc">${rna.sequence.take(200)}${if (rna.length > 200) "…" else ""}</div>
  </div>
  <div>
    <div style="font-weight:600;color:var(--gray);margin-bottom:6px;font-size:13px;">
      Translated Protein (Frame +1)
    </div>
    <div class="seq-display" style="color:#86efac">${protein.sequence}</div>
    <div style="margin-top:8px;display:flex;gap:12px;flex-wrap:wrap;font-size:12px;color:var(--lgray)">
      <span>Length: ${protein.length} aa</span>
      <span>MW: ~${protein.molecularWeight().fmt(1)} Da</span>
      <span>pI: ${protein.isoelectricPoint().fmt(2)}</span>
    </div>
  </div>
</div>"""

        // ── Appendix ──────────────────────────────────────────
        val appendix = """
<hr class="divider">
<div id="appendix">
  <div class="appendix-header">
    <h2>📎 Detailed Appendix</h2>
    <p>Extended analysis data, full tables, and computed values</p>
  </div>

  <div class="section">
    <div class="section-title">A1 — Full Dinucleotide Frequency Table</div>
    ${dataTable(
        listOf("Dinucleotide", "Count", "Observed %", "Expected % (random)"),
        dinucs.entries.sortedByDescending { it.value }.map { (di, f) ->
            val count = (f * (seq.length - 1)).toInt()
            val baseFreqA = (freqs[di[0]] ?: 0.25)
            val baseFreqB = (freqs[di[1]] ?: 0.25)
            val expected  = (baseFreqA * baseFreqB * 100).fmt(3)
            listOf(di, count.toString(), (f * 100).fmt(4), expected)
        },
        listOf("mono", "num", "num", "num")
    )}
  </div>

  <div class="section">
    <div class="section-title">A2 — Codon Usage (Reading Frame +1)</div>
    ${if (seq.length >= 3) {
        val codons = SeqStats.codonUsageTable(seq)
        if (codons.isNotEmpty()) dataTable(
            listOf("Codon", "Amino Acid", "Count"),
            codons.entries.sortedByDescending { it.value }.take(20).map { (codon, count) ->
                listOf(codon, CodonTables.STANDARD.translate(codon).toString(), count.toString())
            },
            listOf("mono", "", "num")
        ) else "<p style='color:var(--lgray);font-size:13px'>Sequence too short for codon analysis.</p>"
    } else "<p style='color:var(--lgray);font-size:13px'>Sequence too short for codon analysis.</p>"}
  </div>

  <div class="section">
    <div class="section-title">A3 — Six-Frame Translation</div>
    ${dataTable(
        listOf("Frame", "Protein Sequence (first 40 aa)", "Length"),
        seq.sixFrameTranslation().entries.map { (frame, prot) ->
            listOf(frame, prot.sequence.take(40) + if (prot.length > 40) "…" else "", "${prot.length} aa")
        },
        listOf("", "mono", "num")
    )}
  </div>

  <div class="section">
    <div class="section-title">A4 — Full Sequence Information</div>
    ${dataTable(
        listOf("Property", "Value"),
        listOf(
            listOf("Sequence ID", seq.id.ifEmpty { "—" }),
            listOf("Description", seq.description.ifEmpty { "—" }),
            listOf("Length", "${seq.length} bp"),
            listOf("GC Content", "${gc.fmt(4)}%"),
            listOf("AT Content", "${(100 - gc).fmt(4)}%"),
            listOf("Melting Temperature", "${tm.fmt(2)} °C"),
            listOf("GC Skew (average)", avgSkew.fmt(6)),
            listOf("Shannon Entropy", entropy.fmt(6)),
            listOf("Linguistic Complexity", (complexity * 100).fmt(4) + "%"),
            listOf("ORFs (≥30 bp)", orfs.size.toString()),
            listOf("Restriction Sites", sites.size.toString()),
            listOf("Translated protein (aa)", protein.length.toString()),
            listOf("Protein MW", "${protein.molecularWeight().fmt(2)} Da"),
            listOf("Protein pI", protein.isoelectricPoint().fmt(3))
        ),
        listOf("", "")
    )}
  </div>
</div>"""

        val viewer3d = if (viewerPath.isNotEmpty()) dna3dSection(seq, viewerPath) else ""
        val body = header + toc + summary + viewer3d + seqViewer + compSection + structSection +
                   orfSection + restrictSection + transSection + appendix

        return page("DNA Report — ${seq.id.ifEmpty { "Sequence" }}", body)
    }

    // ─────────────────────────────────────────────────────────
    // PROTEIN REPORT
    // ─────────────────────────────────────────────────────────

    private fun buildProteinReport(seq: ProteinSequence, viewerPath: String = ""): String {
        val mw       = seq.molecularWeight()
        val pi       = seq.isoelectricPoint()
        val arom     = seq.aromaticity()
        val instab   = seq.instabilityIndex()
        val gravy    = ProteinDescriptors.gravyIndex(seq)
        val aliphat  = ProteinDescriptors.aliphaticIndex(seq)
        val boman    = ProteinDescriptors.bomanIndex(seq)
        val charge7  = ProteinDescriptors.chargeAtPH(seq, 7.0)
        val hydroR   = ProteinDescriptors.hydrophobicRatio(seq)
        val extCoef  = ProteinDescriptors.extinctionCoefficient(seq)
        val aac      = ProteinDescriptors.aminoAcidComposition(seq)
        val tm       = ProteinInteraction.predictTransmembraneHelices(seq)
        val cc       = ProteinInteraction.predictCoiledCoil(seq)
        val iface    = ProteinInteraction.predictInterfaceResidues(seq)
        val stable   = instab < 40.0

        val header = """
<div class="report-header">
  <div>
    <div class="report-type-badge">Protein Analysis Report</div>
    <h1>${seq.id.ifEmpty { "Protein" }}</h1>
    <div class="subtitle">${seq.description.ifEmpty { "Protein sequence physicochemical analysis" }}</div>
  </div>
  <div class="meta">
    Generated by BioKt v2.0<br>
    ${nowString()}<br>
    Length: ${seq.length} aa
  </div>
</div>"""

        val toc = """
<div style="background:var(--surface);border:1px solid var(--border);border-radius:10px;
  padding:16px 20px;margin-bottom:28px;" class="no-print">
  <div style="font-weight:700;color:var(--gray);margin-bottom:10px;font-size:13px;">📋 Contents</div>
  <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:6px;">
    ${tocLink("prot-summary", "Executive Summary")}
    ${tocLink("prot-sequence", "Sequence Viewer")}
    ${tocLink("prot-physchem", "Physicochemical Properties")}
    ${tocLink("prot-composition", "Amino Acid Composition")}
    ${tocLink("prot-structure", "Structural Prediction")}
    ${if (viewerPath.isNotEmpty()) tocLink("viewer3d", "3D Structure Viewer") else ""}
    ${tocLink("prot-appendix", "Detailed Appendix")}
  </div>
</div>"""

        val summary = """
<div class="section" id="prot-summary">
  ${sectionTitle("📊", "Executive Summary")}
  <div class="summary-grid">
    ${statCard("${seq.length}", "Length (aa)")}
    ${statCard("${(mw / 1000).fmt(2)} kDa", "Molecular Weight", "blue")}
    ${statCard(pi.fmt(2), "Isoelectric Point", if (pi > 7) "orange" else "teal")}
    ${statCard("${(arom * 100).fmt(1)}%", "Aromaticity", "purple")}
    ${statCard(gravy.fmt(3), "GRAVY", if (gravy > 0) "orange" else "teal")}
    ${statCard(if (stable) "STABLE" else "UNSTABLE", "Stability", if (stable) "green" else "red")}
    ${statCard("${tm.size}", "TM Helices", if (tm.isNotEmpty()) "orange" else "teal")}
    ${statCard("${cc.size}", "Coiled Coils", if (cc.isNotEmpty()) "purple" else "")}
  </div>

  ${if (stable) alertBox("good", "✅", "Protein Predicted Stable", "Instability index of ${instab.fmt(1)} < 40. Protein is predicted to be stable in vitro.")
    else alertBox("warn", "⚠️", "Protein Predicted Unstable", "Instability index of ${instab.fmt(1)} ≥ 40. May degrade rapidly in vitro.")}
  ${when {
      pi < 6.0 -> alertBox("info", "ℹ️", "Acidic Protein (pI ${pi.fmt(2)})", "Net negative charge at physiological pH. More acidic residues (Asp, Glu) than basic (Lys, Arg, His).")
      pi > 8.5 -> alertBox("info", "ℹ️", "Basic Protein (pI ${pi.fmt(2)})", "Net positive charge at physiological pH. Typical of DNA-binding proteins, histones, and antimicrobial peptides.")
      else     -> alertBox("info", "ℹ️", "Near-Neutral pI (${pi.fmt(2)})", "Isoelectric point close to physiological pH (6–8.5). Minimal net charge at pH 7.4.")
  }}
  ${if (gravy > 0.0) alertBox("warn", "⚠️", "Hydrophobic Protein (GRAVY ${gravy.fmt(3)})", "Positive GRAVY index indicates overall hydrophobicity. May be membrane-associated or poorly soluble.") else ""}
  ${if (tm.isNotEmpty()) alertBox("warn", "⚠️", "${tm.size} Transmembrane Helix(es) Predicted", "Kyte-Doolittle window analysis detected probable TM-spanning helices. Protein may be membrane-embedded.") else ""}
</div>"""

        val seqViewer = """
<div class="section" id="prot-sequence">
  ${sectionTitle("🔬", "Sequence Viewer")}
  <div class="seq-display" style="color:#a5f3fc;font-size:12px;">
    ${seq.sequence.chunked(60).mapIndexed { i, chunk ->
        """<span class="pos">${i * 60 + 1}</span>$chunk"""
    }.joinToString("<br>")}
  </div>
  <div style="margin-top:10px;font-size:12px;color:var(--lgray)">Length: ${seq.length} amino acids</div>
</div>"""

        val physChemSection = """
<div class="section" id="prot-physchem">
  ${sectionTitle("⚗️", "Physicochemical Properties")}
  ${propGrid(listOf(
        Triple("Molecular Weight",    (mw/1000).fmt(3),        "kDa"),
        Triple("Isoelectric Point",   pi.fmt(3),               "pH units"),
        Triple("GRAVY Index",         gravy.fmt(4),            "hydropathicity"),
        Triple("Aromaticity",         (arom*100).fmt(2),       "%"),
        Triple("Instability Index",   instab.fmt(2),           if (stable) "STABLE" else "UNSTABLE"),
        Triple("Aliphatic Index",     aliphat.fmt(2),          "thermostability"),
        Triple("Boman Index",         boman.fmt(3),            "kcal/mol"),
        Triple("Charge at pH 7",      charge7.fmt(3),          "net charge"),
        Triple("Hydrophobic Ratio",   (hydroR*100).fmt(1),     "%"),
        Triple("Extinction Coef",     extCoef.fmt(0),          "M⁻¹cm⁻¹")
  ))}
  <div style="margin-top:16px;">
  ${dataTable(
        listOf("Property", "Value", "Interpretation"),
        listOf(
            listOf("GRAVY", gravy.fmt(3), if (gravy > 0) "Hydrophobic — potential membrane association" else "Hydrophilic — good aqueous solubility expected"),
            listOf("Instability Index", instab.fmt(2), if (stable) "< 40 → Stable in vitro" else "≥ 40 → Likely unstable in vitro"),
            listOf("Boman Index", boman.fmt(3), when { boman > 2.48 -> "High — strong potential for cell penetration / AMPs" ; boman < 1.0 -> "Low — unlikely to bind other proteins strongly" ; else -> "Moderate" }),
            listOf("Aliphatic Index", aliphat.fmt(2), if (aliphat > 80) "High — increased thermostability" else "Moderate thermostability"),
            listOf("Charge at pH 7", charge7.fmt(3), if (charge7 > 0) "Net positive — cationic protein" else if (charge7 < 0) "Net negative — anionic protein" else "Near-zero net charge")
        ),
        listOf("", "num", "")
    )}
  </div>
</div>"""

        val aaColors = mapOf(
            "R" to "#f87171", "K" to "#f87171", "H" to "#f87171",
            "D" to "#60a5fa", "E" to "#60a5fa",
            "G" to "#94a3b8", "A" to "#94a3b8", "V" to "#94a3b8", "L" to "#94a3b8", "I" to "#94a3b8", "P" to "#94a3b8",
            "S" to "#4ade80", "T" to "#4ade80", "C" to "#4ade80", "M" to "#4ade80",
            "F" to "#c084fc", "Y" to "#c084fc", "W" to "#c084fc",
            "N" to "#facc15", "Q" to "#facc15"
        )

        val compSection = """
<div class="section" id="prot-composition">
  ${sectionTitle("📈", "Amino Acid Composition")}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:24px;">
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">All 20 Amino Acids</div>
      ${barChart(
            aac.entries.sortedByDescending { it.value }.associate { (k, v) -> k to v * 100 },
            (if (aac.values.isEmpty()) 0.0 else aac.values.max()!!) * 100
        ) { aa -> aaColors[aa] ?: "#3b82f6" }}
    </div>
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:12px;font-size:13px;">Residue Categories</div>
      ${run {
            val charged = listOf('R','K','H','D','E').map { seq.count(it) }.sum().toDouble() / seq.length * 100
            val polar   = listOf('S','T','N','Q','C','Y').map { seq.count(it) }.sum().toDouble() / seq.length * 100
            val nonpolar= listOf('G','A','V','L','I','P','F','M','W').map { seq.count(it) }.sum().toDouble() / seq.length * 100
            val aromatic= listOf('F','Y','W','H').map { seq.count(it) }.sum().toDouble() / seq.length * 100
            barChart(mapOf("Charged" to charged, "Polar" to polar, "Nonpolar" to nonpolar, "Aromatic" to aromatic), 100.0)
        }}
      <div style="margin-top:16px;">
      ${dataTable(
            listOf("Category", "Residues", "%"),
            listOf(
                listOf("Basic (pos.)", "R, K, H", "${(listOf('R','K','H').map { seq.count(it) }.sum().toDouble() / seq.length * 100).fmt(1)}%"),
                listOf("Acidic (neg.)", "D, E",   "${(listOf('D','E').map { seq.count(it) }.sum().toDouble() / seq.length * 100).fmt(1)}%"),
                listOf("Hydrophobic",   "A,V,I,L,M,F,W,P", "${(hydroR * 100).fmt(1)}%"),
                listOf("Aromatic",      "F, Y, W, H", "${(arom * 100).fmt(1)}%")
            ),
            listOf("", "mono", "num")
        )}
      </div>
    </div>
  </div>
</div>"""

        val structSection = """
<div class="section" id="prot-structure">
  ${sectionTitle("🏗️", "Structural Predictions")}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:20px;">
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:8px;font-size:13px;">Transmembrane Helices</div>
      ${if (tm.isEmpty()) alertBox("good", "✅", "No TM helices", "Protein is predicted to be soluble (non-membrane).")
        else """<div>${tm.mapIndexed { i, r -> """<div style="background:var(--warn);border-radius:6px;padding:8px 12px;margin-bottom:6px;font-size:13px;">
          <strong>TM Helix ${i+1}:</strong> residues ${r.first}–${r.last} (${r.last - r.first + 1} aa)
        </div>"""}.joinToString("")}</div>"""}
    </div>
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:8px;font-size:13px;">Coiled-Coil Segments</div>
      ${if (cc.isEmpty()) alertBox("info", "ℹ️", "No coiled coils", "No heptad repeat patterns detected.")
        else """<div>${cc.mapIndexed { i, seg -> """<div style="background:#ede9fe;border-radius:6px;padding:8px 12px;margin-bottom:6px;font-size:13px;">
          <strong>Segment ${i+1}:</strong> res ${seg.start}–${seg.end}
          &nbsp;score=${seg.score.fmt(3)}&nbsp;P=${seg.probability.fmt(3)}
        </div>"""}.joinToString("")}</div>"""}
    </div>
  </div>
  <div style="font-weight:600;color:var(--gray);margin-bottom:8px;font-size:13px;">
    Interface Propensity (top 10 residues)
  </div>
  ${dataTable(
        listOf("Position", "Residue", "Propensity Score", "Interface?"),
        iface.sortedByDescending { it.propensityScore }.take(10).map { r ->
            listOf(r.position.toString(), r.aminoAcid.toString(),
                r.propensityScore.fmt(4),
                if (r.isLikelyInterface) badge("YES", "pass") else badge("NO", "fail"))
        },
        listOf("num", "mono", "num", "")
    )}
</div>"""

        val appendix = """
<hr class="divider">
<div id="prot-appendix">
  <div class="appendix-header">
    <h2>📎 Detailed Appendix</h2>
    <p>Full descriptor tables, CTD composition, and charge profile</p>
  </div>

  <div class="section">
    <div class="section-title">A1 — Full Amino Acid Composition Table</div>
${run {
        val aaNames = mapOf('A' to "Alanine",'R' to "Arginine",'N' to "Asparagine",'D' to "Aspartate",
            'C' to "Cysteine",'Q' to "Glutamine",'E' to "Glutamate",'G' to "Glycine",
            'H' to "Histidine",'I' to "Isoleucine",'L' to "Leucine",'K' to "Lysine",
            'M' to "Methionine",'F' to "Phenylalanine",'P' to "Proline",'S' to "Serine",
            'T' to "Threonine",'W' to "Tryptophan",'Y' to "Tyrosine",'V' to "Valine")
        val aaCategories = mapOf('R' to "Basic",'K' to "Basic",'H' to "Basic",
            'D' to "Acidic",'E' to "Acidic",
            'A' to "Hydrophobic",'V' to "Hydrophobic",'I' to "Hydrophobic",
            'L' to "Hydrophobic",'M' to "Hydrophobic",'F' to "Hydrophobic",
            'W' to "Hydrophobic",'P' to "Hydrophobic",
            'G' to "Neutral",'S' to "Polar",'T' to "Polar",'C' to "Polar",
            'Y' to "Polar",'N' to "Polar",'Q' to "Polar")
        dataTable(
            listOf("AA", "Name", "Count", "Frequency", "Category"),
            aac.entries.sortedByDescending { it.value }.map { (aa, freq) ->
                val aaChar = aa[0]
                val count = (freq * seq.length).toInt()
                listOf(aa, aaNames[aaChar] ?: "—", count.toString(),
                    "${(freq * 100).fmt(3)}%", aaCategories[aaChar] ?: "—")
            },
            listOf("mono", "", "num", "num", "")
        )
    }}
  </div>

  <div class="section">
    <div class="section-title">A2 — Charge Profile Across pH</div>
    ${dataTable(
        listOf("pH", "Net Charge", "Status"),
        listOf(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.4, 8.0, 9.0, 10.0, 11.0, 12.0).map { ph ->
            val chg = ProteinDescriptors.chargeAtPH(seq, ph)
            listOf(ph.toString(), chg.fmt(3),
                when { chg > 0.5 -> badge("Positive", "pass") ; chg < -0.5 -> badge("Negative", "fail") ; else -> badge("Neutral", "info") })
        },
        listOf("num", "num", "")
    )}
  </div>

  <div class="section">
    <div class="section-title">A3 — Summary Statistics</div>
    ${dataTable(
        listOf("Property", "Value"),
        listOf(
            listOf("Sequence ID", seq.id.ifEmpty { "—" }),
            listOf("Length", "${seq.length} aa"),
            listOf("Molecular Weight", "${mw.fmt(2)} Da (${(mw/1000).fmt(3)} kDa)"),
            listOf("Isoelectric Point (pI)", pi.fmt(4)),
            listOf("GRAVY Index", gravy.fmt(5)),
            listOf("Aromaticity", "${(arom*100).fmt(3)}%"),
            listOf("Instability Index", "${instab.fmt(3)} (${if(stable) "Stable" else "Unstable"})"),
            listOf("Aliphatic Index", aliphat.fmt(3)),
            listOf("Boman Index", boman.fmt(4)),
            listOf("Net Charge at pH 7.0", charge7.fmt(4)),
            listOf("Net Charge at pH 7.4", ProteinDescriptors.chargeAtPH(seq, 7.4).fmt(4)),
            listOf("Hydrophobic Residues", "${(hydroR*100).fmt(2)}%"),
            listOf("Extinction Coefficient", "${extCoef.fmt(0)} M⁻¹cm⁻¹"),
            listOf("Predicted TM Helices", tm.size.toString()),
            listOf("Predicted Coiled Coils", cc.size.toString())
        ),
        listOf("", "")
    )}
  </div>
</div>"""

        val viewer3d = if (viewerPath.isNotEmpty()) protein3dSection(seq, viewerPath) else ""
        val body = header + toc + summary + viewer3d + seqViewer + physChemSection + compSection + structSection + appendix
        return page("Protein Report — ${seq.id.ifEmpty { "Protein" }}", body)
    }

    // ─────────────────────────────────────────────────────────
    // DRUG / MOLECULE REPORT
    // ─────────────────────────────────────────────────────────

    private fun buildDrugReport(mol: Molecule, viewerPath: String = ""): String {
        val desc    = MolDescriptors.calculate(mol)
        val fp      = Fingerprints.morgan(mol)
        val maccs   = Fingerprints.maccs(mol)
        val alerts  = MolDescriptors.checkStructuralAlerts(mol)
        val pharm   = DrugProteinBinding.generatePharmacophore(mol)
        val admet   = desc.admet

        val header = """
<div class="report-header">
  <div>
    <div class="report-type-badge">Drug / Molecule Report</div>
    <h1>${mol.name.ifEmpty { mol.id.ifEmpty { "Molecule" } }}</h1>
    <div class="subtitle">SMILES: <span style="font-family:monospace;font-size:13px;opacity:.9">${mol.smiles}</span></div>
  </div>
  <div class="meta">
    Generated by BioKt v2.0<br>
    ${nowString()}<br>
    Formula: ${mol.molecularFormula()}
  </div>
</div>"""

        val toc = """
<div style="background:var(--surface);border:1px solid var(--border);border-radius:10px;
  padding:16px 20px;margin-bottom:28px;" class="no-print">
  <div style="font-weight:700;color:var(--gray);margin-bottom:10px;font-size:13px;">📋 Contents</div>
  <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:6px;">
    ${tocLink("drug-summary", "Executive Summary")}
    ${tocLink("drug-descriptors", "Molecular Descriptors")}
    ${tocLink("drug-filters", "Drug-Likeness Filters")}
    ${tocLink("drug-admet", "ADMET Profile")}
    ${tocLink("drug-structure", "Structural Analysis")}
    ${tocLink("drug-alerts", "Structural Alerts")}
    ${tocLink("drug-pharm", "Pharmacophore")}
    ${if (viewerPath.isNotEmpty()) tocLink("viewer3d", "3D Structure Viewer") else ""}
    ${tocLink("drug-appendix", "Detailed Appendix")}
  </div>
</div>"""

        // Pass/fail counts
        val filtersPassed = listOf(desc.lipinskiPasses, desc.veberPasses, desc.ghosePasses, desc.eganPasses, desc.mueggeScore >= 7).count { it }

        val summary = """
<div class="section" id="drug-summary">
  ${sectionTitle("📊", "Executive Summary")}
  <div class="summary-grid">
    ${statCard(desc.molecularFormula, "Formula")}
    ${statCard("${desc.molecularWeight.fmt(2)}", "MW (Da)", "blue")}
    ${statCard(desc.logP.fmt(2), "LogP", if (desc.logP in 0.0..3.0) "green" else if (desc.logP > 5.0) "red" else "orange")}
    ${statCard("${desc.topologicalPolarSurfaceArea.fmt(1)} Å²", "TPSA", if (desc.topologicalPolarSurfaceArea <= 90) "green" else if (desc.topologicalPolarSurfaceArea > 140) "red" else "orange")}
    ${statCard("$filtersPassed/5", "Filters Passed", if (filtersPassed >= 4) "green" else if (filtersPassed >= 2) "orange" else "red")}
    ${statCard(if (admet.bbbPenetration) "YES" else "NO", "BBB", if (admet.bbbPenetration) "green" else "orange")}
    ${statCard(if (admet.hergInhibition) "RISK" else "SAFE", "hERG", if (admet.hergInhibition) "red" else "green")}
    ${statCard("${alerts.size}", "Alerts", if (alerts.isEmpty()) "green" else "red")}
  </div>

  ${when {
      filtersPassed >= 4 && alerts.isEmpty() ->
          alertBox("good", "✅", "Drug-Likeness: Excellent", "Passes $filtersPassed/5 drug-likeness filters with no structural alerts. Good candidate for oral drug development.")
      filtersPassed >= 3 && alerts.size <= 1 ->
          alertBox("info", "ℹ️", "Drug-Likeness: Acceptable", "Passes $filtersPassed/5 filters. Minor concerns noted — review alerts and ADMET before proceeding.")
      else ->
          alertBox("warn", "⚠️", "Drug-Likeness: Concerns", "Only $filtersPassed/5 filters passed. Significant modifications may be required for oral bioavailability.")
  }}
  ${if (admet.hergInhibition) alertBox("danger", "❌", "hERG Inhibition Risk", "Structural features (lipophilic + basic nitrogen) suggest possible hERG channel inhibition. Cardiac safety testing required.") else ""}
  ${if (alerts.any { it.severity == "High" }) alertBox("danger", "❌", "High-Severity Structural Alert", "One or more high-severity reactive groups detected (Michael acceptor, nitro group, etc.). Review before proceeding.") else ""}
</div>"""

        val descriptorSection = """
<div class="section" id="drug-descriptors">
  ${sectionTitle("🔬", "Molecular Descriptors")}
  ${propGrid(listOf(
        Triple("Molecular Weight",    desc.molecularWeight.fmt(2),            "Da"),
        Triple("Exact Mass",          desc.exactMass.fmt(4),                  "Da (monoisotopic)"),
        Triple("LogP (Crippen)",      desc.logP.fmt(3),                       "lipophilicity"),
        Triple("TPSA (Ertl)",         desc.topologicalPolarSurfaceArea.fmt(1),"Å²"),
        Triple("H-Bond Donors",       desc.hBondDonors.toString(),            "N-H, O-H"),
        Triple("H-Bond Acceptors",    desc.hBondAcceptors.toString(),         "N, O"),
        Triple("Rotatable Bonds",     desc.numRotatableBonds.toString(),      "flexibility"),
        Triple("Rings",               desc.numRings.toString(),               "(aromatic: ${desc.numAromaticRings})"),
        Triple("Heavy Atoms",         desc.numHeavyAtoms.toString(),          "non-hydrogen"),
        Triple("fsp3",                desc.fractionCSP3.fmt(3),               "sp3 carbon fraction"),
        Triple("Molar Refractivity",  desc.molarRefractivity.fmt(2),          "cm³/mol"),
        Triple("Stereocenters",       desc.numStereocenters.toString(),       "chiral centers")
  ))}
</div>"""

        val filterSection = """
<div class="section" id="drug-filters">
  ${sectionTitle("✅", "Drug-Likeness Filters")}
  ${dataTable(
        listOf("Filter", "Rules", "Result", "Detail"),
        listOf(
            listOf("Lipinski Ro5", "MW≤500, LogP≤5, HBD≤5, HBA≤10",
                passFailBadge(desc.lipinskiPasses),
                "MW=${desc.molecularWeight.fmt(0)}, LogP=${desc.logP.fmt(2)}, HBD=${desc.hBondDonors}, HBA=${desc.hBondAcceptors}"),
            listOf("Veber", "TPSA≤140 Å², RotBonds≤10",
                passFailBadge(desc.veberPasses),
                "TPSA=${desc.topologicalPolarSurfaceArea.fmt(1)}, RotB=${desc.numRotatableBonds}"),
            listOf("Ghose", "MW 160–480, LogP −0.4–5.6, MR 40–130, Atoms 20–70",
                passFailBadge(desc.ghosePasses),
                "MW=${desc.molecularWeight.fmt(0)}, LogP=${desc.logP.fmt(2)}, MR=${desc.molarRefractivity.fmt(1)}"),
            listOf("Egan", "LogP −1–6, TPSA≤150 Å²",
                passFailBadge(desc.eganPasses),
                "LogP=${desc.logP.fmt(2)}, TPSA=${desc.topologicalPolarSurfaceArea.fmt(1)}"),
            listOf("Muegge", "9-rule composite (MW, LogP, TPSA, HBD, HBA, rings…)",
                badge("${desc.mueggeScore}/9", if (desc.mueggeScore >= 7) "pass" else if (desc.mueggeScore >= 5) "warn" else "fail"),
                "${desc.mueggeScore} of 9 rules satisfied")
        ),
        listOf("", "", "", "")
    )}
</div>"""

        val admetSection = """
<div class="section" id="drug-admet">
  ${sectionTitle("💊", "ADMET Profile")}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;">
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Absorption & Distribution</div>
      ${dataTable(
            listOf("Property", "Prediction"),
            listOf(
                listOf("Oral Bioavailability", badge(admet.oralBioavailability, if (admet.oralBioavailability == "High") "pass" else if (admet.oralBioavailability == "Low") "fail" else "warn")),
                listOf("Caco-2 Permeability",  badge(admet.caco2Permeability, if (admet.caco2Permeability == "High") "pass" else "fail")),
                listOf("P-gp Substrate",       badge(if (admet.pgpSubstrate) "Yes" else "No", if (admet.pgpSubstrate) "warn" else "pass")),
                listOf("HIA",                  badge(admet.hia, if (admet.hia.startsWith("High")) "pass" else if (admet.hia.startsWith("Low")) "fail" else "warn")),
                listOf("BBB Penetration",      badge(if (admet.bbbPenetration) "Yes" else "No", if (admet.bbbPenetration) "info" else "pass")),
                listOf("VDss",                 "${admet.vdss.fmt(2)} L/kg"),
                listOf("Plasma Protein Binding", badge(admet.plasmaProteinBinding, "info"))
            ),
            listOf("", "")
        )}
    </div>
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Metabolism, Excretion & Toxicity</div>
      ${dataTable(
            listOf("Property", "Prediction"),
            listOf(
                listOf("CYP3A4 Substrate",  badge(if (admet.cyp3a4Substrate) "Yes" else "No", if (admet.cyp3a4Substrate) "warn" else "pass")),
                listOf("CYP2D6 Substrate",  badge(if (admet.cyp2d6Substrate) "Yes" else "No", if (admet.cyp2d6Substrate) "warn" else "pass")),
                listOf("CYP2C9 Substrate",  badge(if (admet.cyp2c9Substrate) "Yes" else "No", if (admet.cyp2c9Substrate) "warn" else "pass")),
                listOf("CYP3A4 Inhibitor",  badge(if (admet.cyp3a4Inhibitor) "Yes" else "No", if (admet.cyp3a4Inhibitor) "fail" else "pass")),
                listOf("CYP2D6 Inhibitor",  badge(if (admet.cyp2d6Inhibitor) "Yes" else "No", if (admet.cyp2d6Inhibitor) "fail" else "pass")),
                listOf("Half-Life",          badge(admet.halfLife, "info")),
                listOf("Renal Clearance",    badge(admet.renalClearance, "info")),
                listOf("hERG Inhibition",    badge(if (admet.hergInhibition) "Risk" else "Low risk", if (admet.hergInhibition) "fail" else "pass")),
                listOf("Ames Mutagenicity",  badge(admet.amesTest, if (admet.amesTest == "Non-mutagenic") "pass" else "fail")),
                listOf("Oral Toxicity",      badge(admet.oralToxicity, if (admet.oralToxicity.contains("Low")) "pass" else if (admet.oralToxicity.contains("High")) "fail" else "warn")),
                listOf("Skin Sensitisation", badge(if (admet.skinSensitization) "Risk" else "Low risk", if (admet.skinSensitization) "warn" else "pass"))
            ),
            listOf("", "")
        )}
    </div>
  </div>
  ${alertBox("warn", "⚠️", "ADMET Disclaimer", "All ADMET values are rule-based in silico estimates. They are not validated clinical predictions. Do not use for clinical or regulatory decisions.")}
</div>"""

        val structSection = """
<div class="section" id="drug-structure">
  ${sectionTitle("🧪", "Structural Analysis")}
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:20px;">
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Atom Inventory</div>
      ${dataTable(
            listOf("Element", "Count", "Aromatic?"),
            mol.atoms.groupBy { it.symbol }
                .entries.sortedByDescending { it.value.size }
                .map { (sym, atoms) ->
                    listOf(sym, atoms.size.toString(), if (atoms.any { it.isAromatic }) badge("Yes","info") else "No")
                },
            listOf("mono", "num", "")
        )}
    </div>
    <div>
      <div style="font-weight:600;color:var(--gray);margin-bottom:10px;font-size:13px;">Bond Inventory</div>
      ${dataTable(
            listOf("Bond Type", "Count", "Order"),
            mol.bonds.groupBy { it.type }
                .entries.sortedBy { it.key.ordinal }
                .map { (type, bonds) ->
                    listOf(type.name, bonds.size.toString(), type.order.fmt(1))
                },
            listOf("", "num", "num")
        )}
      <div style="margin-top:14px;">
      ${dataTable(
            listOf("Ring Property", "Count"),
            listOf(
                listOf("Total Rings", mol.numRings().toString()),
                listOf("Aromatic Rings", mol.numAromaticRings().toString()),
                listOf("Non-aromatic Rings", (mol.numRings() - mol.numAromaticRings()).toString())
            ),
            listOf("", "num")
        )}
      </div>
    </div>
  </div>
</div>"""

        val alertsSection = """
<div class="section" id="drug-alerts">
  ${sectionTitle("🚨", "Structural Alerts")}
  ${if (alerts.isEmpty())
        alertBox("good", "✅", "No Structural Alerts", "No PAINS or reactive group alerts detected across all ${mol.atoms.size} heavy atoms.")
    else alerts.joinToString("") { alert ->
        val type = when (alert.severity) { "High" -> "danger"; "Moderate" -> "warn"; else -> "info" }
        val icon = when (alert.severity) { "High" -> "❌"; "Moderate" -> "⚠️"; else -> "ℹ️" }
        alertBox(type, icon, "[${alert.severity}] ${alert.name}", alert.description)
    }}
</div>"""

        val pharmSection = """
<div class="section" id="drug-pharm">
  ${sectionTitle("🎯", "Pharmacophore Features")}
  <div style="margin-bottom:10px;font-size:13px;color:var(--lgray)">
    ${pharm.features.size} pharmacophore feature(s) identified.
  </div>
  ${dataTable(
        listOf("Feature Type", "Atom Count", "Atom Indices (first 6)"),
        pharm.features.map { f ->
            listOf(
                badge(f.type, when (f.type) { "HBD" -> "pass"; "HBA" -> "info"; "Aromatic" -> "purple"; "Positive" -> "warn"; "Negative" -> "fail"; else -> "info" }),
                f.atomIndices.size.toString(),
                f.atomIndices.take(6).joinToString(", ")
            )
        },
        listOf("", "num", "mono")
    )}
</div>"""

        val appendix = """
<hr class="divider">
<div id="drug-appendix">
  <div class="appendix-header">
    <h2>📎 Detailed Appendix</h2>
    <p>Full descriptor values, fingerprint density, and atom detail table</p>
  </div>

  <div class="section">
    <div class="section-title">A1 — Complete Descriptor Table</div>
    ${dataTable(
        listOf("Descriptor", "Value", "Notes"),
        listOf(
            listOf("Molecular Formula",        desc.molecularFormula,                        "Hill order"),
            listOf("Exact (Monoisotopic) Mass", "${desc.exactMass.fmt(4)} Da",               ""),
            listOf("Average Molecular Weight",  "${desc.molecularWeight.fmt(4)} Da",          ""),
            listOf("LogP (Wildman-Crippen)",    desc.logP.fmt(4),                            "Lipophilicity"),
            listOf("H-Bond Donors",             desc.hBondDonors.toString(),                 "N-H, O-H groups"),
            listOf("H-Bond Acceptors",          desc.hBondAcceptors.toString(),              "N, O atoms"),
            listOf("TPSA (Ertl 2000)",          "${desc.topologicalPolarSurfaceArea.fmt(2)} Å²",""),
            listOf("Rotatable Bonds",           desc.numRotatableBonds.toString(),           "Non-terminal, non-ring"),
            listOf("Rings (total)",             desc.numRings.toString(),                    ""),
            listOf("Aromatic Rings",            desc.numAromaticRings.toString(),            ""),
            listOf("Heavy Atoms",               desc.numHeavyAtoms.toString(),               "Non-H"),
            listOf("Total Atoms (incl. H)",     desc.numAtoms.toString(),                    ""),
            listOf("Bonds",                     desc.numBonds.toString(),                    ""),
            listOf("Stereocenters",             desc.numStereocenters.toString(),            "Approx."),
            listOf("Fraction sp3 Carbons",      desc.fractionCSP3.fmt(4),                   "fsp3"),
            listOf("Molar Refractivity",        desc.molarRefractivity.fmt(3),              "Wildman-Crippen"),
            listOf("Aromatic Atoms",            desc.numAromaticAtoms.toString(),           ""),
            listOf("Lipinski Ro5",              if (desc.lipinskiPasses) "PASS" else "FAIL", ""),
            listOf("Veber Rules",               if (desc.veberPasses) "PASS" else "FAIL",   ""),
            listOf("Ghose Filter",              if (desc.ghosePasses) "PASS" else "FAIL",   ""),
            listOf("Egan Filter",               if (desc.eganPasses) "PASS" else "FAIL",    ""),
            listOf("Muegge Score",              "${desc.mueggeScore}/9",                     "")
        ),
        listOf("", "", "")
    )}
  </div>

  <div class="section">
    <div class="section-title">A2 — Morgan Fingerprint Density</div>
    <div style="font-size:13px;color:var(--lgray);margin-bottom:10px;">
      Bits set: ${fp.count { it }} / ${fp.size} (${(fp.count { it }.toDouble() / fp.size * 100).fmt(1)}% density)
      &nbsp;|&nbsp; MACCS bits set: ${maccs.count { it }} / 166
    </div>
    <div style="font-family:monospace;font-size:10px;color:var(--lgray);line-height:1.6;word-break:break-all;background:var(--xlgray);padding:10px;border-radius:6px;">
      ${fp.take(256).joinToString("") { if (it) "<span style='color:#3b82f6'>1</span>" else "0" }}…
    </div>
  </div>

  <div class="section">
    <div class="section-title">A3 — Atom Detail Table</div>
    ${dataTable(
        listOf("Index", "Symbol", "Atomic #", "Aromatic", "Charge", "Implicit H", "Total H"),
        mol.atoms.map { a ->
            listOf(a.index.toString(), a.symbol, a.atomicNum.toString(),
                if (a.isAromatic) badge("Yes","info") else "No",
                a.charge.toString(), a.implicitH.toString(), a.totalH.toString())
        },
        listOf("num","mono","num","","","num","num")
    )}
  </div>
</div>"""

        val viewer3d = if (viewerPath.isNotEmpty()) molecule3dSection(mol, viewerPath) else ""
        val body = header + toc + summary + viewer3d + descriptorSection + filterSection +
                   admetSection + structSection + alertsSection + pharmSection + appendix
        return page("Drug Report — ${mol.name.ifEmpty { mol.molecularFormula() }}", body)
    }
    private fun dna3dSection(seq: DNASequence, viewerPath: String): String {
        val viewerFile = if (viewerPath.isNotEmpty()) java.io.File(viewerPath).name else "dna_viewer.html"
        val seqUrl = seq.sequence.take(200)
        return """
<div class="section" id="viewer3d">
  ${sectionTitle("🧬", "3D Structure Viewer")}
  <div style="margin-bottom:10px;font-size:13px;color:var(--lgray)">
    Interactive WebGL viewer — drag to rotate, scroll to zoom, touch supported.
  </div>
  <div class="viewer3d-container">
    <iframe src="$viewerFile?seq=$seqUrl" title="DNA 3D Viewer"
            loading="lazy" sandbox="allow-scripts allow-same-origin"></iframe>
    <div class="viewer3d-overlay">
      <a class="viewer3d-btn full" href="$viewerFile?seq=$seqUrl" target="_blank">⛶ Full Screen</a>
    </div>
  </div>
  <div style="margin-top:10px;font-size:12px;color:var(--lgray)">
    B-form double helix &middot;
    A=<span style="color:#f87171">■</span>
    T=<span style="color:#60a5fa">■</span>
    G=<span style="color:#4ade80">■</span>
    C=<span style="color:#facc15">■</span> &middot;
    Blue/amber backbone &middot; hydrogen bonds shown &middot;
    Controls: Helix / Base Pairs / Surface
  </div>
</div>"""
    }


    private fun protein3dSection(seq: ProteinSequence, viewerPath: String): String {
        val viewerFile = if (viewerPath.isNotEmpty()) java.io.File(viewerPath).name else "protein_viewer.html"
        val seqUrl = seq.sequence.take(300)
        return """
<div class="section" id="viewer3d">
  ${sectionTitle("🧬", "3D Structure Viewer")}
  <div style="margin-bottom:10px;font-size:13px;color:var(--lgray)">
    Predicted secondary structure in 3D. Ribbon / Ball-chain / Surface modes.
  </div>
  <div class="viewer3d-container">
    <iframe src="$viewerFile?seq=$seqUrl" title="Protein 3D Viewer"
            loading="lazy" sandbox="allow-scripts allow-same-origin"></iframe>
    <div class="viewer3d-overlay">
      <a class="viewer3d-btn full" href="$viewerFile?seq=$seqUrl" target="_blank">⛶ Full Screen</a>
    </div>
  </div>
  <div style="margin-top:10px;font-size:12px;color:var(--lgray)">
    α-Helix=<span style="color:#f87171">■</span>
    β-Sheet=<span style="color:#60a5fa">■</span>
    Turn=<span style="color:#facc15">■</span>
    Coil=<span style="color:#4ade80">■</span> &middot;
    Hydrophobicity coloring available &middot;
    Chou-Fasman secondary structure prediction
  </div>
</div>"""
    }


    private fun molecule3dSection(mol: Molecule, viewerPath: String): String {
        val viewerFile = if (viewerPath.isNotEmpty()) java.io.File(viewerPath).name else "molecule_viewer.html"
        val smilesEnc = java.net.URLEncoder.encode(mol.smiles, "UTF-8")
        val nameEnc   = java.net.URLEncoder.encode(mol.name.ifEmpty { mol.id.ifEmpty { "Molecule" } }, "UTF-8")
        return """
<div class="section" id="viewer3d">
  ${sectionTitle("🧪", "3D Structure Viewer")}
  <div style="margin-bottom:10px;font-size:13px;color:var(--lgray)">
    Force-directed 3D embedding from SMILES. Ball-and-stick / Space-fill / Wireframe.
  </div>
  <div class="viewer3d-container">
    <iframe src="$viewerFile?smiles=$smilesEnc&name=$nameEnc" title="Molecule 3D Viewer"
            loading="lazy" sandbox="allow-scripts allow-same-origin"></iframe>
    <div class="viewer3d-overlay">
      <a class="viewer3d-btn full" href="$viewerFile?smiles=$smilesEnc&name=$nameEnc" target="_blank">⛶ Full Screen</a>
    </div>
  </div>
  <div style="margin-top:10px;font-size:12px;color:var(--lgray)">
    CPK: C=<span style="color:#94a3b8">■</span>
    O=<span style="color:#ef4444">■</span>
    N=<span style="color:#60a5fa">■</span>
    S=<span style="color:#facc15">■</span>
    Aromatic=<span style="color:#4ade80">■</span> &middot;
    300-iteration force-directed layout
  </div>
</div>"""
    }


}
