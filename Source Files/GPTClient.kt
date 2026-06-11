package biokt

import kotlin.math.*

// ============================================================
// GPT CLIENT — LLM-powered biological interpretation
//
// Supported providers:
//   • OpenAI   — GPT-4o, GPT-4-turbo, GPT-3.5-turbo
//   • Anthropic — Claude-3.5-Sonnet, Claude-3-Haiku, Claude-3-Opus
//   • Google   — Gemini 1.5 Pro, Gemini 1.5 Flash
//   • Ollama   — Any local model (llama3, mistral, etc.)
//
// Usage:
//   val client = GPTClient.openai(apiKey = "sk-...")
//   val interp = client.interpretDNA(mySeq, mlPredictions)
// ============================================================

// ─────────────────────────────────────────────────────────────
// PROVIDER ENUM
// ─────────────────────────────────────────────────────────────

enum class LLMProvider {
    OPENAI, ANTHROPIC, GOOGLE, OLLAMA
}

// ─────────────────────────────────────────────────────────────
// MODEL CATALOG
// ─────────────────────────────────────────────────────────────

object LLMModels {
    object OpenAI {
        const val GPT4O          = "gpt-4o"
        const val GPT4_TURBO     = "gpt-4-turbo"
        const val GPT4           = "gpt-4"
        const val GPT35_TURBO    = "gpt-3.5-turbo"
        val ALL = listOf(GPT4O, GPT4_TURBO, GPT4, GPT35_TURBO)
    }
    object Anthropic {
        const val CLAUDE_35_SONNET = "claude-sonnet-4-5"
        const val CLAUDE_3_OPUS    = "claude-opus-4-5"
        const val CLAUDE_3_HAIKU   = "claude-haiku-4-5"
        val ALL = listOf(CLAUDE_35_SONNET, CLAUDE_3_OPUS, CLAUDE_3_HAIKU)
    }
    object Google {
        const val GEMINI_15_PRO   = "gemini-1.5-pro"
        const val GEMINI_15_FLASH = "gemini-1.5-flash"
        const val GEMINI_10_PRO   = "gemini-1.0-pro"
        val ALL = listOf(GEMINI_15_PRO, GEMINI_15_FLASH, GEMINI_10_PRO)
    }
    object Ollama {
        const val LLAMA3     = "llama3"
        const val MISTRAL    = "mistral"
        const val CODELLAMA  = "codellama"
        const val MEDLLAMA   = "medllama2"
        val ALL = listOf(LLAMA3, MISTRAL, CODELLAMA, MEDLLAMA)
    }
}

// ─────────────────────────────────────────────────────────────
// LLM RESPONSE
// ─────────────────────────────────────────────────────────────

data class LLMResponse(
    val text: String,
    val model: String,
    val provider: LLMProvider,
    val promptTokens: Int,
    val completionTokens: Int,
    val success: Boolean,
    val error: String? = null
) {
    val totalTokens: Int get() = promptTokens + completionTokens
    fun isError(): Boolean = !success || error != null
}

// ─────────────────────────────────────────────────────────────
// GPT CLIENT
// ─────────────────────────────────────────────────────────────

class GPTClient(
    val provider: LLMProvider,
    val model: String,
    private val apiKey: String,
    private val baseUrl: String = defaultBaseUrl(provider),
    val systemPrompt: String = DEFAULT_SYSTEM_PROMPT,
    val maxTokens: Int = 1200,
    val temperature: Double = 0.3
) {

    companion object {

        const val DEFAULT_SYSTEM_PROMPT = """You are BioKt-GPT, an expert bioinformatics assistant integrated 
into the BioKt v2.0 Kotlin library. You provide concise, scientifically accurate interpretations of 
biological sequences, protein structures, and drug molecules. 
- Use precise scientific terminology
- Keep responses focused and under 400 words unless asked for detail
- Always acknowledge uncertainty in predictions
- Format key findings as bullet points when appropriate"""

        fun defaultBaseUrl(provider: LLMProvider) = when (provider) {
            LLMProvider.OPENAI    -> "https://api.openai.com/v1/chat/completions"
            LLMProvider.ANTHROPIC -> "https://api.anthropic.com/v1/messages"
            LLMProvider.GOOGLE    -> "https://generativelanguage.googleapis.com/v1beta/models"
            LLMProvider.OLLAMA    -> "http://localhost:11434/api/chat"
        }

        /** Factory: OpenAI */
        fun openai(
            apiKey: String,
            model: String = LLMModels.OpenAI.GPT4O,
            maxTokens: Int = 1200,
            temperature: Double = 0.3
        ) = GPTClient(LLMProvider.OPENAI, model, apiKey,
            maxTokens = maxTokens, temperature = temperature)

        /** Factory: Anthropic Claude */
        fun anthropic(
            apiKey: String,
            model: String = LLMModels.Anthropic.CLAUDE_35_SONNET,
            maxTokens: Int = 1200,
            temperature: Double = 0.3
        ) = GPTClient(LLMProvider.ANTHROPIC, model, apiKey,
            maxTokens = maxTokens, temperature = temperature)

        /** Factory: Google Gemini */
        fun google(
            apiKey: String,
            model: String = LLMModels.Google.GEMINI_15_PRO,
            maxTokens: Int = 1200,
            temperature: Double = 0.3
        ) = GPTClient(LLMProvider.GOOGLE, model, apiKey,
            maxTokens = maxTokens, temperature = temperature)

        /** Factory: Ollama local model */
        fun ollama(
            model: String = LLMModels.Ollama.LLAMA3,
            host: String = "localhost",
            port: Int = 11434,
            maxTokens: Int = 1200,
            temperature: Double = 0.3
        ) = GPTClient(LLMProvider.OLLAMA, model, apiKey = "",
            baseUrl = "http://$host:$port/api/chat",
            maxTokens = maxTokens, temperature = temperature)
    }

    // ── Core HTTP call ────────────────────────────────────────

    fun chat(userMessage: String, contextSystemPrompt: String = ""): LLMResponse {
        val effectiveSystem = if (contextSystemPrompt.isNotEmpty())
            "$systemPrompt\n\n$contextSystemPrompt" else systemPrompt
        return when (provider) {
            LLMProvider.OPENAI    -> callOpenAI(effectiveSystem, userMessage)
            LLMProvider.ANTHROPIC -> callAnthropic(effectiveSystem, userMessage)
            LLMProvider.GOOGLE    -> callGoogle(effectiveSystem, userMessage)
            LLMProvider.OLLAMA    -> callOllama(effectiveSystem, userMessage)
        }
    }

    private fun callOpenAI(system: String, user: String): LLMResponse {
        val body = """
{
  "model": "$model",
  "max_tokens": $maxTokens,
  "temperature": $temperature,
  "messages": [
    {"role": "system", "content": ${jsonString(system)}},
    {"role": "user",   "content": ${jsonString(user)}}
  ]
}"""
        return try {
            val raw = httpPost(baseUrl, body, mapOf(
                "Authorization" to "Bearer $apiKey",
                "Content-Type"  to "application/json"
            ))
            parseOpenAIResponse(raw)
        } catch (e: Exception) {
            LLMResponse("", model, provider, 0, 0, false, e.message)
        }
    }

    private fun callAnthropic(system: String, user: String): LLMResponse {
        val body = """
{
  "model": "$model",
  "max_tokens": $maxTokens,
  "system": ${jsonString(system)},
  "messages": [
    {"role": "user", "content": ${jsonString(user)}}
  ]
}"""
        return try {
            val raw = httpPost(baseUrl, body, mapOf(
                "x-api-key"         to apiKey,
                "anthropic-version" to "2023-06-01",
                "Content-Type"      to "application/json"
            ))
            parseAnthropicResponse(raw)
        } catch (e: Exception) {
            LLMResponse("", model, provider, 0, 0, false, e.message)
        }
    }

    private fun callGoogle(system: String, user: String): LLMResponse {
        val url  = "$baseUrl/$model:generateContent?key=$apiKey"
        val body = """
{
  "system_instruction": {"parts": [{"text": ${jsonString(system)}}]},
  "contents": [{"parts": [{"text": ${jsonString(user)}}]}],
  "generationConfig": {
    "maxOutputTokens": $maxTokens,
    "temperature": $temperature
  }
}"""
        return try {
            val raw = httpPost(url, body, mapOf("Content-Type" to "application/json"))
            parseGoogleResponse(raw)
        } catch (e: Exception) {
            LLMResponse("", model, provider, 0, 0, false, e.message)
        }
    }

    private fun callOllama(system: String, user: String): LLMResponse {
        val body = """
{
  "model": "$model",
  "stream": false,
  "options": {"temperature": $temperature, "num_predict": $maxTokens},
  "messages": [
    {"role": "system", "content": ${jsonString(system)}},
    {"role": "user",   "content": ${jsonString(user)}}
  ]
}"""
        return try {
            val raw = httpPost(baseUrl, body, mapOf("Content-Type" to "application/json"))
            parseOllamaResponse(raw)
        } catch (e: Exception) {
            LLMResponse("", model, provider, 0, 0, false, e.message)
        }
    }

    // ── HTTP helper (pure JVM, no external libs) ──────────────

    private fun httpPost(url: String, body: String, headers: Map<String, String>): String {
        val conn = java.net.URL(url).openConnection() as java.net.HttpURLConnection
        conn.requestMethod = "POST"
        conn.doOutput      = true
        conn.connectTimeout = 30_000
        conn.readTimeout    = 60_000
        headers.forEach { (k, v) -> conn.setRequestProperty(k, v) }
        conn.outputStream.use { it.write(body.toByteArray(Charsets.UTF_8)) }
        val code = conn.responseCode
        val stream = if (code in 200..299) conn.inputStream else conn.errorStream
        val response = stream.use { it.readBytes().toString(Charsets.UTF_8) }
        if (code !in 200..299) throw RuntimeException("HTTP $code: $response")
        return response
    }

    // ── Response parsers ──────────────────────────────────────

    private fun parseOpenAIResponse(raw: String): LLMResponse {
        val text    = extractJsonValue(raw, "content")
                   ?: extractJsonValue(raw, "message.content")
                   ?: ""
        val pTok    = extractJsonInt(raw, "prompt_tokens") ?: 0
        val cTok    = extractJsonInt(raw, "completion_tokens") ?: 0
        return LLMResponse(text.trim(), model, provider, pTok, cTok, text.isNotEmpty())
    }

    private fun parseAnthropicResponse(raw: String): LLMResponse {
        val text = extractJsonValue(raw, "text") ?: ""
        val pTok = extractJsonInt(raw, "input_tokens")  ?: 0
        val cTok = extractJsonInt(raw, "output_tokens") ?: 0
        return LLMResponse(text.trim(), model, provider, pTok, cTok, text.isNotEmpty())
    }

    private fun parseGoogleResponse(raw: String): LLMResponse {
        val text = extractJsonValue(raw, "text") ?: ""
        return LLMResponse(text.trim(), model, provider, 0, 0, text.isNotEmpty())
    }

    private fun parseOllamaResponse(raw: String): LLMResponse {
        val text = extractJsonValue(raw, "content") ?: ""
        return LLMResponse(text.trim(), model, provider, 0, 0, text.isNotEmpty())
    }

    // Regex-based JSON extractors (avoids full JSON parser dependency)
    private fun extractJsonValue(json: String, key: String): String? {
        val patterns = listOf(
            """"$key"\s*:\s*"((?:[^"\\]|\\.)*)"""".toRegex(),
            """"text"\s*:\s*"((?:[^"\\]|\\.)*)"""".toRegex()
        )
        val lastKey = key.substringAfterLast(".")
        for (pat in patterns.take(if (key.contains(".")) 2 else 1)) {
            val m = pat.find(json) ?: continue
            return m.groupValues[1]
                .replace("\\n", "\n")
                .replace("\\\"", "\"")
                .replace("\\\\", "\\")
                .replace("\\t", "\t")
        }
        return null
    }

    private fun extractJsonInt(json: String, key: String): Int? =
        """"$key"\s*:\s*(\d+)""".toRegex().find(json)?.groupValues?.get(1)?.toIntOrNull()

    private fun jsonString(s: String): String =
        "\"" + s.replace("\\", "\\\\").replace("\"", "\\\"")
            .replace("\n", "\\n").replace("\r", "").replace("\t", "\\t") + "\""

    // ─────────────────────────────────────────────────────────
    // HIGH-LEVEL BIOLOGICAL INTERPRETATION METHODS
    // ─────────────────────────────────────────────────────────

    /** Interpret a DNA sequence with optional ML predictions */
    fun interpretDNA(
        seq: DNASequence,
        mlPredictions: Map<String, MLPrediction> = emptyMap(),
        includeAnalysis: Boolean = true
    ): LLMResponse {
        val context = buildString {
            append("DNA Sequence Analysis Request\n\n")
            append("Sequence ID: ${seq.id.ifEmpty { "unknown" }}\n")
            append("Length: ${seq.length} bp\n")
            append("GC Content: ${"%.1f".format(seq.gcContent())}%\n")
            append("Melting Temperature: ${"%.1f".format(seq.meltingTemperature())}°C\n")

            if (includeAnalysis) {
                val orfs = seq.findOrfs(minLength = 30)
                append("ORFs found: ${orfs.size}\n")
                if (orfs.isNotEmpty()) {
                    append("Largest ORF: ${orfs.maxBy { it.length }?.let { "${it.length} bp, protein: ${it.protein.sequence.take(30)}" } ?: "—"}\n")
                }
                val sites = RestrictionEnzymes.digest(seq).take(5)
                if (sites.isNotEmpty()) {
                    append("Restriction sites: ${sites.joinToString(", ") { it.enzyme }}\n")
                }
                val complexity = SeqStats.linguisticComplexity(seq)
                append("Sequence complexity: ${"%.2f".format(complexity)}\n")
            }

            if (mlPredictions.isNotEmpty()) {
                append("\nML Predictions:\n")
                mlPredictions.forEach { (task, pred) ->
                    append("  ${task.replace("_", " ").let { s -> if (s.isEmpty()) s else s[0].toUpperCase() + s.drop(1) }}: ${pred.label} (${"%.1f".format(pred.confidence * 100)}% confidence)\n")
                }
            }
        }

        val prompt = """$context

Please provide a biological interpretation of this DNA sequence covering:
1. What type of sequence this likely is (coding, regulatory, non-coding, etc.)
2. Key features and their biological significance
3. Any notable characteristics (GC content implications, structural features)
4. Confidence in the ML predictions if provided
5. Suggested follow-up analyses"""

        return chat(prompt)
    }

    /** Interpret a protein sequence with optional ML predictions */
    fun interpretProtein(
        seq: ProteinSequence,
        mlPredictions: Map<String, MLPrediction> = emptyMap()
    ): LLMResponse {
        val context = buildString {
            append("Protein Sequence Analysis Request\n\n")
            append("Protein ID: ${seq.id.ifEmpty { "unknown" }}\n")
            append("Length: ${seq.length} aa\n")
            append("MW: ${"%.1f".format(seq.molecularWeight())} Da\n")
            append("pI: ${"%.2f".format(seq.isoelectricPoint())}\n")
            append("GRAVY: ${"%.3f".format(ProteinDescriptors.gravyIndex(seq))}\n")
            append("Aliphatic index: ${"%.1f".format(ProteinDescriptors.aliphaticIndex(seq))}\n")
            append("Instability index: ${"%.1f".format(seq.instabilityIndex())} (${if (seq.instabilityIndex() < 40) "STABLE" else "UNSTABLE"})\n")
            append("Boman index: ${"%.3f".format(ProteinDescriptors.bomanIndex(seq))}\n")
            append("Net charge at pH 7: ${"%.2f".format(ProteinDescriptors.chargeAtPH(seq, 7.0))}\n")

            val tmHelices = ProteinInteraction.predictTransmembraneHelices(seq)
            append("Predicted TM helices: ${tmHelices.size}\n")
            val coiledCoils = ProteinInteraction.predictCoiledCoil(seq)
            append("Predicted coiled-coil segments: ${coiledCoils.size}\n")

            append("Sequence (first 60 aa): ${seq.sequence.take(60)}${if (seq.length > 60) "…" else ""}\n")

            if (mlPredictions.isNotEmpty()) {
                append("\nML Predictions:\n")
                mlPredictions.forEach { (task, pred) ->
                    append("  ${task.replace("_", " ").let { s -> if (s.isEmpty()) s else s[0].toUpperCase() + s.drop(1) }}: ${pred.label} (${"%.1f".format(pred.confidence * 100)}%)\n")
                }
            }
        }

        val prompt = """$context

Please interpret this protein sequence covering:
1. Predicted function and protein family based on physicochemical properties
2. Cellular localization clues (membrane, secreted, nuclear, cytoplasmic)
3. Stability, expression, and purification considerations
4. Notable features (hydrophobicity, charge distribution, aromatic residues)
5. Clinical or industrial relevance if apparent"""

        return chat(prompt)
    }

    /** Interpret a drug molecule with ML predictions and ADMET data */
    fun interpretDrug(
        mol: Molecule,
        mlPredictions: Map<String, MLPrediction> = emptyMap(),
        includeAdmet: Boolean = true
    ): LLMResponse {
        val desc = MolDescriptors.calculate(mol)
        val context = buildString {
            append("Drug / Molecule Analysis Request\n\n")
            append("Name: ${mol.name.ifEmpty { mol.id.ifEmpty { "Unknown" } }}\n")
            append("SMILES: ${mol.smiles}\n")
            append("Formula: ${mol.molecularFormula()}\n")
            append("MW: ${"%.2f".format(desc.molecularWeight)} Da\n")
            append("LogP: ${"%.2f".format(desc.logP)}\n")
            append("HBD/HBA: ${desc.hBondDonors}/${desc.hBondAcceptors}\n")
            append("TPSA: ${"%.1f".format(desc.topologicalPolarSurfaceArea)} Å²\n")
            append("Rotatable bonds: ${desc.numRotatableBonds}\n")
            append("Rings: ${desc.numRings} (aromatic: ${desc.numAromaticRings})\n")
            append("Lipinski Ro5: ${if (desc.lipinskiPasses) "PASS" else "FAIL"}\n")
            append("Veber rules: ${if (desc.veberPasses) "PASS" else "FAIL"}\n")
            append("Muegge score: ${desc.mueggeScore}/9\n")

            if (includeAdmet) {
                append("\nADMET Profile:\n")
                append("  Oral bioavailability: ${desc.admet.oralBioavailability}\n")
                append("  BBB penetration: ${if (desc.admet.bbbPenetration) "Yes" else "No"}\n")
                append("  CYP3A4 substrate: ${if (desc.admet.cyp3a4Substrate) "Yes" else "No"}\n")
                append("  hERG inhibition: ${if (desc.admet.hergInhibition) "Risk" else "Low risk"}\n")
                append("  Ames mutagenicity: ${desc.admet.amesTest}\n")
                append("  Oral toxicity: ${desc.admet.oralToxicity}\n")
                append("  Half-life: ${desc.admet.halfLife}\n")
            }

            val alerts = MolDescriptors.checkStructuralAlerts(mol)
            if (alerts.isNotEmpty()) {
                append("\nStructural Alerts:\n")
                alerts.forEach { a -> append("  [${a.severity}] ${a.name}: ${a.description}\n") }
            }

            if (mlPredictions.isNotEmpty()) {
                append("\nML Predictions:\n")
                mlPredictions.forEach { (task, pred) ->
                    append("  ${task.replace("_", " ").let { s -> if (s.isEmpty()) s else s[0].toUpperCase() + s.drop(1) }}: ${pred.label} (${"%.1f".format(pred.confidence * 100)}%)\n")
                }
            }
        }

        val prompt = """$context

Please provide a medicinal chemistry interpretation covering:
1. Drug-likeness assessment and oral bioavailability potential
2. Key pharmacophore features and likely mechanism of action clues
3. Safety/toxicity concerns from the structural alerts and ADMET profile
4. Drug metabolism considerations (CYP interactions, half-life)
5. Optimisation suggestions if drug-likeness concerns are present
6. Therapeutic area potential based on structural features"""

        return chat(prompt)
    }

    /** Ask a free-form bioinformatics question */
    fun ask(question: String): LLMResponse = chat(question)

    /** Compare two sequences or molecules */
    fun compare(description1: String, description2: String, task: String): LLMResponse {
        val prompt = """Compare the following two biological entities in the context of $task:

Entity 1:
$description1

Entity 2:
$description2

Please provide:
1. Key similarities
2. Key differences
3. Functional implications of the differences
4. Which is likely more suitable for $task and why"""
        return chat(prompt)
    }

    /** Batch interpretation — process multiple items and return all responses */
    fun batchInterpret(items: List<Pair<String, String>>): List<LLMResponse> {
        return items.map { (context, question) ->
            try {
                chat("$context\n\n$question")
            } catch (e: Exception) {
                LLMResponse("", model, provider, 0, 0, false, e.message)
            }
        }
    }

    // ─────────────────────────────────────────────────────────
    // STREAMING (for Kotlin coroutines / async contexts)
    // Returns the full response but processes chunk by chunk.
    // In a real coroutine context, replace runBlocking with
    // a proper suspend function.
    // ─────────────────────────────────────────────────────────

    fun chatWithCallback(
        userMessage: String,
        onChunk: (String) -> Unit,
        onComplete: (LLMResponse) -> Unit
    ) {
        // Synchronous fallback — call normally and deliver as one chunk
        val resp = chat(userMessage)
        if (!resp.isError()) {
            resp.text.chunked(80).forEach { chunk ->
                onChunk(chunk)
                Thread.sleep(10) // simulate streaming
            }
        }
        onComplete(resp)
    }
}

// ─────────────────────────────────────────────────────────────
// CONVENIENCE FUNCTIONS
// ─────────────────────────────────────────────────────────────

/** Quick one-liner: interpret DNA with GPT */
fun DNASequence.interpretWith(client: GPTClient, mlPredictions: Map<String, MLPrediction> = emptyMap()): LLMResponse =
    client.interpretDNA(this, mlPredictions)

/** Quick one-liner: interpret protein with GPT */
fun ProteinSequence.interpretWith(client: GPTClient, mlPredictions: Map<String, MLPrediction> = emptyMap()): LLMResponse =
    client.interpretProtein(this, mlPredictions)

/** Quick one-liner: interpret molecule with GPT */
fun Molecule.interpretWith(client: GPTClient, mlPredictions: Map<String, MLPrediction> = emptyMap()): LLMResponse =
    client.interpretDrug(this, mlPredictions)

// ─────────────────────────────────────────────────────────────
// GPT CONFIG BUILDER (fluent API)
// ─────────────────────────────────────────────────────────────

class GPTConfig private constructor() {
    var provider: LLMProvider = LLMProvider.OPENAI
    var model: String = LLMModels.OpenAI.GPT4O
    var apiKey: String = ""
    var maxTokens: Int = 1200
    var temperature: Double = 0.3
    var systemPrompt: String = GPTClient.DEFAULT_SYSTEM_PROMPT
    var ollamaHost: String = "localhost"
    var ollamaPort: Int = 11434

    fun build(): GPTClient = when (provider) {
        LLMProvider.OPENAI    -> GPTClient.openai(apiKey, model, maxTokens, temperature)
        LLMProvider.ANTHROPIC -> GPTClient.anthropic(apiKey, model, maxTokens, temperature)
        LLMProvider.GOOGLE    -> GPTClient.google(apiKey, model, maxTokens, temperature)
        LLMProvider.OLLAMA    -> GPTClient.ollama(model, ollamaHost, ollamaPort, maxTokens, temperature)
    }

    companion object {
        fun create(block: GPTConfig.() -> Unit): GPTClient {
            return GPTConfig().apply(block).build()
        }
    }
}

// ─────────────────────────────────────────────────────────────
// USAGE EXAMPLE (in comments — no side effects at load time)
// ─────────────────────────────────────────────────────────────
/*
  // ── OpenAI ────────────────────────────────────────────────
  val gpt = GPTClient.openai(apiKey = System.getenv("OPENAI_API_KEY"))
  val dna = DNASequence("ATGGCCATTGTAATGGGCCGCTGA")
  val ml  = SequenceMLPredictor.fullAnalysis(dna)
  val resp = gpt.interpretDNA(dna, ml)
  println(resp.text)

  // ── Anthropic Claude ──────────────────────────────────────
  val claude = GPTClient.anthropic(apiKey = System.getenv("ANTHROPIC_API_KEY"),
      model = LLMModels.Anthropic.CLAUDE_35_SONNET)
  val mol  = Molecule("CC(=O)Oc1ccccc1C(=O)O", name="Aspirin")
  val mlD  = DrugMLPredictor.fullAnalysis(mol)
  println(claude.interpretDrug(mol, mlD).text)

  // ── Google Gemini ─────────────────────────────────────────
  val gemini = GPTClient.google(apiKey = System.getenv("GEMINI_API_KEY"))
  val prot   = ProteinSequence("MKALVLLYLLFSSAYS", id="sp1")
  println(gemini.interpretProtein(prot).text)

  // ── Local Ollama (no API key needed) ─────────────────────
  val ollama = GPTClient.ollama(model = "llama3")
  println(ollama.ask("What does GRAVY index tell us about a protein?").text)

  // ── Fluent builder API ────────────────────────────────────
  val client = GPTConfig.create {
      provider    = LLMProvider.ANTHROPIC
      model       = LLMModels.Anthropic.CLAUDE_35_SONNET
      apiKey      = System.getenv("ANTHROPIC_API_KEY")
      maxTokens   = 800
      temperature = 0.2
  }

  // ── Extension function style ─────────────────────────────
  val interp = dna.interpretWith(gpt, ml)
  println(interp.text)
*/
