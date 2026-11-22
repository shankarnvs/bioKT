plugins {
    kotlin("jvm") version "1.9.22"
    application
}

group = "com.biokt"
version = "2.0.0"

repositories {
    mavenCentral()
}

val jmh by sourceSets.creating {
    val main = sourceSets.main.get()
    compileClasspath += main.output
    runtimeClasspath += main.output
}

configurations[jmh.implementationConfigurationName].extendsFrom(configurations.implementation.get())

dependencies {
    implementation(kotlin("stdlib"))
    testImplementation(kotlin("test"))
    testImplementation("org.junit.jupiter:junit-jupiter:5.10.0")
    add(jmh.implementationConfigurationName, "org.openjdk.jmh:jmh-core:1.37")
    add(jmh.annotationProcessorConfigurationName, "org.openjdk.jmh:jmh-generator-annprocess:1.37")
}

kotlin {
    jvmToolchain(17)
    sourceSets {
        named("main") {
            kotlin.setSrcDirs(listOf("."))
            kotlin.include("*.kt", "rice/*.kt")
            kotlin.exclude("BioKtTest.kt")
        }
        named("test") {
            kotlin.setSrcDirs(listOf("."))
            kotlin.include("BioKtTest.kt", "src/test/kotlin/**/*.kt")
        }
    }
}

sourceSets.main {
    resources.srcDir(".")
    resources.include("*_viewer.html")
}

tasks.test {
    useJUnitPlatform()
}

tasks.register<JavaExec>("jmh") {
    group = "verification"
    description = "Run the JVM performance benchmarks"
    classpath = jmh.runtimeClasspath
    javaLauncher.set(javaToolchains.launcherFor { languageVersion.set(JavaLanguageVersion.of(17)) })
    mainClass.set("org.openjdk.jmh.Main")
}

tasks.register<JavaExec>("molecularBindingCheck") {
    group = "verification"
    description = "Compare the albumin pilot scores with the published binding measurements"
    classpath = jmh.runtimeClasspath
    javaLauncher.set(javaToolchains.launcherFor { languageVersion.set(JavaLanguageVersion.of(17)) })
    mainClass.set("biokt.benchmark.AlbuminBindingBenchmark")
}

tasks.register<JavaExec>("interactionChecks") {
    group = "verification"
    description = "Run the interaction benchmark fixture checks"
    classpath = jmh.runtimeClasspath
    javaLauncher.set(javaToolchains.launcherFor { languageVersion.set(JavaLanguageVersion.of(17)) })
    mainClass.set("biokt.benchmark.InteractionBenchmarks")
}

application {
    mainClass.set("biokt.MainKt")
}
