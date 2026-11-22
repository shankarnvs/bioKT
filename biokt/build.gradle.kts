plugins {
    kotlin("jvm") version "1.9.22"
    application
}

group = "com.biokt"
version = "2.0.0"

repositories {
    mavenCentral()
}

dependencies {
    implementation(kotlin("stdlib"))
    testImplementation(kotlin("test"))
    testImplementation("org.junit.jupiter:junit-jupiter:5.10.0")
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

application {
    mainClass.set("biokt.MainKt")
}
