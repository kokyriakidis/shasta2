// Shasta.
#include "AssemblerOptions.hpp"
using namespace shasta;

// Standard library.
#include "stdexcept.hpp"



shasta::AssemblerOptions::AssemblerOptions(int argc, char** argv) :
    CLI::App("Shasta2. Under development.")
{
    allow_config_extras(false);
    set_config("--config", "", "Configuration file.");

    get_formatter()->column_width(20);

    addOptions();

    try {
        parse(argc, argv);
    } catch(const CLI::ParseError& e) {
         exit(e);
         if(e.get_name() == "CallForHelp") {
             isHelp = true;
         } else {
             throw runtime_error("Error parsing options.");
         }
    }
}



void AssemblerOptions::addOptions()
{
    add_option("--input", inputFileNames,
        "Input fasta or fastq files (uncompressed)."
        )->capture_default_str();

    add_option("--output", assemblyDirectory,
        "Assembly directory for output. Must not exist."
        )->capture_default_str();

    add_option("--command", command,
        "Shasta2 command.")
    ->capture_default_str();

    add_option("--memory-mode", memoryMode,
        "Memory mode. Specify disk, 4K, or 2M.")
    ->capture_default_str();

    add_option("--memory-backing", memoryBacking,
        "Memory backing.Specify anonymous or filesystem."
        )->capture_default_str();

    add_option("--threads", threadCount,
        "Number of threads (0 to use hardware_concurrency)."
        )->capture_default_str();

    add_option("--explore-access", exploreAccess,
        "Specifies accessibility of Shasta2 http server.\n"
        "Specify user, local, or unrestricted."
        )->capture_default_str();

    add_option("--port", port,
        "Port number for http server."
        )->capture_default_str();

    add_option("--min-read-length", minReadLength,
        "Read length cutoff."
        )->capture_default_str();

    add_option("--k", k,
        "Marker length"
        )->capture_default_str();

    add_option("--marker-density", markerDensity,
        "Marker density."
        )->capture_default_str();

    add_option("--min-anchor-coverage", minAnchorCoverage,
        "Minimum anchor coverage"
        )->capture_default_str();

    add_option("--max-anchor-coverage", maxAnchorCoverage,
        "Maximum anchor coverage"
        )->capture_default_str();

    add_option("--local-assembly-estimated-offset-ratio",
        localAssemblyOptions.estimatedOffsetRatio,
        "Estimated offset ratio for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-vertex-sampling-rate", localAssemblyOptions.vertexSamplingRate,
        "Vertex sampling rate for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-match-score", localAssemblyOptions.matchScore,
        "Match score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-mismatch-score", localAssemblyOptions.mismatchScore,
        "Mismatch score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-gap-score", localAssemblyOptions.gapScore,
        "Gap score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-skip-bases", localAssemblyOptions.maxSkipBases,
        "Maximum skip in bases for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-drift-ratio", localAssemblyOptions.maxDrift,
        "Maximum drift ratio for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-min-half-band", localAssemblyOptions.minHalfBand,
        "Minimum half band, in markers, for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-min-score-ratio", localAssemblyOptions.minScoreRatio,
        "Score ratio threshold for discarding alignments\n"
        "for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-msa-length", localAssemblyOptions.maxMsaLength,
        "Maximum length of a multiple sequence alignment\n"
        "for local assembly."
        )->capture_default_str();
}



void AssemblerOptions::write(ostream& s) const
{
    s << config_to_str(true,true);
}
