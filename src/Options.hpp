#pragma once

// Shasta.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

// CLI11. Requires package libcli11-dev.
#include "CLI/CLI.hpp"

namespace shasta2 {
    class Options;
}


class shasta2::Options : public CLI::App {
public:

    // Constructor from command line options.
    Options(int argumentCount, char** arguments);

    // Constructor from a configuration file.
    Options(const string& fileName = "");

    // Write a configuration file.
    void write(ostream&) const;

    bool isHelp = false;

    // The names of the input fasta/fastq files.
    vector<string> inputFileNames;

    // An anchor is not generated if its sequence contains an exact repeat
    // consisting of n copies of a unit of length (period) p, if
    // n > maxAnchorRepeatLength[p-1].
    // So for example:
    // - maxAnchorRepeatLength[0] is the maximum allowed length of
    //   a homopolymer run.
    // - maxAnchorRepeatLength[1] is the maximum allowed number of
    //   copies of a repeat with period 2 (e. g. ATATAT).
    //   Note this is the number of copies, not number of bases.
    //   So if maxAnchorRepeatLength[1] is 3, the anchor is not
    //   generated if it contains a 2-base repeat with more than 3 copies
    //   (a total 6 bases).
    vector<uint64_t> maxAnchorRepeatLength = {6, 4, 4, 4, 4};



	// Options with a simple type are defined in OptionsDefine.hpp
	#define SHASTA2_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
		type name = defaultValue;
	#include "OptionsDefine.hpp"
	#undef SHASTA2_OPTION_DEFINE

private:
    void addOptions();
};


