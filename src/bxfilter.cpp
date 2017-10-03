#include "bxfilter.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

namespace opt {

    static std::string bam; // the bam to analyze
    static bool verbose = false;
    static int mapping_quality = 0;
}

static const char* shortopts = "a:hv";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 },
        {"mapping_quality", required_argument, NULL, 'a'}
};


static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools filter in.bam -q X > out.bam\n"
                "Description: Extract all reads from in.bam that satisfy parameters \n"
                "\n"
                "  General options\n"
                "-q, --mapping quality                    Filter read pairs with mapping quality of any read lower than [a]\n"
                "  -v, --verbose                        Set verbose output\n"
                "\n";

static void parseOptions(int argc, char** argv);
static bool CheckConditions(SeqLib::BamRecord &r1, SeqLib::BamRecord &r2);


void runFilter(int argc, char** argv) {
    parseOptions(argc, argv);

    // open the BAM
    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }

    SeqLib::BamWriter writer;
    writer.Open("-");
    writer.SetHeader(reader.Header());
    writer.WriteHeader();
    // loop and filter
    SeqLib::BamRecord r1;
    SeqLib::BamRecord r2;
    size_t count = 0;
    while (reader.GetNextRecord(r1)) {
        if (reader.GetNextRecord(r2)) {
            if (CheckConditions(r1, r2)) {
                writer.WriteRecord(r1);
                writer.WriteRecord(r2);
            }
        }
    }
}

static bool CheckConditions(SeqLib::BamRecord &r1, SeqLib::BamRecord &r2) {
    if (r1.MapQuality() < opt::mapping_quality || r2.MapQuality() < opt::mapping_quality) {
        return false;
    }
    return true;
}

static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc < 1)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
    }

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'a': arg >> opt::mapping_quality; break;
            case 'v': opt::verbose = true; break;
        }
    }

    if (die || help) {
        std::cerr << "\n" << STAT_USAGE_MESSAGE;
        die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }
}