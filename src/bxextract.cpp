#include "bxextract.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

namespace opt {

    static std::string bam; // the bam to analyze
    static bool verbose = false;
    static std::string barcode_file; // file with list of tags to be extracted
}

static const char* shortopts = "hv:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 }
};


static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools extract in.bam list.txt > out.bam\n"
                "Description: Extract all reads from in.bam with barcodes from list.txt \n"
                "\n"
                "  General options\n"
                "  -v, --verbose                        Set verbose output\n"
                "\n";

static void parseOptions(int argc, char** argv);

void fillBarcodeMap(std::unordered_set<std::string> &barcodes_to_filter) {
    std::ifstream in(opt::barcode_file);
    std::string barcode;
    while (in >> barcode) {
        barcodes_to_filter.insert(barcode);
    }
}

void runExtract(int argc, char** argv) {
    parseOptions(argc, argv);

    // open the BAM
    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }

    std::unordered_set<std::string> barcodes_to_filter;
    fillBarcodeMap(barcodes_to_filter);

    SeqLib::BamWriter writer;
    writer.Open("-");
    writer.SetHeader(reader.Header());
    writer.WriteHeader();
    // loop and filter
    SeqLib::BamRecord r;
    size_t count = 0;
    while (reader.GetNextRecord(r)) {
        std::string bx;
        bool tag_present = r.GetZTag("BX", bx);
        if (!tag_present || barcodes_to_filter.count(bx) == 0)
            continue;
        writer.WriteRecord(r);
    }
}

static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc != 3)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
        opt::barcode_file = std::string(argv[2]);
    }

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'v': opt::verbose = true; break;
            case 'h': help = true; break;
        }
    }

    if (die || help) {
        std::cerr << "\n" << STAT_USAGE_MESSAGE;
        die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }
}