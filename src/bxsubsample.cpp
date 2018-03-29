#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <getopt.h>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "bxcommon.h"
#include "bxsubsample.h"


namespace opt {

    static std::string bam; // the bam to split
    static double ratio; // the bam to split
    static std::string out_bam; // unique prefix for output
    static bool verbose = false;
}


static const char* shortopts = "hvro:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { "out-bam",                 required_argument, NULL, 'o' },
        { "ratio",                   required_argument, NULL, 'r' },
        { "verbose",                 no_argument, NULL, 'v' },
        { NULL, 0, NULL, 0 }
};


static const char *SUBSAMPLE_USAGE_MESSAGE =
        "Usage: bxtools subsample <BAM> -r <float> -o <out-BAM> \n"
                "Description: subsample bam files by excluding (1-r)*total_barcodes\n"
                "\n"
                "  General options\n"
                "  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
                "  -h, --help                           Display this help and exit\n"
                "  -o, --out-bam                        Output bam-file\n"
                "  -r, --ratio                          Output bam-file\n"
                "\n";

void parseSubsampleOptions(int argc, char** argv) {

    bool die = false;

    if (argc < 2)
        die = true;
    else
        opt::bam = std::string(argv[1]);

    bool help = false;
    std::stringstream ss;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'o': arg >> opt::out_bam; break;
            case 'r': arg >> opt::ratio; break;
            case 'v': opt::verbose = true; break;
        }
    }

    if (die || help) {
        std::cerr << "\n" << SUBSAMPLE_USAGE_MESSAGE;
        die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }
}

void fillBarcodeSet(std::unordered_set<std::string> &barcodes) {
    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }
    SeqLib::BamRecord r;
    while (reader.GetNextRecord(r)) {
        std::string bx;
        r.GetTag("BX", bx);
        if (bx.empty()) {
            continue;
        } else {
            barcodes.insert(bx);
        }
    }
    reader.Close();
}

void runSubsample(int argc, char** argv) {
    parseSubsampleOptions(argc, argv);
    std::unordered_set<std::string> barcodes;
    fillBarcodeSet(barcodes);
    int total_barcodes = barcodes.size();
    int target_barcodes = total_barcodes * opt::ratio;
    std::cout << target_barcodes << " out of " << total_barcodes << " will be kept" << std::endl;
    std::unordered_set<std::string> barcodes_to_keep;
    int i = 0;
    for (auto barcode : barcodes) {
        barcodes_to_keep.insert(barcode);
        i++;
        if (i == target_barcodes) {
            break;
        }
    }
    // opeen the BAM
    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }
    SeqLib::BamWriter writer;
    writer.Open(opt::out_bam);
    writer.SetHeader(reader.Header());
    writer.WriteHeader();
    SeqLib::BamRecord r;

    while (reader.GetNextRecord(r)) {
        std::string bx;
        r.GetTag("BX", bx);
        if (bx.empty()) {
            writer.WriteRecord(r);
        } else {
            if (barcodes_to_keep.find(bx) != barcodes_to_keep.end()) {
                writer.WriteRecord(r);
            }
        }
    }

    writer.Close();
    reader.Close();
}
