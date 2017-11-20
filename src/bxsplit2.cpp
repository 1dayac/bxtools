#include "bxsplit2.h"

#include "bxcommon.h"
#include <string>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>
#include <sys/stat.h>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"


namespace opt {

    static std::string bam; // the bam to split
    static std::string out_folder; // unique prefix for output
    static bool verbose = false;
}

static const char* shortopts = "hvo:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { "out-folder",              required_argument, NULL, 'o' },
        { "verbose",                 no_argument, NULL, 'v' },
        { NULL, 0, NULL, 0 }
};

static const char *SPLIT2_USAGE_MESSAGE =
        "Usage: bxtools split-by-ref <BAM/SAM/CRAM> -o out_folder \n"
                "Description: Split / count a BAM into multiple BAMs, one BAM per unique BX tag\n"
                "\n"
                "  General options\n"
                "  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
                "  -h, --help                           Display this help and exit\n"
                "  -o, --out-folder                     Folder to store output\n"
                "\n";

void parseSplit2Options(int argc, char** argv) {

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
            case 'o': arg >> opt::out_folder; break;
            case 'v': opt::verbose = true; break;
        }
    }

    if (die || help) {
        std::cerr << "\n" << SPLIT2_USAGE_MESSAGE;
        die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }
}

void runSplitByReference(int argc, char** argv) {
    parseSplit2Options(argc, argv);

    // opeen the BAM
    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }

    int err_code = mkdir(opt::out_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (err_code)
        exit(EXIT_FAILURE);

    // make a collection of writers
    std::unordered_map<std::string, std::set<std::string>> tags;
    // loop and write
    SeqLib::BamRecord r;
    size_t count = 0;
    bool hit = false;
    while (reader.GetNextRecord(r)) {
        std::string bx;
        r.GetTag("BX", bx);
        if (bx.empty()) {
            continue;
        } else {
            hit = true;
        }
        tags[r.ChrName()].insert(bx);
    }

    for (auto chr : tags) {
        std::string out_file = opt::out_folder + "/" + chr.first + ".txt";
        std::ofstream out(out_file.c_str(), std::ofstream::out);
        for (auto tag : chr.second) {
            out << tag << "\n";
        }
    }

}