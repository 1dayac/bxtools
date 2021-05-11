//
// Created by dmm2017 on 7/24/18.
//
#include "bxamfilter.h"
#include "bxcommon.h"
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <vector>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

namespace opt {

    static std::string bam; // the bam to analyze
    static bool verbose = false;
    static std::string output_bam;
}


static const char* shortopts = "hv:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 }
};

static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools amfilter in.bam out.bam\n"
        "Description: Filter read mappings with AM=0 and reads with wrong AM tag \n"
        "\n"
        "  General options\n"
        "  -v, --verbose                        Set verbose output\n"
        "\n";

static void parseOptions(int argc, char** argv);

void runAmFilter(int argc, char** argv) {
    parseOptions(argc, argv);

    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }
    SeqLib::BamWriter writer;
    writer.Open(opt::output_bam);
    writer.SetHeader(reader.Header());
    writer.WriteHeader();
    // loop and filter
    SeqLib::BamRecord r1;

    std::unordered_set<std::string> records_to_discard;
    std::string tag;
    std::string barcode;

    int distance_diff = 15000;
    std::unordered_map<std::string, std::stack<std::pair<std::string, int>>> current_reads;
    while (reader.GetNextRecord(r1)) {
        std::string read_id = r1.Qname();
        int pos = r1.Position();
        r1.GetTag("AM", tag);
        r1.GetTag("BX", barcode);
        if (tag == "0") {
            records_to_discard.insert(read_id);
        } else {
            if (current_reads[barcode].size() == 0) {
                current_reads[barcode].push(std::make_pair(read_id, pos));
            } else if (current_reads[barcode].size() == 1) {
                if (abs(pos - current_reads[barcode].top().second) < distance_diff) {
                    current_reads[barcode].push(std::make_pair(read_id, pos));
                } else {
                    records_to_discard.insert(current_reads[barcode].top().first);
                    current_reads[barcode] = std::stack<std::pair<std::string, int>>();
                    current_reads[barcode].push(std::make_pair(read_id, pos));
                }
            } else {
                if (abs(pos - current_reads[barcode].top().second) < distance_diff) {
                    current_reads[barcode].push(std::make_pair(read_id, pos));
                } else {
                    current_reads[barcode] = std::stack<std::pair<std::string, int>>();
                    current_reads[barcode].push(std::make_pair(read_id, pos));
                }
            }
        }
    }

    for (auto it : current_reads) {
        if (it.second.size() == 1) {
            records_to_discard.insert(it.first);
        }
    }

    reader.Close();
    SeqLib::BamReader reader2;

    if (!reader2.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }

    while (reader2.GetNextRecord(r1)) {
        if (!records_to_discard.count(r1.Qname())) {
            writer.WriteRecord(r1);
        }
    }
    writer.Close();
}

static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc != 3)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
        opt::output_bam = std::string(argv[2]);
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