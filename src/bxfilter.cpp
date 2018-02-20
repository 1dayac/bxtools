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
    static double max_soft_clipping = 1.0;
    static double max_hard_clipping = 1.0;
    static bool filter_bad = false;
}

static const char* shortopts = "c:q:s:hvb";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 },
        {"mapping_quality", optional_argument, NULL, 'q'},
        {"max_soft_clipping", optional_argument, NULL, 's'},
        {"max_hard_clipping", optional_argument, NULL, 'c'},
        {"filter_bad", optional_argument, NULL, 'b'}
};


static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools filter in.bam -q X > out.bam\n"
                "Description: Extract all reads from in.bam that satisfy parameters \n"
                "\n"
                "  General options\n"
                "-q, --mapping_quality                    Filter read pairs with mapping quality of any read lower than [q]\n"
                "-s, --max_soft_clipping                  Filter read pairs with any read with portion of soft clipped pairs more than [s]\n"
                "-c, --max_hard_clipping                  Filter read pairs with any read with portion of hard clipped pairs more than [c]\n"
                "-b, --filter_bad                         Filter read pairs that don't satisfy given conditions\n"
                "  -v, --verbose                          Set verbose output\n"
                "\n";

static void parseOptions(int argc, char** argv);
static bool CheckConditions(const std::vector<SeqLib::BamRecord> &records);
static bool AdditionalChecks(const SeqLib::BamRecord &record);



static bool AdditionalChecks(const SeqLib::BamRecord &record) {
    if (record.NumMatchBases() < 50) {
        return false;
    }
    if (record.AlignmentPosition() > 25 && record.Length() - record.AlignmentEndPosition() > 25) {
        return false;
    }
    return true;
}

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
    size_t count = 0;
    std::cerr << "max-soft-clipping " << opt::max_soft_clipping << std::endl;
    std::cerr << "filter bad " << opt::filter_bad << std::endl;

    std::string read_id = "";
    std::vector<SeqLib::BamRecord> bam_records;
    while (reader.GetNextRecord(r1)) {

        if (read_id == r1.Qname()) {
            bam_records.push_back(r1);
        } else {
            if (CheckConditions(bam_records)) {
                for (const auto& record : bam_records) {
                    writer.WriteRecord(record);
                }
            }
            bam_records.clear();
            bam_records.push_back(r1);
            read_id = r1.Qname();
        }
    }
    if (CheckConditions(bam_records)) {
        for (const auto& record : bam_records) {
            writer.WriteRecord(record);
        }
    }
}

static bool CheckConditions(const std::vector<SeqLib::BamRecord> &records) {
    if (!opt::filter_bad) {
        for (const auto &record : records) {
            if (record.MapQuality() < opt::mapping_quality) {
                return false;
            }
        }
        for (const auto &record : records) {
            if (record.NumSoftClip()/(double)record.Length() > opt::max_soft_clipping) {
                return false;
            }
        }

        for (const auto &record : records) {
            if (record.NumHardClip()/(double)record.Length() > opt::max_hard_clipping) {
                return false;
            }
        }
        return true;
    } else {
        for (const auto &record : records) {
            if (record.NumMatchBases()/(double)record.Length() < 0.1) {
                std::cerr << "Filtered: bad mapping quality" << std::endl;
                return true;
            }
        }
        for (const auto &record : records) {
            if (record.NumSoftClip()/(double)record.Length() > opt::max_soft_clipping) {
                std::vector<std::string> chrom;
                for (const auto &record2 : records) {
                    chrom.push_back(record2.ChrName());
                }
                bool allAreEqual =
                        find_if(chrom.begin() + 1,
                                chrom.end(),
                                bind1st(std::not_equal_to<std::string>(), chrom.front())) == chrom.end();
                if (!allAreEqual) {
                    return false;
                }

                for (const auto &record2 : records) {
                    if (AdditionalChecks(record2)) {
                        return false;
                    }
                }

                std::cerr << "Filtered: soft clips" << std::endl;
                return true;
            }
        }

//        for (const auto &record : records) {
//            if (record.NumHardClip()/(double)record.Length() > opt::max_hard_clipping) {
//               std::cerr << "Filtered: hard clips" << std::endl;
//                return true;
//            }
//        }
        return false;
    }

}

static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc < 2)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
    }

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'q': arg >> opt::mapping_quality; break;
            case 'v': opt::verbose = true; break;
            case 's': arg >> opt::max_soft_clipping; break;
            case 'c': arg >> opt::max_hard_clipping; break;
            case 'b' : opt::filter_bad = true; break;
        }
    }

    if (die || help) {
        std::cerr << "\n" << STAT_USAGE_MESSAGE;
        die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }
}