#include "bxfindsv.hpp"
#include "SeqLib/BamReader.h"

namespace opt {
    static std::vector<std::string> bams; // the bam to analyze
    static bool verbose = false;
}



static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools findsv bam1.bam bam2.bam ... bamN.bam\n"
                "Description: Find low-covered SVs  \n"
                "\n";


static void parseOptions(int argc, char** argv);


void runFindSV(int argc, char** argv) {
    parseOptions(argc, argv);
    SeqLib::BamReader reader;
    std::unordered_map<int, std::unordered_map<int, int>> sv_map;
    for (auto bam : opt::bams) {
        if (!reader.Open(bam)) {
            std::cerr << "Failed to open bam: " << bam << std::endl;
            exit(EXIT_FAILURE);
        }

        SeqLib::BamRecord r;
        int count = 0;
        while (reader.GetNextRecord(r)) {
            count++;
            if (count % 1000000 == 0) {
                std::cout << count << " records processed." << std::endl;
            }
            if (!r.MappedFlag() || r.NumSoftClip() > 0 || r.NumHardClip() > 0) {
                continue;
            }
            int start = r.Position();
            for (auto c_data : r.GetCigar()) {
                if (c_data.Type() == 'M') {
                    start += c_data.Length();
                    continue;
                }
                if (c_data.Type() == 'D') {
                    sv_map[start][c_data.Length()]++;
                    start += c_data.Length();
                    continue;
                }
                break;
            }
        }

    }
    std::ofstream out("out.txt", std::ofstream::out);
    for (auto pos : sv_map) {
        for (auto sv : pos.second) {
            if (sv.second > 5) {
                out << pos.first << " " << sv.first << " " << sv.second << std::endl;
            }
        }
    }
}


static void parseOptions(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        opt::bams.push_back(argv[i]);
    }
}
