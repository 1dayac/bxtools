//
// Created by dmm2017 on 7/24/18.
//
#include "bxbamtofastq.h"
#include "bxcommon.h"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <vector>
#include "SeqLib/BamReader.h"

namespace opt {

    static std::string bam; // the bam to analyze
    static bool verbose = false;
    static std::string output_folder;
}

static const char* shortopts = "hv:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 }
};

static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools bamtofastq in.bam <folder_to_output_reads (should exist)>\n"
                "Description: Convert bam-file to fastq file and keep BX tag \n"
                "\n"
                "  General options\n"
                "  -v, --verbose                        Set verbose output\n"
                "\n";

static void parseOptions(int argc, char** argv);



void processReadPair(std::vector<SeqLib::BamRecord> &records, std::ofstream &first, std::ofstream &second) {
    std::string first_read = "";
    std::string second_read = "";
    std::string first_qual = "";
    std::string second_qual = "";
    std::string bx;

    if (records.size() == 0) {
        return;
    }
    std::string read_name = records[0].Qname();
    bool tag_present = records[0].GetZTag("BX", bx);

    for (auto record : records) {
        auto cigar = record.GetCigar();
        int start_offset = 0;
        int end_offset = 0;
        if (cigar.size() != 0 && cigar.front().Type() == 'H') {
            start_offset = cigar.front().Length();
        }
        if (cigar.size() != 0 && cigar.back().Type() == 'H') {
            end_offset = cigar.back().Length();
        }

        int total_length = cigar.TotalLength();

        for (auto cigar_field : cigar) {
            if (cigar_field.Type() == 'D')
                total_length -= cigar_field.Length();
        }
        if (record.FirstFlag()) {

            if (first_read.size() < total_length) {
                first_read.resize(total_length, '?');
                first_qual.resize(total_length, '?');
            }

            first_read.replace(start_offset, total_length - start_offset - end_offset, record.Sequence());
            first_qual.replace(start_offset, total_length - start_offset - end_offset, record.Qualities());
        } else {
            if (second_read.size() < total_length) {
                second_read.resize(total_length, '?');
                second_qual.resize(total_length, '?');
            }
            second_read.replace(start_offset, total_length - start_offset - end_offset, record.Sequence());
            second_qual.replace(start_offset, total_length - start_offset - end_offset, record.Qualities());
        }
    }

    first << "@" << read_name << (tag_present ? " BX:Z:" + bx : "")  << std::endl;
    first << first_read << std::endl;
    first << "+" << std::endl;
    first << first_qual << std::endl;

    second << "@" << read_name << (tag_present ? " BX:Z:" + bx : "") << std::endl;
    second << second_read << std::endl;
    second << "+" << std::endl;
    second << second_qual << std::endl;
}



void runBamToFastq(int argc, char** argv) {
    parseOptions(argc, argv);

    std::string basename = opt::bam.substr(opt::bam.rfind("/") == std::string::npos ? 0 : opt::bam.rfind("/") + 1,
                                           opt::bam.length() - (opt::bam.rfind("/") == std::string::npos ? 0 : opt::bam.rfind("/")) - 5);
    std::string left_fastq = opt::output_folder + "/" + basename + "_R1.fastq";
    std::ofstream out(left_fastq, std::ofstream::out);
    std::string right_fastq = opt::output_folder + "/" + basename + "_R2.fastq";
    std::ofstream out2(right_fastq, std::ofstream::out);

    SeqLib::BamReader reader;
    if (!reader.Open(opt::bam)) {
        std::cerr << "Failed to open bam: " << opt::bam << std::endl;
        exit(EXIT_FAILURE);
    }
    SeqLib::BamRecord r;
    std::vector<SeqLib::BamRecord> records;
    std::string current_name = "";
    while (reader.GetNextRecord(r)) {

        if (current_name != r.Qname()) {
            processReadPair(records, out, out2);
            records.clear();
            records.push_back(r);
            current_name = r.Qname();
        } else {
            records.push_back(r);
        }
    }
    processReadPair(records, out, out2);
    out.close();
    out2.close();
}

static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc != 3)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
        opt::output_folder = std::string(argv[2]);
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