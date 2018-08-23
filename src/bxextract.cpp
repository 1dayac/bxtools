#include "bxextract.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "dirent.h"


namespace opt {

    static std::string bam; // the bam to analyze
    static bool verbose = false;
    static std::string folder_with_barcode_files; // file with list of tags to be extracted
    static std::string folder_with_small_bams;
}

static const char* shortopts = "hv:";
static const struct option longopts[] = {
        { "help",                    no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 }
};


static const char *STAT_USAGE_MESSAGE =
        "Usage: bxtools extract in.bam <folder_with_barcode_lists> <folder_with_bams>\n"
                "Description: Extract all reads from in.bam with barcodes from list.txt \n"
                "\n"
                "  General options\n"
                "  -v, --verbose                        Set verbose output\n"
                "\n";

static void parseOptions(int argc, char** argv);

void fillBarcodeMap(std::unordered_map<std::string, std::vector<std::string>> &barcodes_to_filter, std::unordered_map<std::string, SeqLib::BamWriter> &writers) {


    DIR *dirp = opendir(opt::folder_with_barcode_files.c_str());
    dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
        std::string filename(dp->d_name);
        std::cout << filename << std::endl;
        if (filename == "." || filename == ".." )
            continue;

        std::string path = opt::folder_with_barcode_files + "/" + filename;
        std::ifstream in(path);
        std::string barcode;
        std::vector<std::string> barcodes;
        int count = 0;
        while (in >> barcode) {
            barcodes.push_back(barcode);
            count++;
            if (count == 100)
                break;
        }
        if (barcodes.size() > 5) {
            writers[filename.substr(0,filename.length() - 4)] = SeqLib::BamWriter();
            for (auto barcode : barcodes) {
                barcodes_to_filter[barcode].push_back(filename.substr(0,filename.length() - 4));
            }
        }
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



    std::unordered_map<std::string, std::vector<std::string>> barcodes_to_filter;
    std::unordered_map<std::string, SeqLib::BamWriter> writers;
    std::unordered_map<std::string, std::vector<std::shared_ptr<SeqLib::BamRecord> > > records;
    fillBarcodeMap(barcodes_to_filter, writers);

    for (auto& writer : writers) {
        writer.second.Open(opt::folder_with_small_bams + writer.first + ".bam");
        writer.second.SetHeader(reader.Header());
        writer.second.WriteHeader();
        writer.second.Close();
    }


    std::vector<std::shared_ptr<SeqLib::BamRecord>> all_records;
    // loop and filter
    SeqLib::BamRecord r;
    size_t count = 0;
    while (reader.GetNextRecord(r)) {
        count++;
        if (count % 1000 == 0)
            std::cout << count << " alignments are processed" << std::endl;
        std::string bx;
        bool tag_present = r.GetZTag("BX", bx);
        if (!tag_present)
            continue;
        all_records.push_back(std::make_shared<SeqLib::BamRecord>(r));
        for (auto ids : barcodes_to_filter[bx]) {
            records[ids].push_back(all_records.back());
        }
        if (count % 100000 == 0) {
            for (auto& writer : writers) {
                writer.second.Open(opt::folder_with_small_bams + writer.first + ".bam", "ba");
                for (auto& rec : records[writer.first]) {
                    writer.second.WriteRecord(*rec);
                }
                writer.second.Close();
                records[writer.first].clear();
            }
        }
    }
    for (auto& writer : writers) {
        writer.second.Open(opt::folder_with_small_bams + writer.first + ".bam", "ba");
        for (auto& rec : records[writer.first]) {
            writer.second.WriteRecord(*rec);
        }
        writer.second.Close();
        records[writer.first].clear();
    }


}


static void parseOptions(int argc, char** argv) {

    bool die = false;
    bool help = false;

    if (argc != 4)
        die = true;
    else {
        opt::bam = std::string(argv[1]);
        opt::folder_with_barcode_files = std::string(argv[2]);
        opt::folder_with_small_bams = std::string(argv[3]);
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