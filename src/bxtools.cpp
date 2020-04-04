/* bxtools - Tools for analyzing 10X genomics data
 * Copyright 2016 Jeremiah Wala
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the MIT license
 */

#include <iostream>
#include <bxsplit.h>
#include <bxsplit2.h>
#include <bxstats.h>
#include <bxtile.h>
#include <bxrelabel.h>
#include <bxconvert.h>
#include <bxmol.h>
#include <bxgroup.h>
#include <bxextract.h>
#include <bxfilter.h>
#include <bxsubsample.h>
#include <bxbamtofastq.h>
#include <bxfindsv.hpp>

static const char *USAGE_MESSAGE =
"Program: bxtools \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: snowman <command> [options]\n\n"
"Commands:\n"
"           split          Split a BAM into multiple BAMs, one per BX tag\n"
"           bamtofastq     Extract reads from bam file with BX barcode\n"
"           stats          Collect BX-level statistics across a BAM\n"
"           tile           Collect BX-level coverage in tiles or regions along genome\n"
"           group          Group together BX tags into molecules\n"
"           relabel        Move BX barcodes from BX tags (e.g. BX:TAATACG) to qname_TAATACG\n"
"           mol            Output BED with footprint of each molecule (from MI tag)\n"
"           convert        Flip the BX tag and chromosome, so as to allow for a BX-sorted and indexable BAM\n"
"           extract        Extract reads from BAM-file with given barcodes\n"
"           filter         Filter reads from BAM-file by quality\n"
"           split-by-ref   Create list of barcodes for each reference sequence \n"
"           subsample      Create list of barcodes for each reference sequence \n"

        "\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {

  if (argc <= 1) {
    std::cerr << USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << USAGE_MESSAGE;
      return 0;
    } else if (command == "split") {
      runSplit(argc -1, argv + 1);
    } else if (command == "stats") {
      runStat(argc -1, argv + 1);
    } else if (command == "tile") {
      runTile(argc -1, argv + 1);
    } else if (command == "relabel") {
      runRelabel(argc -1, argv + 1);
    } else if (command == "convert"){
      runConvert(argc -1, argv + 1);
    } else if (command == "group") {
      runGroup(argc -1, argv + 1);
    } else if (command == "mol") {
      runMol(argc -1, argv + 1);
    } else if (command == "extract") {
      runExtract(argc -1, argv + 1);
    } else if (command == "filter") {
      runFilter(argc - 1, argv + 1);
    } else if (command == "split-by-ref") {
      runSplitByReference(argc - 1, argv + 1);
    } else if (command == "subsample") {
      runSubsample(argc - 1, argv + 1);
    } else if (command == "bamtofastq") {
      runBamToFastq(argc - 1, argv + 1);
    } else if (command == "findsv") {
      runFindSV(argc - 1, argv + 1);
    }
    else {
      std::cerr << USAGE_MESSAGE;
      return 0;
    }
  }

  return 0;

}
