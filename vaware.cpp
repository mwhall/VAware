#include <string>
#include <mutex>
#include <thread>
#include <iostream> 
#include <fstream>
#include "vaware_config.h"
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>

using namespace seqan;

// ==========================================================================
// Typedefs, structs
// ==========================================================================

typedef StringSet<IupacString> IupacStringSet;                    // container for strings, restricted to IupacStrings
typedef Align<IupacString> IupacAlign;                            // alignment object
typedef Row<IupacAlign>::Type IupacRow;
typedef StringSet<Dna5String> Dna5StringSet;
typedef Score<int, ScoreMatrix<Iupac, Default> > IupacScoringScheme;  //Scoring scheme

namespace seqan {
// Define the degenerate scoring matrix

// We have to create a new specialization of the ScoringMatrix_ class
// for the DNA alphabet.  For this, we first create a new tag.
struct UserDefinedMatrix {};

// Then, we specialize the class ScoringMatrix_ for the Iupac alphabet.
// NOTE: N is treated as a miss when it's on the template, but not the primer
template <>
struct ScoringMatrixData_<int, Iupac, UserDefinedMatrix>
{
    enum
    {
        VALUE_SIZE = ValueSize<Iupac>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData()
    {
        // The user defined data table.
        static int const _data[TAB_SIZE] =
        { //U  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
            1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, //U
            0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, //A
            0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, //C
            0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, //M
            0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, //G
            0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, //R
            0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, //S
            0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, //V
            1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, //T
            1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, //W
            1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, //Y
            1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, //H
            1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, //K
            1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, //D
            1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, //B
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1  //N
        };
        return _data;
    }

};
} //End define degenerate scoring matrix

// ==========================================================================
// Global variables
// ==========================================================================

std::mutex m;

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function count3PrimeMistmatches()
// --------------------------------------------------------------------------

unsigned int count3PrimeMismatches(IupacAlign& trimmedAlign, IupacScoringScheme& scoreScheme, int k) {

    unsigned int numMismatches = 0;

    Gaps<IupacString> templateRow = row(trimmedAlign, 0);
    Gaps<IupacString> primerRow = row(trimmedAlign, 1);

    for (int i = length(templateRow) - 1; i >= length(templateRow) - k; --i)
    {
        if (score(scoreScheme, templateRow[i], primerRow[i]) == 0) {
            ++numMismatches;
        }
    }

    return numMismatches;

}

// --------------------------------------------------------------------------
// Function printResultsRow()
// --------------------------------------------------------------------------

void printPreamble(CharString& seqFileName, CharString& alignFileName,
                   CharString& insertFileName, CharString& taxonFilter,
                   IupacString& fwdPrimer, IupacString& revPrimer, 
                   bool excludePrimers, int k, int nThreads)
{
// Print to stdout the selected parameters
    std::cout << "########################################################" << std::endl;
    std::cout << "# vaware Version: " << VAware_VERSION << std::endl;
    std::cout << "# Reference File (input): " << seqFileName << std::endl;
    std::cout << "# Alignment File (output): ";
    (alignFileName != "NULL") ? (std::cout << alignFileName) : (std::cout << "not output");
    std::cout << std::endl;
    std::cout << "# Insert File (output): ";
    (insertFileName != "NULL") ? (std::cout << insertFileName) : (std::cout << "not output");
    std::cout << std::endl;
    std::cout << "# Taxon Filter: ";
    (taxonFilter != "NULL") ? (std::cout << taxonFilter) : (std::cout << "no filter");
    std::cout << std::endl;
    std::cout << "# Forward Primer (5'-3'): " << fwdPrimer << std::endl;
    std::cout << "# Reverse Primer (5'-3'): " << revPrimer << std::endl;
    std::cout << "# Exclude Primers from Insert Output: ";
    (excludePrimers) ? (std::cout << "yes") : (std::cout << "no");
    std::cout << std::endl;
    std::cout << "# k Mismatch Distance from 3' End: " << k << std::endl;
    std::cout << "# Number of Threads: " << nThreads << std::endl;
    std::cout << "########################################################" << std::endl;
}

void printResultsRow(CharString& label, AlignmentStats& fwdAlignStats, AlignmentStats& revAlignStats,
                     IupacAlign& fwdAlign, IupacAlign& revAlign, IupacScoringScheme& scoreScheme,
                     int k, bool excludePrimers, std::ofstream& insertOut, std::ofstream& alignOut)
{
//     - Number of gaps within k bp of the 3' end


    //SILVA label parsing:
    // Split up the label into ID and taxonomy string
    // Note: SILVA has the host (?) as the last entry in the tax string, so we cut that off
    std::string strLabel = toCString(label);
    std::size_t firstSpace = strLabel.find_first_of(" ");
    std::string silvaId = strLabel.substr(0, firstSpace);
    std::string taxonomy = strLabel.substr(firstSpace + 1, strLabel.find_last_of(";") - firstSpace - 1);

    unsigned int templateLength = unclippedLength(row(fwdAlign, 0));
    //Forward primer stats:
    unsigned int numFwdPosScores = fwdAlignStats.numPositiveScores;
    unsigned int numFwdNegScores = fwdAlignStats.numNegativeScores;
    unsigned int numFwdGaps = fwdAlignStats.numGaps;
    unsigned int fwdLeadGaps = countLeadingGaps(row(fwdAlign, 1));
    unsigned int fwdBeginPos = clippedBeginPosition(row(fwdAlign, 0));
    unsigned int fwdEndPos = clippedEndPosition(row(fwdAlign, 0));
    unsigned int fwdEndMismatches = count3PrimeMismatches(fwdAlign, scoreScheme, k);

    //Reverse primer stats:
    unsigned int numRevPosScores = revAlignStats.numPositiveScores;
    unsigned int numRevNegScores = revAlignStats.numNegativeScores;
    unsigned int numRevGaps = revAlignStats.numGaps;
    unsigned int revLeadGaps = countLeadingGaps(row(revAlign, 1));
    unsigned int revBeginPos = clippedBeginPosition(row(revAlign, 0));
    unsigned int revEndPos = clippedEndPosition(row(revAlign, 0));
    unsigned int revEndMismatches = count3PrimeMismatches(revAlign, scoreScheme, k);

    IupacString fwdTemplate = source(row(fwdAlign, 0));
    //Check to ensure that the substring we are trying to find can exist
    // under both primer-included and primer-excluded trimming schemes
    bool printSeq = false;
    if ((fwdBeginPos < templateLength - revBeginPos) &&
        (fwdEndPos < templateLength - revEndPos)) {
        printSeq = true;
    }

    // TODO: Proper reporting of either a fwd or rev primer fail
    if ((numFwdPosScores == 0) || (numRevPosScores == 0)) {
        printSeq = false;   
    }

    if (printSeq) {
        Infix<IupacString >::Type primerInsert;
        if (excludePrimers) {
            primerInsert = infix(fwdTemplate, fwdEndPos, templateLength - revEndPos);
        } else {
            primerInsert = infix(fwdTemplate, fwdBeginPos, templateLength - revBeginPos);
        }
        std::cout << silvaId << "\t" << taxonomy << "\t";
        std::cout << numFwdPosScores << "\t" << numFwdNegScores << "\t";
        std::cout << numFwdGaps << "\t" << fwdEndMismatches << "\t";
        std::cout << fwdBeginPos << "\t" << fwdEndPos << "\t";
        std::cout << numRevPosScores << "\t" << numRevNegScores << "\t";
        std::cout << numRevGaps << "\t" << revEndMismatches << "\t";
        std::cout << templateLength - revBeginPos << "\t" << templateLength - revEndPos << "\t";

        if (excludePrimers) {
            std::cout << templateLength - fwdEndPos - revEndPos << std::endl;
        } else {
            std::cout << templateLength - fwdBeginPos - revBeginPos << std::endl;
        }

        insertOut << ">" << silvaId << " " << taxonomy << std::endl;
        insertOut << primerInsert << std::endl;

        alignOut << silvaId << " " << taxonomy << std::endl;
        alignOut << "Forward Primer: " << std::endl << fwdAlign;
        alignOut << "Reverse Primer: " << std::endl << revAlign;

    // Alignment may be in reverse position, or bad/no primer match
    // No insert possible, so report NAs
    } else {
         std::cout << silvaId << "\t" << taxonomy << "\t";
        std::cout << "NA\tNA\tNA\tNA\tNA\tNA\t";
        std::cout << "NA\tNA\tNA\tNA\tNA\tNA\tNA" << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function alignPrimer()
// --------------------------------------------------------------------------

void alignPrimer(bool forwardPrimer, IupacString& templateSequence, IupacString& primer, 
                 IupacScoringScheme& scoreScheme, IupacAlign& align, IupacStringSet& sequences)
{
    // If forward primer, all good. If reverse primer, need template to be the - strand
    // Maintain primers in 5'-3' position
    if (forwardPrimer) {
        appendValue(sequences, templateSequence);
    } else {
        // Make a copy so we don't permanently mutate templateSequence
        IupacString revCompTemplate(templateSequence);
        reverseComplement(revCompTemplate);
        appendValue(sequences, revCompTemplate);
    }

    appendValue(sequences, primer);
    setStrings(align, sequences);

    int score = globalAlignment(align, scoreScheme, 
                                AlignConfig<true, true, true, true>(), 
                                LinearGaps());

}

// --------------------------------------------------------------------------
// Function trimAlignment()
// --------------------------------------------------------------------------

void trimAlignment(IupacAlign& align)
{
    //Assumes that the 0th sequence is the template, 1st sequence is the primer
    uint64_t leadingGaps {countLeadingGaps(row(align, 1))};
    uint64_t trailingGaps {countTrailingGaps(row(align, 1))};
    uint64_t alignLength {length(row(align, 1))};

    setClippedBeginPosition(row(align, 0), leadingGaps);
    setClippedBeginPosition(row(align, 1), leadingGaps);
    setClippedEndPosition(row(align, 0), alignLength - trailingGaps);
    setClippedEndPosition(row(align, 1), alignLength - trailingGaps);

}

void processSequence(CharString id, IupacString templateSeq,
                     IupacString fwdPrimer, IupacString revPrimer,
                     IupacScoringScheme& scoringScheme, int k, bool excludePrimers,
                     std::ofstream& alignOut, std::ofstream& insertOut)
{
        //Align forward

        IupacStringSet fwdSequences;
        IupacAlign fwdAlign;
        alignPrimer(true, templateSeq, fwdPrimer, scoringScheme, fwdAlign, fwdSequences);

        trimAlignment(fwdAlign);

        AlignmentStats fwdAlignStats;
        computeAlignmentStats(fwdAlignStats, fwdAlign, scoringScheme);

        //Align reverse

        IupacStringSet revSequences;
        IupacAlign revAlign;
        alignPrimer(false, templateSeq, revPrimer, scoringScheme, revAlign, revSequences);

        trimAlignment(revAlign);

        AlignmentStats revAlignStats;
        computeAlignmentStats(revAlignStats, revAlign, scoringScheme);
        //Lock before printing, to avoid lines printing on top of one another
        m.lock();
        printResultsRow(id, fwdAlignStats, revAlignStats, fwdAlign, revAlign, 
                        scoringScheme, k, excludePrimers, insertOut, alignOut);
        m.unlock();

}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{

    // ----------------------------------------------------------------------
    // Argument parsing
    // ----------------------------------------------------------------------

    ArgumentParser parser("vaware");

    // Version and Descriptions
    setVersion(parser, VAware_VERSION);
    setShortDescription(parser, "PCR primer alignment tool");
    setDate(parser, "January 2018");
    addDescription(parser, "Tool to analyse the alignment of specified PCR "
                           "primer pairs against a specified sequence file");
    addUsageLine(parser, "\\fI-i INPUT_FASTA\\fP \\fI-f FWD_PRIMER_SEQ\\fP \\fI-r REV_PRIMER_SEQ\\fP [\\fIOPTIONS\\fP]");

    // Threading
    addOption(parser, ArgParseOption(
        "t", "nthreads", "Number of threads to use.",
        ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "nthreads", 1);

    // IO
    addOption(parser, ArgParseOption(
        "i", "input-file", "Path to the input reference FASTA file.",
        ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "i");

    addOption(parser, ArgParseOption(
        "a", "align-file", "Path to the optional alignment output file.",
        ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "align-file", "NULL");
    addOption(parser, ArgParseOption(
        "o", "insert-file", "Path to the optional insert output FASTA file.",
        ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "insert-file", "NULL");

    // Primers
    addOption(parser, ArgParseOption(
        "f", "forward-primer", "Forward primer sequence in 5'-3' orientation (can include degeneracy).",
        ArgParseArgument::STRING, "STRING"));
    setRequired(parser, "f");
    addOption(parser, ArgParseOption(
        "r", "reverse-primer", "Reverse primer sequence in 5'-3' orientation (can include degeneracy).",
        ArgParseArgument::STRING, "STRING"));
    setRequired(parser, "r");
    addOption(parser, ArgParseOption(
        "x", "exclude-primers", "Exclude primer sequences from insert output."));

    // Other options
    addOption(parser, ArgParseOption(
        "k", "mismatch-distance", "Number of nucleotides from the 3' end to consider ``3' mismatches''.",
        ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "mismatch-distance", 4);
    addOption(parser, ArgParseOption(
        "m", "match", "Optional text that must appear in taxonomic classification to be included in results. Default: Include all taxa.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "match", "NULL");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    int nThreads = 1;
    getOptionValue(nThreads, parser, "nthreads");

    IupacString fwdPrimer;
    getOptionValue(fwdPrimer, parser, "forward-primer");

    IupacString revPrimer;
    getOptionValue(revPrimer, parser, "reverse-primer");

    bool excludePrimers = false;
    getOptionValue(excludePrimers, parser, "exclude-primers");

    int k = 4;
    getOptionValue(k, parser, "mismatch-distance");

    CharString taxonFilter;
    getOptionValue(taxonFilter, parser, "match");

    CharString alignFileName;
    std::ofstream alignOut(nullptr);
    getOptionValue(alignFileName, parser, "align-file");

    if (alignFileName != "NULL") {
        alignOut.open(toCString(alignFileName));
    }

    CharString insertFileName;
    std::ofstream insertOut(nullptr);
    getOptionValue(insertFileName, parser, "insert-file");

    if (insertFileName != "NULL") {
        insertOut.open(toCString(insertFileName));
    }

    CharString seqFileName;
    getOptionValue(seqFileName, parser, "input-file");

    // Before doing anything, export the settings to stdout
    printPreamble(seqFileName, alignFileName, insertFileName, taxonFilter,
                  fwdPrimer, revPrimer, excludePrimers, k, nThreads);

    SeqFileIn seqFileIn(toCString(seqFileName));

    if (!open(seqFileIn, toCString(seqFileName)))
    {
        std::cerr << "ERROR: Could not open the reference file." << std::endl;
        return 1;
    }

    //Set the global alignment parameters

    //Increased gap open penalty to avoid gaps in primers
    int const gapOpenScore = -5;
    int const gapExtendScore = -1;

    IupacScoringScheme scoringScheme(gapExtendScore, gapOpenScore);
    //Set the alignment score based on Iupac nucleotide degeneracy matrix
    setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());

    std::cout << "SILVA ID\tTaxonomy\t";
    std::cout << "FP Matches\tFP Mismatches\tFP Gaps\tFP 3' Mismatches\tFP Align Begin\tFP Align End\t";
    std::cout << "RP Matches\tRP Mismatches\tRP Gaps\tRP 3' Mismatches\tRP Align Begin\tRP Align End\t";
    std::cout << "Insert Length" << std::endl;

    std::thread myThreads[nThreads];

    CharString id;
    std::string strId;
    IupacString templateSeq;

    int i = 0;

    while (!atEnd(seqFileIn))
    {
        //If there is a taxon filter, continually read until the id contains the search string
        if (taxonFilter != "NULL") {
            id = NULL;
            strId = toCString(id);
            while (strId.find(toCString(taxonFilter)) == std::string::npos && !atEnd(seqFileIn)) {
                readRecord(id, templateSeq, seqFileIn);
                strId = toCString(id);
            }
            if (atEnd(seqFileIn)) {
                break;
            }
        //If there is no filter, just read the record
        } else {
            readRecord(id, templateSeq, seqFileIn);
        }

        myThreads[i] = std::thread(processSequence, id, templateSeq, fwdPrimer, revPrimer, std::ref(scoringScheme), k, excludePrimers, std::ref(alignOut), std::ref(insertOut));

        if (i == nThreads - 1) {
            for (i = 0; i < nThreads; ++i) {
                myThreads[i].join();
            }
            i = 0;
        } else {
            ++i;
        }
    }
    //Clean up unexecuted threads
    for (int j = 0; j < i; ++j) {
        myThreads[j].join();
    }

    if (alignOut) {
        alignOut.close();
    }

    return 0;
}
