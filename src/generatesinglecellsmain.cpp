#include "utils.h"
#include "arg_parser.h"
#include "singlecellgeneration.h"
#include <string>

int main(int argc, char **argv) {
    int numSCS = 1000;
    double read_depth = .05;
    double alpha_fp = .001;
    std::string out_dir;
    std::string in_dir;
    int k = 10; //number of segments
    int m = 2; //number of samples
    double cna_error = 0;
    double plsOne = .3;
    double minusOne = .3;
    double copyNeutral = .4;

    lemon2::ArgParser ap(argc, argv);
    ap.refOption("num_cells", "The number of cells to simulate with the single cell generation (default: 1000)", numSCS,
                 false)
            .refOption("read_depth", "The read_depth for the single cell generation (default: .05)", read_depth, false)
            .refOption("alpha_fp", "The sequencing error for single cell generation (default .001)", alpha_fp, false)
            .refOption("out_dir", "The output directory for single cell generation (default: results)", out_dir, false)
            .refOption("in_dir", "The input directory of files for single cell generation", in_dir)
            .refOption("k", "Number of segments (default: 10)", k, false)
            .refOption("m", "Number of samples (default: 2)", m, false)
            .refOption("e", "The error rate for CNA data (default: 0)", cna_error, false)
            .refOption("add", "The proportion of copy number errors adding an allele (default .3)", plsOne, false)
            .refOption("sub", "The proportion of copy number errors subtracting an allele (default .3)", minusOne, false)
            .refOption("neutral", "The proportion of copy number errors that are total copy number neutral (default .4)", copyNeutral, false);

    ap.parse();


    try {
        std::cerr << "Generating single cell data..." << std::endl;
        for (int i = 0; i < m; i++) {
            std::cerr << "starting sample " << i << std::endl;
            SingleCell sc(numSCS, read_depth, alpha_fp, out_dir, k, m, cna_error);
            sc.main(std::cout, in_dir, i, cna_error, plsOne, minusOne, copyNeutral);
        }
    }
    catch (std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
