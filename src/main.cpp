/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * main file
 *
 * For theory please have a look at and also cite, if you have used rasbhari in your publications:
 *
 * - Hahn L, Leimeister C-A, Ounit R, Lonardi S, Morgenstern B (2016)
 * rasbhari: Optimizing Spaced Seeds for Database Searching, Read Mapping and Alignment-Free Sequence Comparison.
 * PLoS Comput Biol 12(10):e1005107. doi:10.1371/journal.pcbi.1005107
 *
 * - B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
 * Estimating evolutionary distances between genomic sequences from spaced-word matches
 * Algorithms for Molecular Biology 10, 5. (http://www.almob.org/content/10/1/5/abstract)
 *
 *
 * @author: Lars Hahn - 23.08.2017, Georg-August-Universitaet Goettingen
 * @version: 1.4.0 08/2017
 */
#include <iostream>

#include "rasbhari.hpp"
#include "rasbopt.hpp"
#include "rasbimp.hpp"

void SecurityMessage(std::string errmsg, char* argv[], int pos);

int main(int argc, char *argv[]){
    rasb_opt::Quiet = true;
    rasb_opt::Silent = false;
    rasb_opt::Sens = true;
    rasb_opt::OutFile = "";
    double TmpDb;
    bool LengSet = false;

    bool Exit = false;
    std::string tmp;
    if(argc <= 1){
        SecurityMessage("param", argv, 0);
        Exit = true;
    }
    for(int i = 1; i < argc; i++){
        tmp = argv[i];
        if(tmp == "--version"){
            SecurityMessage("version", argv, 0);
            Exit = true;
        }
        else if(tmp == "--help"){
            SecurityMessage("param", argv, 0);
            Exit = true;
        }
    }
    if(Exit){
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        std::string parse = argv[i];
        switch (parse[1]) {
        case 'D':
        case 'd':
            if (i < argc - 1) {
                rasb_opt::parse_length(argv[i+1]);
                i++;
                LengSet = true;
            }
            break;
        case 'H':
        case 'h':
            if (i < argc - 1) {
                rasb_opt::H = atoi(argv[i + 1]);
                i++;
            }
            break;
        case 'M':
        case 'm':
            if (i < argc - 1) {
                rasb_opt::Size = atoi(argv[i + 1]);
                i++;
            }
            break;
        case 'Q':
        case 'q':
            if (i < argc - 1) {
                TmpDb = atof(argv[i + 1]);
                if(0 < TmpDb && TmpDb < 1){
                    rasb_opt::Q = TmpDb;
                }
                i++;
            }
            break;
        case 'P':
        case 'p':
            if (i < argc - 1) {
                TmpDb = atof(argv[i + 1]);
                if(0 < TmpDb && TmpDb < 1){
                    rasb_opt::P = TmpDb;
                }
                i++;
            }
            break;
        case 'S':
        case 's':
            if (i < argc - 1) {
                rasb_opt::SeqLength = atoi(argv[i + 1]);
                i++;
            }
            break;
        case 'W':
        case 'w':
            if (i < argc - 1) {
                rasb_opt::Weight = atoi(argv[i + 1]);
                i++;
            }
            break;
        case '-':
            parse = argv[i];
            if(parse == "--forcesens"){
                rasb_opt::Forcesens = true;
            }
            else if(parse == "--nosens"){
                rasb_opt::Sens = false;
            }
            else if(parse == "--notquiet") {
                rasb_opt::Quiet = false;
            }
            else if(parse == "--opt-oc"){
                if (i < argc - 1) {
                    rasb_opt::OptOc = atoi(argv[i + 1]);
                    i++;
                    if(rasb_opt::OptOc < 1){
                        rasb_opt::OptOc = 1;
                    }
                }
            }
            else if(parse == "--opt-sens"){
                if (i < argc - 1) {
                    rasb_opt::OptSens = atoi(argv[i + 1]);
                    i++;
                    if(rasb_opt::OptSens < 1){
                        rasb_opt::OptSens = 1;
                    }
                }
            }
            else if(parse == "--outfile") {
                if (i < argc - 1) {
                    rasb_opt::OutFile = argv[i + 1];
                    i++;
                }
            }
            else if(parse == "--pattern") {
                if (i < argc - 1) {
                    rasb_opt::InFile = argv[i + 1];
                    i++;
                }
            }
            else if(parse == "--permut"){
                if (i < argc - 1) {
                    rasb_opt::Limit = atoi(argv[i + 1]);
                    i++;
                }
            }
            else if(parse == "--silent"){
                rasb_opt::Silent = true;
            }
            else if(parse == "--variance") {
                rasb_opt::Oc = false;
            }
            else {
                SecurityMessage("parsing", argv, i);
            }
            break;
        default:
            SecurityMessage("parsing", argv, i);
            break;
        }
    }

    if(LengSet == false){
        if(rasb_opt::Oc){
            rasb_opt::MinDontcare = rasb_opt::Weight+3;
            rasb_opt::MaxDontcare = 2*(rasb_opt::Weight+3);
        }
        else{
            rasb_opt::MinDontcare = (1/rasb_opt::P)*rasb_opt::Weight;
        }
    }

    std::cout << "\n====================================" << std::endl;
    std::cout << "Parameter for variance calculation:" << std::endl;
    std::cout << "====================================\n" << std::endl;
    std::cout << "S                           = " << rasb_opt::SeqLength << std::endl;
    std::cout << "p                           = " << rasb_opt::P << std::endl;
    std::cout << "q                           = " << rasb_opt::Q << std::endl;
    std::cout << "H (for sensitivity)         = " << rasb_opt::H << std::endl;
    std::cout << "Pattern don't care          = " << rasb_opt::MinDontcare << " - " << rasb_opt::MaxDontcare << std::endl;
    std::cout << "Pattern weight              = " << rasb_opt::Weight << std::endl;
    std::cout << "Patternset                  = " << rasb_opt::Size << " Pattern" << std::endl;
    std::cout << "Calculation type            :";
    if(rasb_opt::Oc){
        std::cout << " Overlap complexity\n\n" << std::endl;
    }
    else{
        std::cout << " Variance\n\n" << std::endl;
    }

    std::cout << "Hillclimbing ?";
    if(rasb_opt::Limit > 0){
        std::cout << " [Y]          : " << rasb_opt::Limit << " permutations" << std::endl;
    }
    else{
        std::cout << " [N]" << std::endl;
    }
    std::cout << "Iterative    ?";
    if(rasb_opt::OptOc > 1){
        std::cout << " [Y]          : " << rasb_opt::OptOc << " iterations" << std::endl;
    }
    else{
        std::cout << " [N]" << std::endl;
    }
    std::cout << "Sensitivity  ?";
    if(rasb_opt::Sens){
        std::cout << " [Y]          : " << rasb_opt::OptSens << " loops" << std::endl;
    }
    else{
        std::cout << " [N]" << std::endl;
    }
    std::cout << "\n\n" << std::endl;

    rasbhari RasbSet;

    if(rasb_opt::Sens == false){
        if(rasb_opt::OptOc == 1){
            if(rasb_opt::Oc){
                RasbSet = rasb_implement::hillclimb_oc();
            }
            else{
                RasbSet = rasb_implement::hillclimb_var();
            }
        }
        else{
            if(rasb_opt::Oc){
                RasbSet = rasb_implement::hillclimb_oc_iterative();
            }
            else{
                RasbSet = rasb_implement::hillclimb_var_iterative();
            }
        }
    }
    else{
        if(rasb_opt::Oc){
            RasbSet = rasb_implement::hillclimb_sens_oc();
        }
        else{
            RasbSet = rasb_implement::hillclimb_sens_var();
        }
    }
    RasbSet.to_file();
}

/**
 * This function holds neccessary information like the help and version info
 * or the error message for wrong parameters.
 *
 * @param errmsg        The Error Message type.
 *
 * @param argv          Pointer to an array of char arrays, for the arg-list.
 *
 * @param pos           Integer indicating which parameter is wrong.
 */
void SecurityMessage(std::string errmsg, char* argv[], int pos){
    if (errmsg == "param") {
        std::cerr << "This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight."<< std::endl;
        std::cerr << "It is possible to improve your patternset and read patterns from a file.\n" << std::endl;
        std::cerr << "Usage:\n\t" << argv[pos] << " \t\t\t Print this help.\n" << std::endl;
        std::cerr << "\t" << argv[pos] << " {options}\t\t (E.g.: " << argv[pos] << " -m 10 -w 8 -d 6-15 -H 64 --permut 25000)\n" << std::endl;
        std::cerr << "Options:\n" << std::endl;
        std::cerr << "\t\t --nosens: \t\t Deactivate the sensitivity calculation.\n" << std::endl;
        std::cerr << "\t\t --forcesens: \t\t Force the sensitivity calculation ON YOUR OWN RISK!\n" << std::endl;
        std::cerr << "\t\t --variance: \t\t Change calculation from overlap complexity to variance.\n" << std::endl;
        std::cerr << "\t\t --opt-oc [int]: \t Creates [int] times new patternsets and tries to optimize them to best variance/oc; after modifying by '--permut [int]'." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: opt-oc = 100\n" << std::endl;
        std::cerr << "\t\t --opt-sens [int]: \t Creates [int] times new patternsets and tries to optimize them to best sensitivity; after variance/oc optimization." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: opt-sens = 5000\n" << std::endl;
        std::cerr << "\t\t --notquiet: \t\t Show each step of the improving mode.\n" << std::endl;
        std::cerr << "\t\t --permut [int]: \t Selects [int] times a specific pattern and tries to modify it randomly by permutation." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: permut = 25000\n" << std::endl;
        std::cerr << "\t=== Variance Parameters ===" << std::endl;
        std::cerr << "\t\t -S [int]: \t\t Sequence length of the dataset." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: S = 10000\n" << std::endl;
        std::cerr << "\t\t -q [double]: \t\t Background match probability; 0 < q <= 1." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: q = 0.25\n" << std::endl;
        std::cerr << "\t\t -p [double]: \t\t Match probability for a pair of 'homologous' positions; 0 < q <= p <= 1." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: p = 0.75\n" << std::endl;
        std::cerr << "\t\t -H [int]: \t\t Length of a possible homologue region on a dataset." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: H = 64\n" << std::endl;
        std::cerr << "\t=== Pattern Parameter ====" << std::endl;
        std::cerr << "\t\t -m [int]: \t\t Number of patterns, afterwards creating autopatternset." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: m = 10\n" << std::endl;
        std::cerr << "\t\t -d [int]: \t\t Number of don't care positions, afterwards creating autopatternset." << std::endl;
        std::cerr << "\t\t    [int]-[int]: \t Min. and max. number of don't care positions, afterwards creating autopatternset." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: d = (1/p)*weight \t\t\t{variance}" << std::endl;
        std::cerr << "\t\t\t\t\t\t default: d = [weight+3]-[2*(weight+3)] \t{OC}" << std::endl;
        std::cerr << "\t\t\t\t\t\t patternlength = weight + don't care\n" << std::endl;
        std::cerr << "\t\t -w [int]: \t\t Pattern weight, afterwards creating autopatternset." << std::endl;
        std::cerr << "\t\t\t\t\t\t default: w = 8\n" << std::endl;
        std::cerr << "\t\t --pattern <File>: \t Reading pattern from <File> in pattern format with '0' and '1', seperated by ','|' '|'.'|';'|'\\n'|'\\t'.\n\n" << std::endl;
        std::cerr << "\t\t --outfile <File>: \t Save the best pattern, its variance/oc and norm_variance/oc into <File>.\n" << std::endl;
        std::cerr << "\t=== Additional Parameters ====" << std::endl;
        std::cerr << "\t\t --version: \t\t Print the program version.\n" << std::endl;
        std::cerr << "\t\t --help: \t\t Print this help.\n" << std::endl;
        return;
    }
    if (errmsg == "parsing") {
        printf("%c[1;31mError", 27);
        printf("%c[0m", 27);
        std::cerr << " while parsing " << pos / 2 + 1 << ". argument, unknown option '" << argv[pos] << "'!\n" << std::endl;
        return;
    }
    if (errmsg == "version") {
        std::cout << "rasbhari, Version 1.4.0 - (c) 2017 Lars Hahn" << std::endl;
        std::cout << "This program is released under GPLv3." << std::endl << std::endl;
        std::cout << "This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight." << std::endl;
        std::cout << "It is possible to improve your patternset and read patterns from a file." << std::endl;
    }
}