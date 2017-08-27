/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * sensitivity memory observation namespace file
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
 */
#include "sensmem.hpp"

double sensitivity_memory::MemVal = 0;

/**
 * A function, which organizes the checking of the memory usage and interacts
 * with the user, if a lot of memory is currently used.
 * 
 * @return  Possible user decision, if the calculation has to be stop or not.
 */
bool sensitivity_memory::check_memory(){
    std::string Input;
    uint32_t Stry;
    bool Answer, InputErr;
    MemVal = available_memory();
    Answer = false;
    Stry = 0;
    InputErr = false;
    if(MemVal< 0.25){
        security_message("lowmemory");
        while(Stry < 3){
            Stry++;
            InputErr = false;
            std::cin >> Input;
            switch(Input[0]){
                case 'Y':
                case 'y':
                    rasb_opt::Forcesens = true;
                    Stry = 4;
                    Answer = false;
                    break;
                case 'N':
                case 'n':
                    Stry = 4;
                    Answer = true;
                    break;
                default:
                    security_message("inputconf");
                    Answer = true;
                    InputErr = true;
                    break;
            }
        }
    }

    if(InputErr && Stry < 4){
        security_message("inputerr");
    }
    return Answer;
}


/**
 * Calculates the amount of available memory in total and returns the amount as
 * percentage.
 *
 * @return  The percentage of available memory in total.
 */
double sensitivity_memory::available_memory(){
    struct sysinfo MemInfo;
    double IsAvailable;
    uint64_t RamTotal, RamFree;

    sysinfo (&MemInfo);

    RamTotal = MemInfo.totalram;
    RamFree = MemInfo.freeram;
    RamTotal *= MemInfo.mem_unit;
    RamFree *= MemInfo.mem_unit;

    IsAvailable = (double)RamFree/(double)RamTotal;
    return IsAvailable;
}


/**
 * Prints warnings or error messages to std::cerr if some problems were detected
 * during the sensitivity calculation. 
 *
 * @param errmsg        Indicator, which error massage has to be thrown.
 */
void sensitivity_memory::security_message(std::string errmsg){
    if (errmsg == "bitmode") {
        printf("%c[1;33m", 27);
        std::cerr << "A patternlength is over 63, leaving bitmode ..." << std::endl;
        std::cerr << "Using your pattern conditions it is not possible to calculate the sensitivity!" << std::endl;
        std::cerr << "Deactivating sensitivity calculation\n" << std::endl;
        printf("%c[0m", 27);
        return;
    }
    else if (errmsg == "lowmemory"){
        std::cerr << "\rSensitivity Calculation consumes a lot of memory (>75% RAM is used in total). If you want to continue, ";
        printf("%c[1;33m", 27);
        std::cerr << "PLEASE SAVE ALL RELEVANT DATA!";
        printf("%c[0m", 27);
        std::cout << std::endl;
        std::cout << std::endl << "Do you want to continue calculating the sensitivity on your own risk (slow system with swap)? [Y/N]\t";
    }
    else if (errmsg == "inputconf"){
        std::cerr << "Please enter 'Y' to continue with the sensitivity calculating or 'N' to abort the calculation!" << std::endl;
    }
    else if (errmsg == "inputerr"){
        std::cerr << "Could not read answer in three attempts, abort sensitivity calculation!" << std::endl;
    }
    else if (errmsg == "memerror"){
        std::cerr << "There is not enough memory! Aborting sensitivity calculation!" << std::endl;
    }
}