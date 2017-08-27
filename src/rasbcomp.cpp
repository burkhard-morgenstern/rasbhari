/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari computation namespace file
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
#include "rasbcomp.hpp"

/**
 * Variable, used for sensitivity calculation. If there is a hugh memory usage
 * and the users wishes to continue, there wont be any question about the
 * memory usage in this run.
 */
bool rasbhari_compute::DoOnce = false;

/**
 * Calculates for a pair of pattern the variance Var(N). Here it is Pat1 != Pat2
 * and it is 2*Weight-Overlap = n(P1,P2,s). For theory, please have a look at 
 * rasbhari, L. Hahn et al.
 *
 * @param Pat1          The first (upper) pattern.
 *
 * @param Pat2          The second pattern, which is shifted.
 *
 * @return              The pairwise value of Var(N)
 */
double rasbhari_compute::pair_coef_var(pattern & Pat1, pattern & Pat2){
    double CoEf = 0;
    double LengthMean = (Pat1.length() + 1 + Pat2.length()) / 2;
    for(int s = -(int)Pat2.length()+1; s < (int)Pat1.length(); s++){
        unsigned Overlap = Pat1.get_overlap(Pat2,s);
        CoEf += (rasb_opt::SeqLength - LengthMean + 1)*(pow(rasb_opt::P, 2*rasb_opt::Weight-Overlap) - pow(rasb_opt::P, 2 * rasb_opt::Weight))+(rasb_opt::SeqLength - LengthMean + 1)*(rasb_opt::SeqLength - LengthMean)*(pow(rasb_opt::Q, 2*rasb_opt::Weight-Overlap) - pow(rasb_opt::Q, 2 * rasb_opt::Weight));
    }
    return CoEf;
}
/**
 * Calculates for a pair of pattern the variance Var(N). Here it is Pat1 == Pat2
 * and it is 2*Weight-Overlap = n(P1,P2,s). For theory, please have a look at 
 * rasbhari, L. Hahn et al.
 *
 * @param Pat1          The first (upper) pattern.
 *
 * @param Pat2          The second pattern, which is shifted.
 *
 * @return              The pairwise value of Var(N)
 */
double rasbhari_compute::pair_coef_var_sym(pattern & Pat1, pattern & Pat2){
    double CoEf = 0;
    double LengthMean = (Pat1.length() + 1 + Pat2.length()) / 2;
    for(int s = 0; s < (int)Pat1.length(); s++){
        unsigned Overlap = Pat1.get_overlap(Pat2,s);
        CoEf += (rasb_opt::SeqLength - LengthMean + 1)*(pow(rasb_opt::P, 2*rasb_opt::Weight-Overlap) - pow(rasb_opt::P, 2 * rasb_opt::Weight))+(rasb_opt::SeqLength - LengthMean + 1)*(rasb_opt::SeqLength - LengthMean)*(pow(rasb_opt::Q, 2*rasb_opt::Weight-Overlap) - pow(rasb_opt::Q, 2 * rasb_opt::Weight));
    }
    return CoEf;
}

/**
 * Calculates for a pair of pattern the overlap complexity.
 * For theory, please have a look at rasbhari, L. Hahn et al.
 *
 * @param Pat1          The first (upper) pattern.
 *
 * @param Pat2          The second pattern, which is shifted.
 *
 * @return              The pairwise value of Var(N)
 */
double rasbhari_compute::pair_coef_oc(pattern & Pat1, pattern & Pat2){
    double CoEf = 0;
    for(int s = -(int)Pat2.length()+1; s < (int)Pat1.length(); s++){
        CoEf += (uint64_t)1 << Pat1.get_overlap(Pat2,s);
    }
    return (double) CoEf;
}


/**
 * Interface function to create the right pattern format for the speed functions.
 *
 * @param Pattern       Contains the patternset, for which the sensitivity
 *                          has to be calculated.
 *
 * @return              Sensitivity of the patterset.
 */
double rasbhari_compute::sensitivity(patternset & Pattern){
    if(Pattern.max_weight() > 63){
        sensitivity_memory::security_message("bitmode");
        rasb_opt::Sens = false;
        rasb_opt::OptSens = 1;
        if(rasb_opt::OptOc > 1){
            rasb_opt::ImproveMode = 2;
        }
        else{
            rasb_opt::ImproveMode = 1;
        }
        return -1;
    }
    std::vector<char> Clean;
    double SensVal = -1;
    unsigned PSize = Pattern.size();
    if(!DoOnce){
        std::cout << "\rCalculating sensitivity ...";
        std::cout.flush();
    }
    char** Pats = new char*[PSize];   //set of seeds
    for (unsigned i = 0; i < PSize; i++) {
        unsigned Length = Pattern[i].length();
        Pats[i] = new char[Length + 1];
    }

    for (unsigned i = 0; i < PSize; i++) {
        unsigned Length = Pattern[i].length();
        for (unsigned j = 0; j < Length; j++) {
            if (Pattern[i][j] == 1) {
                Pats[i][j] = '1';
            }
            else {
                Pats[i][j] = '0';
            }
        }
        Pats[i][Length] = '\0';
    }
    try{
        SensVal = speedsens::MULTIPLE_SENSITIVITY2(Pats, rasb_opt::Size, rasb_opt::H, rasb_opt::P);
        if(!DoOnce){
            std::cout << "\r" << std::string(80,' ') << "\r";
            std::cout.flush();
            DoOnce = true;
        }
    }
    catch(std::bad_alloc){
        if(!DoOnce){
            std::cout << "\r" << std::string(80,' ') << "\r";
            std::cout.flush();
            DoOnce = true;
        }
        std::cout << std::endl;
        sensitivity_memory::security_message("memerror");
        rasb_opt::Sens = false;
        rasb_opt::OptSens = 1;
        SensVal = -1;
    }
    return SensVal;
}