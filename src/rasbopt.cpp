/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari options namespace file
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
#include "rasbopt.hpp"

/**
 * Contains all necessary parameters for the rasbhari pattern sets.
 * These values are specific only for the rasbhari patterns.
 * Normal patterns or patternsets do not depend on these parameters.
 */
namespace rasb_opt{
    std::string OutFile = "rasbhari_set.pat";
    std::string InFile;
    double P = 0.75;
    double Q = 0.25;
    long int Seed;
    unsigned Size = 10;
    unsigned Weight = 8;
    unsigned SeqLength = 10000;
    unsigned Limit = 25000;
    unsigned H = 64;
    unsigned OptOc = 100;
    unsigned OptSens = 5000;
    unsigned MinDontcare = 10;
    unsigned MaxDontcare = 10;
    unsigned ImproveMode = 0;
    bool Improve = false;
    bool Forcesens = false;
    bool Oc = true;
    bool Quiet = true;
    bool Sens = false;
    bool SetSeed = false;
    bool Silent = true;

    /**
     * Parses minimal and maximal length from a character array.
     *
     * @param Str       Character array containing maximal and minimal
     *                       DC-positions in form of 'a' or 'a-b'.       
     */
    void parse_length(const char* Str){
        std::string Length = std::string(Str);
        unsigned Pos = (unsigned) Length.size();
        unsigned Leng = Pos;

        for(unsigned i = 0; i < Leng; i++){
            if(std::isdigit(Length[i]) == false){
                Pos = i;    
            }
        }

        if(Pos < Leng){
            char *Tmp1, *Tmp2;
            Tmp1 = new char[Pos];
            Tmp2 = new char[Leng-Pos-1];
            for(uint32_t i = 0; i < Leng; i++){
                if(i < Pos){
                    Tmp1[i] = Length[i];
                }
                if(i > Pos){
                    Tmp2[i-Pos-1] = Length[i];
                }
            }
            MinDontcare = std::atoi(Tmp1);
            MaxDontcare = std::atoi(Tmp2);
            delete Tmp1;
            delete Tmp2;
            
        }
        else{
            MinDontcare = std::atoi(Str);
            MaxDontcare = std::atoi(Str);
        }
    }
}