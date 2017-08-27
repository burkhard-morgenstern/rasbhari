/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari options namespace header file
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
#ifndef RASBOPT_HPP_
#define RASBOPT_HPP_

#include <cstdlib>
#include <string>

/**
 * Contains all necessary parameters for the rasbhari pattern sets.
 * These values are specific only for the rasbhari patterns.
 * Normal patterns or patternsets do not depend on these parameters.
 */
namespace rasb_opt{
    extern std::string OutFile;
    extern std::string InFile;
    extern double P;
    extern double Q;
    extern long int Seed;
    extern unsigned Size;
    extern unsigned Weight;
    extern unsigned MinDontcare;
    extern unsigned MaxDontcare;
    extern unsigned SeqLength;
    extern unsigned Limit;
    extern unsigned H;
    extern unsigned OptOc;
    extern unsigned OptSens;
    extern unsigned ImproveMode;
    extern bool Improve;
    extern bool Forcesens;
    extern bool Oc; 
    extern bool Quiet;
    extern bool Sens;
    extern bool SetSeed;
    extern bool Silent;

    void parse_length(const char* Str);
};
#endif
