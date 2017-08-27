/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * sensitivity computation namespace header file
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
#ifndef SPEEDSENS_HPP_
#define SPEEDSENS_HPP_

#include <cstring>

#include "sensmem.hpp"
#include "rasbopt.hpp"

/**
 * This namespace holds the sensitivity caluclation implemented by Ilie&Ilie.
 * The functions in the sensitivity computation namespace file, speedsens.cpp,
 * are directly extracted from their spaced seed optimisation tool, SpEED.
 * Only commands for calculating the currently used memory and some user
 * interacting points have been added.
 */
namespace speedsens{
    inline long long BIN_REVERSED_TO_INT2(char *s);
    double MULTIPLE_SENSITIVITY2(char** SEEDS, int NO_SEEDS, long long N, double P);
};
#endif