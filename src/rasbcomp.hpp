/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari computation namespace header file
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
#ifndef RASBCOMP_HPP_
#define RASBCOMP_HPP_

#include <cstdint>
#include <string>
#include "pattern.hpp"
#include "patternset.hpp"
#include "rasbopt.hpp"
#include "speedsens.hpp"
#include "sensmem.hpp"

/**
 * The rasbhari_compute namespace holds the actual computations for the
 * OC/variance optimisation and sensitivity reduction. This is the
 * implementation of the basic algorithms for overlap complexity, introduced
 * by Ilie&Ilie, and for the variance, introduced in the paper mentioned above
 * by B. Morgenstern and L. Hahn. 
 */
namespace rasbhari_compute{
    double pair_coef_var(pattern & Pat1, pattern & Pat2);
    double pair_coef_var_sym(pattern & Pat1, pattern & Pat2);
    double pair_coef_oc(pattern & Pat1, pattern & Pat2);
    double sensitivity(patternset & Pattern);
    extern bool DoOnce;
};
#endif