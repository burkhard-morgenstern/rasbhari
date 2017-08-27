/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * sensitivity memory observation namespace header file
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
#ifndef SENSMEM_HPP_
#define SENSMEM_HPP_

#include <cstdint>
#include <iostream>
#include <string>
#include "sys/sysinfo.h"
#include "rasbopt.hpp"

/**
 * Sensitivity memory observation namespace; it checks if there is enough
 * free memory for the sensitivity calculation without entering possible swap.
 * If there is not enough memory, a user interaction is started once.
 */
namespace sensitivity_memory{
    bool check_memory();
    double available_memory();
    void security_message(std::string errmsg);
    extern double MemVal;
};

#endif