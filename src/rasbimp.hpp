/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari implementation namespace header file
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
#ifndef RASBIMP_HPP_
#define RASBIMP_HPP_

#include "rasbhari.hpp"
#include "rasbopt.hpp"

/**
 * The implementation namespace comprises all necessary functions in a way of an
 * 'interface', such that the user does only have to specify some parameters
 * and call the needed function, so that the final rasbhari pattern set will be
 * returned with its OC/Variance for given parameters. The patternset can be 
 * extracted. 
 */
namespace rasb_implement{
    rasbhari hillclimb_var(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    rasbhari hillclimb_oc(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    rasbhari hillclimb_var_iterative(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    rasbhari hillclimb_oc_iterative(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    rasbhari hillclimb_sens_var(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    rasbhari hillclimb_sens_oc(unsigned Size = rasb_opt::Size, unsigned Weight = rasb_opt::Weight, unsigned MinDontCare = rasb_opt::MinDontcare, unsigned MaxDontCare = rasb_opt::MaxDontcare);
    
    rasbhari _hillclimb(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
    rasbhari _hillclimb_iterative(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
    rasbhari _hillclimb_sens(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
};

#endif