/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari object header file
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
#ifndef RASBHARI_HPP_
#define RASBHARI_HPP_

#include <ios>
#include <fstream>
#include "patternset.hpp"
#include "pattern.hpp"
#include "rasbcomp.hpp"
#include "rasbopt.hpp"

/**
 * An object/instance of the rasbhari class represents and contains depending
 * on the actions an optimised patternsets with its OC/Variance score, the
 * sensitivity and other necessary information.
 * Iterators can be used to iterate over the rasbhari patternset.
 */
class rasbhari{
    public:
        rasbhari();
        rasbhari(const rasbhari &RasbObj);
        rasbhari(unsigned Size, unsigned Weight, unsigned DontCare);
        rasbhari(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
        rasbhari(std::string PatternFile);

        void push_back(pattern & Pat);
        void push_back(pattern && Pat);

        void calculate();
        double calculate_pair(size_t Idx1, size_t Idx2);
        void calculate_sensitivity();
        void update(unsigned Idx);
        bool climb_hill();
        bool hill_climbing(unsigned Limit = rasb_opt::Limit);
        bool iterate_hill_climbing(unsigned Limit = rasb_opt::Limit, unsigned Iteration = rasb_opt::OptOc);
        bool climb_hill_sensitivity(unsigned Limit = rasb_opt::Limit, unsigned Iteration = rasb_opt::OptOc, bool InitSens = true);
        bool hill_climbing_sensitivity(unsigned Limit = rasb_opt::Limit, unsigned Iteration = rasb_opt::OptOc, unsigned Loop = rasb_opt::OptSens);

        patternset pattern_set() const;
        patternset & pattern_set();
        double score() const;
        double sensitivity() const;
        unsigned size() const;
        void print();
        void to_file(std::string OutFile = rasb_opt::OutFile);

        pattern operator[](size_t Idx) const;
        pattern & operator[](size_t Idx);
        typedef std::vector<pattern>::iterator iterator;
        iterator begin();
        iterator end();

    private:
        void _debug();
        void _make_pattern_list();
        void _adjust_coef_mat();
        void _check_pattern_number(unsigned & Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
        double _max_pat_no(unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare);
        double _binom_coef(unsigned n, unsigned k);

        std::vector< std::vector<double> > _CoefMat;
        std::vector<pattern> _PatternList;
        patternset _RasbhariPattern;
        double _RasbhariScore;
        double _RasbhariSensitivity;
        unsigned long _PatNo;
};

#endif