/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * pattern object header file
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
#ifndef PATTERN_HPP_
#define PATTERN_HPP_

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

/**
 * An object/instance of the pattern class represents a binary pattern that
 * consists of match positions and don't-care positions. Each pattern can have
 * a score for comparing these patterns against others.
 * Lists of don't-care positions and match positons can be requested.
 * Random patterns can be initialised depending on passed parameter 
 * configurations.
 * Iterators can be used to iterate symbol wise over a pattern.
 */
class pattern{
    public:
        pattern(uint64_t Idx = 0);
        pattern(const pattern & Pat);
        pattern(unsigned Weight, unsigned DontCare,unsigned Idx = 0);
        pattern(unsigned Weight, unsigned DontCare, uint64_t Seed, unsigned Idx = 0);
        pattern(std::string & StrPat);
        pattern(std::string && StrPat);
        pattern(std::vector<char> & StrPat);
        pattern(std::vector<unsigned char> & StrPat);

        void push_back(std::string && Chars);
        void push_back(std::string & Chars);
        void push_back(char C);
        void push_back(unsigned char C);

        std::string to_string() const;
        std::vector<unsigned> match_pos() const;
        std::vector<unsigned> dc_pos() const;
        bool is_match(unsigned Pos) const;
        void bit_swap(unsigned PosA, unsigned PosB);
        void random_swap();
        void random_swap(uint64_t Seed);

        void random(unsigned Weight, unsigned DontCare);
        void random(unsigned Weight, unsigned DontCare, uint64_t Seed);

        unsigned char operator[](size_t Idx) const;
        unsigned char & operator[](size_t Idx);
        bool operator<(const pattern & P) const;
        bool operator>(const pattern & P) const;
        bool operator==(const pattern & P) const;
        bool operator!=(const pattern & P) const;

        double score() const;
        uint64_t bits() const;
        unsigned weight() const;
        unsigned length() const;
        unsigned dontcare() const;
        unsigned get_overlap(const pattern & P, int Shift) const;
        unsigned idx() const;
        void set_score(double Scr);
        void set_idx(unsigned Idx);

        typedef std::vector<unsigned char>::iterator iterator;
        iterator begin();
        iterator end();


    private:
        void _check_pattern();

        std::vector<unsigned char> _VectorPattern;
        std::vector<unsigned> _MatchPos;
        std::vector<unsigned> _DCPos;
        double _Score;
        uint64_t _BitPattern;
        unsigned _Idx;
        bool _IsBit;
};
#endif