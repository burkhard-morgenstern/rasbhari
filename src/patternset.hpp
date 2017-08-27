/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * pattern set object header file
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
#ifndef PATTERNSET_HPP_
#define PATTERNSET_HPP_

#include <algorithm>
#include <assert.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include "pattern.hpp"

/**
 * An object/instance of the pattern set class represents a set of binary
 * pattern. This set is dynamically; patterns can be pushed into the set.
 * Necessary parameter like mini/maximal weight / dc positions can be returned.
 *
 * Special operators and iterators are defined for an easier usage of a pattern
 * set. Patterns pushed into the set can either be patterns or strings or vectors
 * of characters.
 *
 * In addition, only important pattern set parameters can be passed and a random
 * pattern set will be created.
 */
class patternset{
    public:
        patternset();
        patternset(const patternset & PatSet);
        patternset(unsigned Size, unsigned Weight, unsigned DontCare, bool Uniq = false);
        patternset(unsigned Size, unsigned MinW, unsigned MaxW, unsigned MinDc, unsigned MaxDc, bool Uniq = false);
        patternset(std::string InputFile);

        void push_back(pattern & Pat);
        void push_back(pattern && Pat);
        void push_back(std::string & Pat);
        void push_back(std::string && Pat);
        void push_back(std::vector<char> & Pat);
        void push_back(std::vector<unsigned char> & Pat);
        void random(unsigned Size, unsigned Length, unsigned Weight, bool Uniq = false);
        void random(unsigned Size, unsigned MinLength, unsigned MaxLength, unsigned Weight, bool Uniq = false);
        void random(unsigned Size, unsigned MinLength, unsigned MaxLength, unsigned MinWeight, unsigned MaxWeight, bool Uniq = false);
        void sort();
        bool is_uniq(const pattern & Pat) const;
        void random_swap_uniq(unsigned Idx);
        void to_file(std::string OutFile = "patternset.pat");

        pattern operator[](size_t Idx) const;
        pattern & operator[](size_t Idx);
        bool operator<(const patternset & P) const;
        bool operator>(const patternset & P) const;

        unsigned max_weight() const;
        unsigned min_weight() const;
        unsigned weight() const;
        unsigned max_dontcare() const;
        unsigned min_dontcare() const;
        unsigned dontcare() const;
        unsigned max_length() const;
        unsigned min_length() const;
        unsigned length() const;
        unsigned size() const;

        double score() const;
        void set_score(double Score);


        typedef std::vector<pattern>::iterator iterator;
        iterator begin();
        iterator end();

    private:
        template<typename T>
        void _adjust(T & MinValue, T & MaxValue);
        std::vector<std::string> _split(std::string &S);
        std::vector<pattern> _PatternSet;
        double _Score;
        unsigned _MinWeight;
        unsigned _MaxWeight;
        unsigned _MinDontCare;
        unsigned _MaxDontCare;
};
template<typename T>
void patternset::_adjust(T & MinValue, T & MaxValue){
    if(MaxValue < MinValue){
        std::swap(MaxValue,MinValue);
    }
}
#endif