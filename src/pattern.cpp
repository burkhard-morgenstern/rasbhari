/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * pattern object file
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
#include "pattern.hpp"

/**
 * The empty-constructor, can be used, if symbols should be pushed into the
 * pattern instance.
 */
pattern::pattern(uint64_t Idx):  _Score(1), _BitPattern(0), _Idx(Idx), _IsBit(false){
}
/**
 * The copy-constructor, can be used, if an already existing pattern instance
 * shall be copied.
 *
 * @param Pat       The pattern instance that should be copied.
 */
pattern::pattern(const pattern & Pat){
    _VectorPattern = Pat._VectorPattern;
    _MatchPos = Pat._MatchPos;
    _DCPos = Pat._DCPos;
    _Score = Pat._Score;
    _BitPattern = Pat._BitPattern;
    _Idx = Pat._Idx;
    _IsBit = Pat._IsBit;
}
/**
 * The default-constructer; the number of don't care positions and match
 * match positions can be passed. A possible index can passed, too, if needed.
 *
 * @param Weight        The weight of a pattern, i.e. the number of '1'.
 *
 * @param DontCare      The number of don't care positions in each pattern,
 *                          i.e. the number of '0'.
 *
 * @param Index         The index, probabliy used within a set of different
 *                          patterns.
 */
pattern::pattern(unsigned Weight, unsigned DontCare, unsigned Idx): _Score(1), _BitPattern(0), _Idx(Idx),_IsBit(false){
    random(Weight, DontCare);
}
/**
 * The seed-default-constructer; the number of don't care positions and match
 * match positions can be passed. A possible index can passed, too, if needed.
 * A seed for the random generator can be passed; can be used for debugging
 * by always generating patterns from the same seed.
 *
 * @param Weight        The weight of a pattern, i.e. the number of '1'.
 *
 * @param DontCare      The number of don't care positions in each pattern,
 *                          i.e. the number of '0'.
 *
 * @param Index         The index, probabliy used within a set of different
 *                          patterns.
 *
 * @param Seed          The seed that is used for the random generator.
 */
pattern::pattern(unsigned Weight, unsigned DontCare, uint64_t Seed, unsigned Idx): _Score(1), _BitPattern(0), _Idx(Idx), _IsBit(false){
    random(Weight, DontCare, Seed);
}
/**
 * The string-constructor; the pattern is created from a string representing a
 * binary pattern.
 *
 * @param StrPat        The string pattern transformed into a pattern instance.
 */
pattern::pattern(std::string & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
    _VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
    _check_pattern();
}
/**
 * The string-constructor; the pattern is created from a string representing a
 * binary pattern.
 *
 * @param StrPat        The string pattern transformed into a pattern instance.
 */
pattern::pattern(std::string && StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
    _VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
    _check_pattern();
}
/**
 * The vector-constructor; the pattern is created from a vector representing a
 * binary pattern with characters.
 *
 * @param StrPat        The vector pattern transformed into a pattern instance.
 */
pattern::pattern(std::vector<char> & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
    _VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
    _check_pattern();
}
/**
 * The vector-constructor; the pattern is created from a vector representing a
 * binary pattern with unsigned characters.
 *
 * @param StrPat        The vector pattern transformed into a pattern instance.
 */
pattern::pattern(std::vector<unsigned char> & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
    _VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
    _check_pattern();
}


/**
 * An entire string, representing a binary pattern (part), can be used to push
 * symbols into the pattern.
 *
 * @param Chars         String containing symbols to be appended to the pattern.
 */
void pattern::push_back(std::string && Chars){
    push_back(Chars);
}
/**
 * An entire string, representing a binary pattern (part), can be used to push
 * symbols into the pattern.
 *
 * @param Chars         String containing symbols to be appended to the pattern.
 */
void pattern::push_back(std::string & Chars){
    for(unsigned char C : Chars){
        push_back(C);
    }
}

/**
 * Symbols can be pushed onto the pattern. These symbols will be automatically
 * transformed into the correct representation.
 *
 * @param Char          A symbol, that should be appended to the pattern.
 */
void pattern::push_back(char C){
    push_back((unsigned char)C);
}

/**
 * Symbols can be pushed onto the pattern. These symbols will be automatically
 * transformed into the correct representation.
 *
 * @param Char          A symbol, that should be appended to the pattern.
 */
void pattern::push_back(unsigned char C){
    switch(C){
        case '1':
        case 'X':
        case 'x':
        case '#':
        case '*':
            C = 1;
            break;
        case '0':
        case 'O':
        case 'o':
        case '-':
            C = 0;
            break;
    }
    if(C != 1 && C != 0){
        std::cerr << "Illegal character in pattern!\nFormat: match = {1,X,x,#,*,} | don't care = {0,O,o,-,}" << std::endl;
        std::exit(-1);
    }
    _VectorPattern.push_back(C);
    if(C == 1){
        _MatchPos.push_back(_VectorPattern.size()-1);
    }
    else{
        _DCPos.push_back(_VectorPattern.size()-1);
    }
    _BitPattern <<= 1;
    _BitPattern |= C;
}

/**
 * Checks, if the pattern is in correct format. This is used, if an entire
 * vector or string is used in the constuctor.
 */
void pattern::_check_pattern(){
    unsigned Ctr = 0;

    if(_VectorPattern.size() < 64){
        _IsBit = true;
    }

    for(auto & C : _VectorPattern){
        switch(C){
            case '1':
            case 'X':
            case 'x':
            case '#':
            case '*':
                C = 1;
                break;
            case '0':
            case 'O':
            case 'o':
            case '-':
                C = 0;
                break;
        }
        if(C != 1 && C != 0){
            std::cerr << "Illegal character in pattern!\nFormat: match = {1,X,x,#,*,} | don't care = {0,O,o,-,}" << std::endl;
            std::exit(-1);
        }
        if(C == 1){
            _MatchPos.push_back(Ctr);
        }
        else{
            _DCPos.push_back(Ctr);
        }
        _BitPattern <<= 1;
        _BitPattern |= C;
        Ctr++;
    }
}

/**
 * Creates a random pattern from passed pattern parameters.
 *
 * @param Weight        The weight of the pattern, i.e. the number of '1's
 *
 * @param DontCare      The number of don't care positions, i.e. the number of
 *                          '0's.
 */
void pattern::random(unsigned Weight, unsigned DontCare){
    std::random_device RandomBit;
    random(Weight, DontCare, RandomBit());
}
/**
 * Creates a random pattern from passed pattern parameters using a passed
 * integer seed for the random generator.
 *
 * @param Weight        The weight of the pattern, i.e. the number of '1's
 *
 * @param DontCare      The number of don't care positions, i.e. the number of
 *                          '0's.
 *
 * @param Seed          The seed used for the random generator to create a
 *                          random intialised pattern instance.
 */
void pattern::random(unsigned Weight, unsigned DontCare, uint64_t Seed){
    unsigned Length = DontCare + Weight;
    if(Length < 64){
        _IsBit = true;
    }
    if(Weight < 2){
        std::cerr << "Illegal value for weight!\nMinimum allowed weight is 2!" << std::endl;
        std::exit(-1);
    }
    std::default_random_engine BitRandom(Seed);
    std::mt19937 BitGenerator(BitRandom());
    std::uniform_int_distribution<unsigned> BitPosition(1, Length-2);

    _MatchPos.clear();
    _DCPos.clear();
    _BitPattern = 0;
    _VectorPattern = std::vector<unsigned char>(Length,0);
    
    _VectorPattern[0] = 1;
    _VectorPattern[Length-1] = 1;
    _BitPattern = ((uint64_t) 1 << (Length-1)) | (uint64_t)1;
    Weight -= 2;

    while(Weight > 0){
        unsigned Pos = BitPosition(BitGenerator);
        if(_VectorPattern[Pos] == 0){
            Weight--;
            _VectorPattern[Pos] = 1;
            _BitPattern |= (uint64_t)1 << (Length -1 - Pos);
        }
    }
    for(unsigned i = 0; i < _VectorPattern.size(); i++){
        if(_VectorPattern[i] == 1){
            _MatchPos.push_back(i);
        }
        else{
            _DCPos.push_back(i);
        }
    }
}

/**
 * Transforms a pattern instance into a human readable string format.
 *
 * @return              The string representation of the pattern instance.
 */
std::string pattern::to_string() const{
    std::string Pattern(_VectorPattern.begin(), _VectorPattern.end());
    std::transform(Pattern.begin(), Pattern.end(), Pattern.begin(),[](unsigned char C) { return (C+'0'); });
    return Pattern;
}

/**
 * Shows, if the Pos-th position in the pattern is a match position, i.e.
 * if it is a '1'.
 *
 * @param Pos           Requested position to be checked for being a match.
 *
 * @return              True, if the requested position is a match position.
 */
bool pattern::is_match(unsigned Pos) const{
    return _VectorPattern[Pos] == 1;
}

/**
 * Returns a list containing indices, where a match position is found.
 *
 * @return              Integer list of all match positions.
 */
std::vector<unsigned> pattern::match_pos() const{
    return _MatchPos;
}
/**
 * Returns a list containing indices, where a don't-care position is found.
 *
 * @return              Integer list of all don't care positions.
 */
std::vector<unsigned> pattern::dc_pos() const{
    return _DCPos;
}

/**
 * Performs a bit swap for two distinct positions.
 * One position is a match position, the other one is a don't care position.
 * 
 * @param PosA          The first position to be swapped.
 *
 * @param PosB          The second position to be swapped.
 */
void pattern::bit_swap(unsigned PosA, unsigned PosB){
    std::swap(_VectorPattern[PosA],_VectorPattern[PosB]);
    if(_IsBit && _VectorPattern[PosA] !=_VectorPattern[PosB]){
        _BitPattern ^= ((uint64_t)1 << (_VectorPattern.size()-PosA-1)) | ((uint64_t)1 << (_VectorPattern.size()-PosB-1));
        if(_VectorPattern[PosA] != 0){
            std::swap(PosA,PosB);
        }
        _MatchPos.erase(std::find(_MatchPos.begin(), _MatchPos.end(), PosA));
        _DCPos.erase(std::find(_DCPos.begin(), _DCPos.end(), PosB));
        _MatchPos.push_back(PosB);
        _DCPos.push_back(PosA);
        std::sort(_MatchPos.begin(), _MatchPos.end());
        std::sort(_DCPos.begin(), _DCPos.end());
    }
}

/**
 * Performs a random bit swap for two distinct positions.
 * One position is a match position, the other one is a don't care position.
 */
void pattern::random_swap(){
    std::random_device RandomBit;
    random_swap(RandomBit());
}

/**
 * Performs a random bit swap for two distinct positions using a passed seed
 * for the random generator used for determining the positions.
 * One position is a match position, the other one is a don't care position.
 */
void pattern::random_swap(uint64_t Seed){
    if(_DCPos.size() != 0 && _MatchPos.size() > 2){
        std::default_random_engine BitRandom(Seed);
        std::mt19937 BitGenerator(BitRandom());
        std::uniform_int_distribution<unsigned> SwapPosition(1, _MatchPos.size()-2);
        unsigned MatchPosition = _MatchPos[SwapPosition(BitGenerator)];
        SwapPosition = std::uniform_int_distribution<unsigned>(0, _DCPos.size()-1);
        unsigned DCPosition = _DCPos[SwapPosition(BitGenerator)];
        bit_swap(MatchPosition, DCPosition);
    }
}


/**
 * Random access operator, that returns the symbol from a specific
 * position in the pattern.
 *
 * @param Idx           The index number of a symbol, that sould be returned.
 *
 * @return              The symbol from the Idx-th position; constant
 */
unsigned char pattern::operator[](size_t Idx) const{
    return _VectorPattern[Idx];
}
/**
 * Random access operator, that returns the symbol from a specific
 * position in the pattern.
 *
 * @param Idx           The index number of a symbol, that sould be returned.
 *
 * @return              The symbol from the Idx-th position; r-value
 */
unsigned char & pattern::operator[](size_t Idx){
    return _VectorPattern[Idx];
}
/**
 * Compares the score of the pattern with another one and returns true,
 * if this pattern has the lower score.
 *
 * @param P             The pattern that is used for score comparison
 *
 * @return              True, if this pattern has a lower score
 */
bool pattern::operator<(const pattern & P) const{
    return _Score < P._Score;
}
/**
 * Compares the score of the pattern with another one and returns true,
 * if this pattern has the higher score.
 *
 * @param P             The pattern that is used for score comparison
 *
 * @return              True, if this pattern has a higher score
 */
bool pattern::operator>(const pattern & P) const{
    return _Score > P._Score;
}
/**
 * Compares another pattern with this pattern if they are identical.
 *
 * @param P             The pattern that is used for comparison.
 *
 * @return              True, if the other pattern is identical to this one.
 */
bool pattern::operator==(const pattern & P) const{
    if(_IsBit){
        return _BitPattern == P._BitPattern;
    }
    return _VectorPattern == P._VectorPattern;
}
/**
 * Compares another pattern with this pattern if they are not identical.
 *
 * @param P             The pattern that is used for comparison.
 *
 * @return              True, if the other pattern is not identical to this one.
 */
bool pattern::operator!=(const pattern & P) const{
    if(_IsBit){
        return _BitPattern != P._BitPattern;
    }
    return _VectorPattern != P._VectorPattern;
}

/**
 * Returns the score of the pattern, that was set before.
 *
 * @return              Double value, the pattern score.
 */
double pattern::score() const{
    return _Score;
}
/**
 * Returns the bit-integer representation of the pattern.
 * Available only, if the length is short than 64 positions.
 *
 * @return              Bit-Integer pattern representation
 */
uint64_t pattern::bits() const{
    return _BitPattern;
}
/**
 * Returns the weight of the pattern
 *
 * @return              The pattern weight, i.e. the number of '1's.
 */
unsigned pattern::weight() const{
    return _MatchPos.size();
}
/**
 * Returns the lenght of the pattern
 *
 * @return              The pattern length, the sum of weight and DC positions.
 */
unsigned pattern::length() const{
    return _VectorPattern.size();
}
/**
 * Returns the number of don't care positons of the pattern
 *
 * @return              The number of don't care positions, i.e. the number of
 *                          '1's.
 */
unsigned pattern::dontcare() const{
    return _VectorPattern.size() - _MatchPos.size();
}
/**
 * Returns the index of the pattern, that was set before.
 *
 * @return              Index value, that was set before.
 */
unsigned pattern::idx() const{
    return _Idx;
}
/**
 * Sets the score of the pattern, probably used in a pattern set.
 *
 * @param Scr           Double value, the pattern score.
 */
void pattern::set_score(double Scr){
    _Score = Scr;
}
/**
 * Sets the index for the pattern, probably used in a pattern set.
 *
 * @Param Idx           Index value for the pattern
 */
void pattern::set_idx(unsigned Idx){
    _Idx = Idx;
}


/**
 * Returns the overlap for this pattern with another pattern, i.e. the
 * number of shared match positions for a fixed shift.
 *
 * E.g:  p1:   10011101
 *       p2:     11100101
 *                ++  +
 * The overlap is three for a shift=2
 *
 * @param P             The pattern used for calculating the overlap
 *
 * @param OldShift      The number of shifts performed on the pattern.
 */
unsigned pattern::get_overlap(const pattern & P, int OldShift) const{
    unsigned Overlap = 0;
    int Shift =  (int)length() - ((int)P.length()+OldShift);
    if(P._IsBit && _IsBit){
        uint64_t ShiftPattern = bits();
        uint64_t BitPat = P.bits();
        if(Shift < 0){
            std::swap(ShiftPattern, BitPat);
            Shift = std::abs(Shift);
        }
        ShiftPattern >>= Shift;
        ShiftPattern &= BitPat;
        unsigned ShiftSize = std::log2(ShiftPattern)+1;
        for(unsigned i = 0; i < ShiftSize; i++){
            Overlap += (ShiftPattern & (uint64_t)1);
            ShiftPattern >>= 1;
        }
    }
    else{
        auto PatAIt = _VectorPattern.begin();
        auto PatBIt = P._VectorPattern.begin();
        Shift += (int)P._VectorPattern.size()-(int)_VectorPattern.size();
        if(Shift > 0){
            PatBIt += Shift;
        }
        else{
            PatAIt -= Shift;
        }
        for(;PatAIt != _VectorPattern.end() && PatBIt != P._VectorPattern.end(); PatBIt++, PatAIt++){
            Overlap += *PatAIt * *PatBIt;
        }
    }
    return Overlap;
}


/**
 * Begin iterator that can be used to iterate over the pattern.
 *
 * @return              The begin iterator of the pattern symbol vector.
 */
pattern::iterator pattern::begin(){
    return _VectorPattern.begin();
}
/**
 * End iterator that can be used to iterate over the pattern.
 *
 * @return              The end iterator of the pattern symbol vector.
 */
pattern::iterator pattern::end(){
    return _VectorPattern.end();
}