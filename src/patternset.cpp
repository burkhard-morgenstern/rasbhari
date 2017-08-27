/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * pattern set object file
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
#include "patternset.hpp"

/**
 * The empty-constructor, can be used, if patterns should be pushed into the
 * pattern set instance.
 */
patternset::patternset():_Score(1),_MinWeight(std::numeric_limits<unsigned>::max()),_MaxWeight(0),
    _MinDontCare(std::numeric_limits<unsigned>::max()),_MaxDontCare(0){    
}
/**
 * The copy-constructor, can be used, if an already existing pattern set instance
 * shall be copied.
 *
 * @param PatSet       The pattern set instance that should be copied.
 */
patternset::patternset(const patternset & PatSet){
    _PatternSet = PatSet._PatternSet;
    _Score = PatSet._Score;
    _MinWeight = PatSet._MinWeight;
    _MaxWeight = PatSet._MaxWeight;
    _MinDontCare = PatSet._MinDontCare;
    _MaxDontCare = PatSet._MaxDontCare;
}
/**
 * The short-default-constructer; here all patterns have the same length and the
 * same weight. A random pattern set will be created.
 *
 * @param Size          The number of patterns.
 *
 * @param Weight        The weight of each pattern.
 *
 * @param DontCare      The number of don't care positions in each pattern.
 *
 * @param Uniq          Boolean; true, if each pattern has to be a unique
 *                          pattern.
 */
patternset::patternset(unsigned Size, unsigned Weight, unsigned DontCare, bool Uniq):
    _Score(1),_MinWeight(Weight),_MaxWeight(Weight),_MinDontCare(DontCare),
    _MaxDontCare(DontCare){
    random(Size, Weight, DontCare);
}
/**
 * The long-default-constructer; here patterns can have different lengths and
 * weights. These values are evenly distributed if minimal and maximum values
 * are different.
 *
 * @param Size          The number of patterns.
 *
 * @param MinW          The minimal weight for a pattern.
 *
 * @param MaxW          The maximal weight for a pattern.
 *
 * @param MinD          The minimal number of don't care positions for patterns.
 *
 * @param MaxD          The maximal number of don't care positions for patterns.
 *
 * @param Uniq          Boolean; true, if each pattern has to be a unique
 *                          pattern.
 */
patternset::patternset(unsigned Size, unsigned MinW, unsigned MaxW, 
    unsigned MinD, unsigned MaxD, bool Uniq):_Score(1),_MinWeight(MinW),
    _MaxWeight(MaxW),_MinDontCare(MinD),_MaxDontCare(MaxD){
    _adjust(_MinWeight,_MaxWeight);
    _adjust(_MinDontCare,_MaxDontCare);
    random(Size, MaxW, MinW, MaxD, MinD);
}
/**
 * The file-constructor, can be used, if the patterns should be read from a
 * given file.
 *
 * @param InputFile   String, containing the name for a file containing a
 *                          pattern set.
 */
patternset::patternset(std::string InputFile):_Score(1),_MinWeight(std::numeric_limits<unsigned>::max()),_MaxWeight(0),
    _MinDontCare(std::numeric_limits<unsigned>::max()),_MaxDontCare(0){
    std::ifstream Input(InputFile);
    assert(Input.is_open());
    while(!Input.eof()){
        std::string Line;
        std::getline(Input,Line);
        if(Line[0] != '#'){
            std::vector<std::string> Splits = _split(Line);
            for(auto & Str : Splits){
                push_back(Str);
            }
        }
    }
    Input.close();
}


/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 *
 * @param Pat           The pattern pushed into the set of patterns; r-value.
 */
void patternset::push_back(pattern & Pat){
    _MaxDontCare = std::max(_MaxDontCare, Pat.dontcare());
    _MinDontCare = std::min(_MinDontCare, Pat.dontcare());
    _MaxWeight = std::max(_MaxWeight, Pat.weight());
    _MinWeight = std::min(_MinWeight, Pat.weight());
    _PatternSet.push_back(Pat);
}
/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 *
 * @param Pat           The pattern pushed into the set of patterns; l-value.
 */
void patternset::push_back(pattern && Pat){
    push_back(Pat);
}
/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 * A pattern instance is created from a string.
 *
 * @param Pat           The string pattern, transformed into a pattern instance
 *                          and then is pushed into the set of patterns; r-value.
 */
void patternset::push_back(std::string & Pat){
    push_back(pattern(Pat));
}
/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 * A pattern instance is created from a string.
 *
 * @param Pat           The string pattern, transformed into a pattern instance
 *                          and then is pushed into the set of patterns; l-value.
 */
void patternset::push_back(std::string && Pat){
    push_back(pattern(Pat));
}
/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 * A pattern instance is created from character vector.
 *
 * @param Pat           The vector pattern, transformed into a pattern instance
 *                          and then is pushed into the set of patterns; r-value.
 */
void patternset::push_back(std::vector<char> & Pat){
    push_back(pattern(Pat));
}
/**
 * Pushes a pattern into the set of pattern and adjusts the
 * mini-/maximal values for weight and DC positions.
 * A pattern instance is created from character vector.
 *
 * @param Pat           The vector pattern, transformed into a pattern instance
 *                          and then is pushed into the set of patterns; r-value.
 */
void patternset::push_back(std::vector<unsigned char> & Pat){
    push_back(pattern(Pat));
}

/**
 * Creates a random set with given, passed parameter.
 * Here all patterns have the same weight and length.
 *
 * @param Size          The number of patterns.
 *
 * @param Weight        The weight of each pattern.
 *
 * @param DontCare      The number of don't care positions in each pattern.
 *
 * @param Uniq          Boolean; true, if each pattern has to be a unique
 *                          pattern.
 */
void patternset::random(unsigned Size, unsigned Weight, unsigned DontCare, bool Uniq){
    random(Size, Weight, Weight, DontCare, DontCare, Uniq);
}
/**
 * Creates a random set with given, passed parameter.
 * Here all patterns have the same number of DC positions.
 *
 * @param Size          The number of patterns.
 *
 * @param MinWeight     The minimal weight for a pattern.
 *
 * @param MaxWeight     The maximal weight for a pattern.
 *
 * @param DontCare      The number of don't care positions in each pattern.
 *
 * @param Uniq          Boolean; true, if each pattern has to be a unique
 *                          pattern.
 */
void patternset::random(unsigned Size, unsigned MinWeight, unsigned MaxWeight,
    unsigned DontCare, bool Uniq){
    random(Size, MinWeight, MaxWeight, DontCare, DontCare, Uniq);
}
/**
 * Creates a random set with given, passed parameter.
 *
 * @param Size          The number of patterns.
 *
 * @param MinWeight     The minimal weight for a pattern.
 *
 * @param MaxWeight     The maximal weight for a pattern.
 *
 * @param MinDontCare   The minimal number of don't care positions for patterns.
 *
 * @param MaxDontCare   The maximal number of don't care positions for patterns.
 *
 * @param Uniq          Boolean; true, if each pattern has to be a unique
 *                          pattern.
 */
void patternset::random(unsigned Size, unsigned MinWeight, unsigned MaxWeight,
    unsigned MinDontCare, unsigned MaxDontCare, bool Uniq){
    _PatternSet.clear();
    if(Size == 0){
        std::cerr << "Illegal value for pattern number!\nSize has to be greater than 0!" << std::endl;
        std::exit(-1);
    }
    _adjust(MinWeight,MaxWeight);
    _adjust(MinDontCare,MaxDontCare);
    if(Size >= 2){
        _PatternSet.push_back(pattern(MinWeight,MinDontCare));
        double StepW = (MaxWeight - MinWeight + 1)/(double)(Size);
        double StepD = (MaxDontCare - MinDontCare + 1)/(double)(Size);
        double CurrentW = MinWeight + StepW;
        double CurrentD = MinDontCare + StepD;
        unsigned PatSize = 0;
        char TryCtr = 0;
        while(PatSize < Size-2){
            pattern Pat((unsigned)CurrentW, (unsigned)CurrentD,(uint64_t)PatSize+1);
            TryCtr = 0;
            while(Uniq == false && is_uniq(Pat) == false && TryCtr++ < 100){
                Pat.random((unsigned)CurrentW, (unsigned)CurrentD);
            }
            CurrentW += StepW;
            CurrentD += StepD;
            if(TryCtr < 100){
                _PatternSet.push_back(Pat);
                PatSize++;
            }
            else if(TryCtr >= 100 && (unsigned)CurrentW >= MaxWeight && (unsigned)CurrentD >= MaxDontCare){
                break;
            }
        }
        TryCtr = 0;
        pattern Pat = pattern(MaxWeight,MaxDontCare,_PatternSet.size());
        while(Uniq == false && is_uniq(Pat) == false && TryCtr++ < 100){
            Pat.random(MaxWeight,MaxDontCare);
        }
        if(TryCtr < 100){
            _PatternSet.push_back(Pat);
        }
    }
    else{
        _PatternSet.push_back(pattern((MinWeight + MaxWeight)/2,(MinDontCare + MaxDontCare)/2));
    }
}

/**
 * Sorts the set of pattern according to the score saved in each pattern.
 */
void patternset::sort(){
    std::sort(_PatternSet.begin(), _PatternSet.end());
}

/**
 * Checks, if a pattern is pairwise unique with all patterns occuring
 * in the set. If the pattern is part of the set, the identical comparison
 * is ignored.
 *
 * @param Pat           The pattern that is compared to the set
 *
 * @return              Returns if this pattern is unique to the set.
 */
bool patternset::is_uniq(const pattern & Pat) const{
    return (std::find_if(_PatternSet.begin(),_PatternSet.end(), [& Pat](const pattern & P){
        if( & Pat != & P){
            return (Pat == P);
        }
        return false;
    }) == _PatternSet.end());
}

/**
 * Performs a random swap on a pattern an checks, if the resulting pattern
 * is unique in the set; accepts swap if unique.
 *
 * @param Idx           The index of the pattern in the set that should be
 *                          swapped randomly.
 */
void patternset::random_swap_uniq(unsigned Idx){
    pattern PatSave = _PatternSet[Idx];
    for(unsigned i = 0; i < (size()+10)*(size()+10); i++){
        _PatternSet[Idx].random_swap();
        if(is_uniq(_PatternSet[Idx]) == true){
            return;
        }
        _PatternSet[Idx] = PatSave;
    }
}

/**
 * Stores the pattern set into a file.
 *
 * @param OutFile       The name of the file, in which the patterns should be
 *                      written into.
 */
void patternset::to_file(std::string OutFile){
    std::ofstream Output(OutFile);
    for(auto & Pat : _PatternSet){
        Output << Pat.to_string() << std::endl;
    }
    Output.close();
}


/**
 * Random access operator, that returns the pattern from a specific
 * position in the set.
 *
 * @param Idx           The index number of a pattern, that sould be returned.
 *
 * @return              The pattern from the Idx-th position; constant
 */
pattern patternset::operator[](size_t Idx) const{
    return _PatternSet[Idx];
}
/**
 * Random access operator, that returns the pattern from a specific
 * position in the set.
 *
 * @param Idx           The index number of a pattern, that sould be returned.
 *
 * @return              The pattern from the Idx-th position; r-value
 */
pattern & patternset::operator[](size_t Idx){
    return _PatternSet[Idx];
}

/**
 * Compares the score of the pattern set with another set and returns true,
 * if this pattern set has the lower score.
 *
 * @param P             The patternset that is used for score comparison
 *
 * @return              True, if this pattern set has a lower score
 */
bool patternset::operator<(const patternset & P) const{
    return _Score < P._Score;
}
/**
 * Compares the score of the pattern set with another set and returns true,
 * if this pattern set has the higher score.
 *
 * @param P             The patternset that is used for score comparison
 *
 * @return              True, if this pattern set has a higher score
 */
bool patternset::operator>(const patternset & P) const{
    return _Score > P._Score;
}


/**
 * Returns the maximum weight in this set.
 *
 * @return              Integer value, the maximum weight
 */
unsigned patternset::max_weight() const{
    return _MaxWeight;
}
/**
 * Returns the minimum weight in this set.
 *
 * @return              Integer value, the minimum weight
 */
unsigned patternset::min_weight() const{
    return _MinWeight;
}
/**
 * Returns the average weight in this set.
 *
 * @return              Integer value, the average weight
 */
unsigned patternset::weight() const{
    unsigned MeanWeight = 0;
    for(auto & Pat : _PatternSet){
        MeanWeight += Pat.weight();
    }
    return MeanWeight/_PatternSet.size();
}

/**
 * Returns the maximum number of don't care positions in this set.
 *
 * @return              Integer value, the maximum DC position number.
 */
unsigned patternset::max_dontcare() const{
    return _MaxDontCare;
}
/**
 * Returns the minimum number of don't care positions in this set.
 *
 * @return              Integer value, the minimum DC position number.
 */
unsigned patternset::min_dontcare() const{
    return _MinDontCare;
}
/**
 * Returns the average number of don't care positions in this set.
 *
 * @return              Integer value, the average DC position number.
 */
unsigned patternset::dontcare() const{
    unsigned MeanDC = 0;
    for(auto & Pat : _PatternSet){
        MeanDC += Pat.dontcare();
    }
    return MeanDC/_PatternSet.size();
}

/**
 * Returns the maximum pattern length in this set.
 *
 * @return              Integer value, the maximum pattern length
 */
unsigned patternset::max_length() const{
    return _MaxWeight + _MaxDontCare;
}
/**
 * Returns the minimum pattern length in this set.
 *
 * @return              Integer value, the minimum pattern length
 */
unsigned patternset::min_length() const{
    return _MinWeight + _MinDontCare;
}
/**
 * Returns the average pattern length in this set.
 *
 * @return              Integer value, the average pattern length
 */
unsigned patternset::length() const{
    unsigned MeanLength = 0;
    for(auto & Pat : _PatternSet){
        MeanLength += Pat.length();
    }
    return MeanLength/_PatternSet.size();
}

/**
 * The number of patterns in this set.
 *
 * @return              The number of patterns.
 */
unsigned patternset::size() const{
    return _PatternSet.size();
}

/**
 * The score saved to this set (has to be set before).
 *
 * @return              The pattern set score
 */
double patternset::score() const{
    return _Score;
}

/**
 * Stores a score for this pattern set. Can be used to compared different set of
 * patterns.
 *
 * @param Score         Double value, saved to this pattern set instance.
 */
void patternset::set_score(double Score){
    _Score = Score;
}


/**
 * Begin iterator that can be used to iterate over the set of patterns.
 *
 * @return              The begin iterator of the pattern set vector.
 */
patternset::iterator patternset::begin(){
    return _PatternSet.begin();
}
/**
 * End iterator that can be used to iterate over the set of patterns.
 *
 * @return              The end iterator of the pattern set vector.
 */
patternset::iterator patternset::end(){
    return _PatternSet.end();
}


/**
 * Function that is used to parse patterns from a file/string.
 * 
 * @param S             String, possibly containing multiple patterns
 *
 * @return              A vector of pattern strings.
 */
std::vector<std::string> patternset::_split(std::string &S){
    std::vector<std::string> Splits;
    std::string::iterator SPos = S.begin();
    for(auto Iter = S.begin(); Iter != S.end(); Iter++){
        switch(*Iter){
            case ';':
            case ',':
            case '.':
            case ' ':
            case '\t':
                std::string Tmp(SPos,Iter);
                if(Tmp.size() != 0){
                    Splits.push_back(Tmp);
                }
                SPos = Iter+1;
                break;
        }
    }
    std::string Tmp(SPos,S.end());
    if(Tmp.size() != 0){
        Splits.push_back(Tmp);
    }
    return Splits;
}