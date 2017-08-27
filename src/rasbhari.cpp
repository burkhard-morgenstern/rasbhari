/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari object file
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
#include "rasbhari.hpp"

/**
 * The empty-constructor, can be used, if patterns should be pushed into the
 * rasbhari instance for optimising.
 */
rasbhari::rasbhari():_RasbhariScore(0),_RasbhariSensitivity(-1), _PatNo(0){
}
/**
 * The copy-constructor, can be used, if an already existing rasbhari instance
 * shall be copied.
 *
 * @param RasbObj       The rasbhari instance that should be copied.
 */
rasbhari::rasbhari(const rasbhari &RasbObj){
    _CoefMat = RasbObj._CoefMat;
    _PatternList = RasbObj._PatternList;
    _RasbhariPattern = RasbObj._RasbhariPattern;
    _RasbhariScore = RasbObj._RasbhariScore;
    _RasbhariSensitivity = RasbObj._RasbhariSensitivity;
    _PatNo = RasbObj._PatNo;
}
/**
 * The std-constructor-1, can be used, if patterns should be created randomly.
 * Each pattern has the same weight and dc-position number.
 *
 * @param Size          The number of patterns.
 *
 * @param Weight        The number of match positions, i.e. the weight.
 *
 * @param DontCare      The number of don't-care positions
 */
rasbhari::rasbhari(unsigned Size, unsigned Weight, unsigned DontCare):_RasbhariScore(0),_RasbhariSensitivity(-1), _PatNo(0){
    _check_pattern_number(Size, Weight, DontCare, DontCare);
    _RasbhariPattern = patternset(Size,Weight,DontCare,true);
    _make_pattern_list();
    _adjust_coef_mat();
    calculate();
    calculate_sensitivity();
}
/**
 * The std-constructor-2, can be used, if patterns should be created randomly.
 * Each pattern has the same weight and dc-position number.
 *
 * @param Size          The number of patterns.
 *
 * @param Weight        The number of match positions, i.e. the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions
 *
 * @param MaxDontCare   The maximal number of don't-care positions
 */
rasbhari::rasbhari(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare):_RasbhariScore(0),_RasbhariSensitivity(-1), _PatNo(0){\
    _check_pattern_number(Size, Weight, MinDontCare, MaxDontCare);
    _RasbhariPattern = patternset(Size,Weight,Weight,MinDontCare,MaxDontCare,true);
    _make_pattern_list();
    _adjust_coef_mat();
    calculate();
    calculate_sensitivity();
}
/**
 * The file-constructor, can be used, if the patterns should be read from a
 * given file.
 *
 * @param PatternFile   String, containing the name for a file containing a
 *                          pattern set.
 */
rasbhari::rasbhari(std::string PatternFile):_RasbhariScore(0),_RasbhariSensitivity(-1), _PatNo(0){\
    _RasbhariPattern = patternset(PatternFile);
    _make_pattern_list();
    _adjust_coef_mat();
    calculate();
    calculate_sensitivity();
}


/**
 * The push-bach functions allows to dynamically add patterns to a set, if it
 * is wished. E.g. the improvement with a higher number of pattern shall be
 * investigated.
 * 
 * @param Pat           The pattern to be pushed to the patternset; l-value
 */
void rasbhari::push_back(pattern && Pat){
    push_back(Pat);
}
/**
 * The push-bach functions allows to dynamically add patterns to a set, if it
 * is wished. E.g. the improvement with a higher number of pattern shall be
 * investigated.
 * 
 * @param Pat           The pattern to be pushed to the patternset; r-value
 */
void rasbhari::push_back(pattern & Pat){
    if(Pat.weight() != _RasbhariPattern.weight()){
        return;
    }
    _PatternList.push_back(Pat);
    _RasbhariPattern.push_back(Pat);
    _RasbhariPattern[_RasbhariPattern.size()-1].set_idx(_RasbhariPattern.size()-1);
    _PatternList[_PatternList.size()-1].set_idx(_PatternList.size()-1);
    _adjust_coef_mat();
    _RasbhariScore = -1;
    _RasbhariSensitivity = -1;
    rasb_opt::Size++;
    rasb_opt::MinDontcare = std::min(Pat.dontcare(),_RasbhariPattern.min_dontcare());
    rasb_opt::MaxDontcare = std::max(Pat.dontcare(),_RasbhariPattern.max_dontcare());
}

/**
 * Since often needed a function that resets the coefficient matrix for oc/var.
 */
void rasbhari::_adjust_coef_mat(){
    _CoefMat = std::vector< std::vector<double> >(_RasbhariPattern.size(),std::vector<double>(_RasbhariPattern.size(),0));
}
/**
 * Since the order of the the rasbhari pattern changes due to the changing
 * contribution by random permutations, another list is needed, holding
 * the correct order for the coefficient matrix.
 */
void rasbhari::_make_pattern_list(){
    _PatternList = std::vector<pattern>(_RasbhariPattern.size());
    for(unsigned i = 0; i < _RasbhariPattern.size(); i++){
        _RasbhariPattern[i].set_idx(i);
        _PatternList[i] = _RasbhariPattern[i];
    }
}
/**
 * For the entire set, the total coef_matrix is set up and pattern contribution
 * were calculated. This part is only the initial step. Changes on patterns
 * are only for single patterns, thus only one column and row changes!
 */
void rasbhari::calculate(){
    for(unsigned i = 0; i < _PatternList.size(); i++){
        double CoefPat = 0;
        for(unsigned j = 0; j < _PatternList.size(); j++){
            _CoefMat[i][j] = calculate_pair(i,j);
            CoefPat += _CoefMat[i][j];
        }
        _PatternList[i].set_score(CoefPat);
    }
    _RasbhariScore = 0;
    for(unsigned i = 0; i < _PatternList.size(); i++){
        for(unsigned j = i; j < _PatternList.size(); j++){
            _RasbhariScore += _CoefMat[i][j];
        }
    }
    for(auto & Pat : _RasbhariPattern){
        Pat.set_score(_PatternList[Pat.idx()].score());
    }
}

/**
 * Calculates for a specific pair either the overlap complexity or the variance.
 * 
 * @param Idx1          The index of the first pattern of a pattern pair.
 *
 * @param Idx2          The index of the second pattern of a pattern pair.
 *
 * @return              The OC/variance of the pattern pair.
 */
double rasbhari::calculate_pair(size_t Idx1, size_t Idx2){
    if(rasb_opt::Oc){
        return rasbhari_compute::pair_coef_oc(_PatternList[Idx1],_PatternList[Idx2]);
    }
    if(Idx1 == Idx2){
        return rasbhari_compute::pair_coef_var_sym(_PatternList[Idx1],_PatternList[Idx2]);
    }
    return rasbhari_compute::pair_coef_var(_PatternList[Idx1],_PatternList[Idx2]);    
}
/**
 * If whished, calculates for the entire patternset the sensitivity.
 */
void rasbhari::calculate_sensitivity(){
    if(rasb_opt::Sens){
        _RasbhariSensitivity = rasbhari_compute::sensitivity(_RasbhariPattern);
    }
}

/**
 * Changes on one pattern will change the coef_matrix in only one row and one
 * column, therefore the OC/variance coefficient only needs to be updated!
 *
 * @param Idx           The index of the pattern that was permutated randomly.
 */
void rasbhari::update(unsigned Idx){
    double CoefPat = 0;
    for(unsigned i = 0; i < _PatternList.size(); i++){
        _RasbhariScore -= _CoefMat[i][Idx];
        double PatternScore = _PatternList[i].score() - _CoefMat[i][Idx];
        _CoefMat[i][Idx] = calculate_pair(i,Idx);
        _CoefMat[Idx][i] = _CoefMat[i][Idx];
        CoefPat += _CoefMat[i][Idx];
        _PatternList[i].set_score(PatternScore+_CoefMat[i][Idx]);
        _RasbhariScore += _CoefMat[i][Idx];
    }
    _PatternList[Idx].set_score(CoefPat);
    for(auto & Pat : _RasbhariPattern){
        Pat.set_score(_PatternList[Pat.idx()].score());
    }
}
/**
 * The actual optimising step. The pattern with highest contribute is 
 * investigated, a random permutation is performed and afterwards it is checked,
 * if the OC/variance is optimised. If not, the change is undone, otherwise
 * it is accepted.
 *
 * @return              Returns if the permutation was succesfull or not.
 */
bool rasbhari::climb_hill(){
    unsigned PatIdx = _PatNo%_RasbhariPattern.size();
    unsigned OrigIdx = _RasbhariPattern[PatIdx].idx();

    double LastCoef = _RasbhariScore;
    pattern LastPattern = _RasbhariPattern[PatIdx];
    std::vector< std::vector<double> > LastCMat = _CoefMat;
    std::vector<pattern> LastPatList = _PatternList;

    _RasbhariPattern.random_swap_uniq(PatIdx);
    _PatternList[OrigIdx] = _RasbhariPattern[PatIdx];
    update(OrigIdx);
    
    if(_RasbhariScore < LastCoef){
        _PatNo = 0;
        _RasbhariPattern.sort();
        return true;
    }

    std::swap(_RasbhariPattern[PatIdx], LastPattern);
    std::swap(_CoefMat,LastCMat);
    std::swap(_PatternList,LastPatList);
    for(auto & Pat : _RasbhariPattern){
        Pat.set_score(_PatternList[Pat.idx()].score());
    }
    _RasbhariScore = LastCoef;
    _PatNo++;
    return false;
}
/**
 * The actuall hillclimbing process. For a specific number, Limit,
 * the optimising step is done to optimise the set.
 *
 * @param Limit         The number of permutation tries executed on the pattern
 *                          set.
 *
 * @return              Returns if an coefficient improvement took place or not.
 */
bool rasbhari::hill_climbing(unsigned Limit){
    double ScoreBest = _RasbhariScore, InitialScore = _RasbhariScore;
    _PatNo = 0;
    unsigned ModeSave = rasb_opt::ImproveMode, Ctr = 0;
    if(rasb_opt::ImproveMode < 1 || (rasb_opt::OptOc <= 1 && rasb_opt::Sens == false)){
        rasb_opt::ImproveMode = 1;
    }
    for(unsigned i = 0; i < Limit; i++){
        if(!rasb_opt::Silent && rasb_opt::ImproveMode == 1){
            std::cout << "\rStep " << i << "/" << Limit << "  Improvement +" << Ctr;
            std::cout.flush();
        }
        climb_hill();
        if(ScoreBest > _RasbhariScore && !rasb_opt::Silent && rasb_opt::ImproveMode == 1){
            ScoreBest = _RasbhariScore;
            Ctr++;
            if(!rasb_opt::Silent && !rasb_opt::Quiet && rasb_opt::ImproveMode == 1){
                std::cout << std::endl;
                print();
            }
        }
    }
    if(!rasb_opt::Silent && rasb_opt::ImproveMode == 1){
        std::cout << "\r"<< std::string(80,' ') << "\r";
        std::cout << "Number of improvements: " << Ctr << std::endl << std::endl;
    }
    rasb_opt::ImproveMode = ModeSave;
    return InitialScore > _RasbhariScore;
}   
/**
 * The iterative hillclimbing process. For a specific number, Limit,
 * the optimising step is done to optimise each set of Iteration-many pattern
 * sets.
 * Thus, in each of the Iteration-many steps, a random patternset is generated,
 * the best of all is taken.
 *
 * @param Limit         The number of permutation tries executed on the pattern
 *                          set.
 *
 * @param Iteration     The number of random intial pattern sets.
 *
 * @return              Returns if an coefficient improvement took place or not.
 */
bool rasbhari::iterate_hill_climbing(unsigned Limit, unsigned Iteration){
    patternset HillClimbBest = _RasbhariPattern;
    std::vector<pattern> PatListBest = _PatternList;
    std::vector< std::vector<double> > CoefBest = _CoefMat;
    double ScoreBest = _RasbhariScore, InitialScore = ScoreBest;
    unsigned ModeSave = rasb_opt::ImproveMode, Ctr = 0;
    if(rasb_opt::ImproveMode < 2 || rasb_opt::Sens == false){
        rasb_opt::ImproveMode = 2;
    }
    if(rasb_opt::Sens == true && Iteration == 0){
        Iteration = 1;
    }
    for(unsigned i = 0; i < Iteration; i++){
        if(!rasb_opt::Silent && rasb_opt::ImproveMode == 2){
            std::cout << "\rStep " << i << "/" << Iteration << "  Improvement +" << Ctr;
            std::cout.flush();
        }
        _RasbhariPattern = patternset(rasb_opt::Size,rasb_opt::Weight,rasb_opt::Weight,rasb_opt::MinDontcare,rasb_opt::MaxDontcare,true);
        _make_pattern_list();
        _adjust_coef_mat();
        calculate();
        hill_climbing(Limit);
        if(ScoreBest > _RasbhariScore){
            Ctr++;
            if(!rasb_opt::Silent && !rasb_opt::Quiet && rasb_opt::ImproveMode == 2){
                std::cout << std::endl;
                print();
            }
            ScoreBest = _RasbhariScore;
            std::swap(HillClimbBest,_RasbhariPattern);
            std::swap(PatListBest,_PatternList);
            std::swap(CoefBest,_CoefMat);
        }
    }
    if(!rasb_opt::Silent && rasb_opt::ImproveMode == 2){
        std::cout << "\r"<<  std::string(80,' ') << "\r";
        std::cout << "Number of improvements: " << Ctr << std::endl << std::endl;
    }
    rasb_opt::ImproveMode = ModeSave;
    _RasbhariScore = ScoreBest;
    std::swap(_RasbhariPattern,HillClimbBest);
    std::swap(_PatternList,PatListBest);
    std::swap(_CoefMat,CoefBest);
    return InitialScore > _RasbhariScore;
}
/**
 * The iterative hillclimbing process with sensitivity calculation. For a
 * specific number, Limit, the optimising step is done to optimise each set of
 * Iteration-many pattern sets.
 * Thus, in each of the Iteration-many steps, a random patternset is generated,
 * the best of all is taken.
 * For the final set, the sensitivity is calculated and the set will be
 * accepted, if the sensitivity is improved.
 *
 * @param Limit         The number of permutation tries executed on the pattern
 *                          set.
 *
 * @param Iteration     The number of random intial pattern sets.
 *
 * @param InitSens      True, if the sensitivity shall be calculated
 *                          (eventually not for the initial random set!)
 *
 * @return              Returns if an coefficient improvement took place or not.
 */
bool rasbhari::climb_hill_sensitivity(unsigned Limit, unsigned Iteration, bool InitSens){
    if(InitSens){
        calculate_sensitivity();
    }
    double SensBest = _RasbhariSensitivity;
    iterate_hill_climbing(Limit,Iteration);
    calculate_sensitivity();
    return SensBest < _RasbhariSensitivity;
}
/**
 * The iterative hillclimbing process with sensitivity optimisation. For a
 * specific number, Limit, the optimising step is done to optimise each set of
 * Iteration-many pattern sets.
 * Thus, in each of the Iteration-many steps, a random patternset is generated,
 * the best of all is taken.
 * These steps are don Loop-many times, the sensitivity of all sets is
 * calculated and the set with highest sensitivity will be returned.
 *
 * @param Limit         The number of permutation tries executed on the pattern
 *                          set.
 *
 * @param Iteration     The number of random intial pattern sets.
 *
 * @param Loop          The number of pattern sets for which the sensitivity
 *                          should be calculated. Each pattern set is the best
 *                          of Iteration-many optimised pattern sets.
 *
 * @return              Returns if an coefficient improvement took place or not.
 */
bool rasbhari::hill_climbing_sensitivity(unsigned Limit, unsigned Iteration, unsigned Loop){
    patternset HillClimbBest = _RasbhariPattern;
    std::vector<pattern> PatListBest = _PatternList;
    std::vector< std::vector<double> > CoefBest = _CoefMat;
    double SensBest = _RasbhariSensitivity, InitialSens = _RasbhariSensitivity;
    unsigned Ctr = 0;
    unsigned ModeSave = rasb_opt::ImproveMode;
    if(rasb_opt::ImproveMode < 3 && rasb_opt::Sens){
        rasb_opt::ImproveMode = 3;
    }
    for(unsigned i = 0; i < Loop; i++){
        if(!rasb_opt::Silent && rasb_opt::ImproveMode == 3){
            std::cout << "\rStep " << i << "/" << Loop << "  Improvement +" << Ctr;
            std::cout.flush();
        }
        _RasbhariPattern = patternset(rasb_opt::Size,rasb_opt::Weight,rasb_opt::Weight,rasb_opt::MinDontcare,rasb_opt::MaxDontcare,true);
        _make_pattern_list();
        _adjust_coef_mat();
        calculate();
        climb_hill_sensitivity(Limit, Iteration);
        if(SensBest < _RasbhariSensitivity){
            Ctr++;
            if(!rasb_opt::Silent && !rasb_opt::Quiet && rasb_opt::ImproveMode == 3){
                std::cout << std::endl;
                print();
            }
            SensBest = _RasbhariSensitivity;
            std::swap(HillClimbBest,_RasbhariPattern);
            std::swap(PatListBest,_PatternList);
            std::swap(CoefBest,_CoefMat);
        }
    }
    if(!rasb_opt::Silent && rasb_opt::ImproveMode == 3){
        std::cout << "\r" << std::string(80,' ') << "\r";
        std::cout << "Number of improvements: " << Ctr << std::endl << std::endl;
    }
    rasb_opt::ImproveMode = ModeSave;
    _RasbhariSensitivity = SensBest;
    std::swap(_RasbhariPattern,HillClimbBest);
    std::swap(_PatternList,PatListBest);
    std::swap(_CoefMat,CoefBest);
    return InitialSens < _RasbhariSensitivity;
}


/**
 * Returns the pattern set optimised by rasbhari conditions.
 *
 * @return              The optimised pattern set
 */
patternset rasbhari::pattern_set() const{
    return _RasbhariPattern;
}
/**
 * Returns the pattern set (Reference) optimised by rasbhari conditions.
 *
 * @return              The optimised pattern set
 */
patternset & rasbhari::pattern_set(){
    return _RasbhariPattern;
}

/**
 * Returns the pattern set score, either OC or variance.
 *
 * @return              Double value, the rasbhari score (OC/var)
 */
double rasbhari::score() const{
    return _RasbhariScore;
}
/**
 * Returns the pattern set sensitivity.
 *
 * @return              Double value, the rasbhari sensitivity
 */
double rasbhari::sensitivity() const{
    return _RasbhariSensitivity;
}

/**
 * Returns the pattern set size.
 *
 * @return              Pattern set size.
 */
unsigned rasbhari::size() const{
    return _RasbhariPattern.size();
}


/**
 * Prints all relevant information about the pattern set to std::cout.
 */
void rasbhari::print(){
    for(auto Pat : _RasbhariPattern){
        std::cout << Pat.idx() << " " << Pat.to_string() << std::endl;
    }
    std::cout << "rasbhari coefficient        : " << _RasbhariScore << std::endl;
    std::cout << "rasbhari coefficient (norm) : " << _RasbhariScore / ((_RasbhariPattern.size()*(_RasbhariPattern.size()+1))/(double)2) << std::endl;
    if(rasb_opt::Sens){
        std::cout << "rasbhari sensitivity        : " << _RasbhariSensitivity << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
}
/**
 * Prints all relevant information about the pattern set into a file.
 *
 * @param OutFile       String, containing the output file-name for pattern set.
 */
void rasbhari::to_file(std::string OutFile){
    if(OutFile.size() == 0){
        return;
    }
    _RasbhariPattern.to_file(OutFile);
    std::ofstream Output(OutFile, std::ios::app);
    Output << "#rasbhari coefficient        : " << _RasbhariScore << std::endl;
    Output << "#rasbhari coefficient (norm) : " << _RasbhariScore / ((_RasbhariPattern.size()*(_RasbhariPattern.size()+1))/(double)2) << std::endl;
    if(rasb_opt::Sens){
        Output << "#rasbhari sensitivity        : " << _RasbhariSensitivity << std::endl;
    }
    Output.close();
    if(!rasb_opt::Silent){
        std::cout << "rasbhari set written to file:\n    '" << OutFile << "'" << std::endl;
    }
}


/**
 * Random access operator, that returns the pattern from a specific
 * position in the rasbhari set.
 *
 * @param Idx           The index number of a pattern, that sould be returned.
 *
 * @return              The pattern from the Idx-th position; constant
 */
pattern rasbhari::operator[](size_t Idx) const{
    return _PatternList[Idx];
}
/**
 * Random access operator, that returns the pattern from a specific
 * position in the rasbhari set.
 *
 * @param Idx           The index number of a pattern, that sould be returned.
 *
 * @return              The pattern from the Idx-th position; r-value
 */
pattern & rasbhari::operator[](size_t Idx){
    return _PatternList[Idx];
}

/**
 * An begin iterator is returned, that can be used to iterate over the pattern
 * set.
 *
 * @return              Begin-Iterator for the pattern set.
 */
rasbhari::iterator rasbhari::begin(){
    return _PatternList.begin();
}
/**
 * An end iterator is returned, that can be used to iterate over the pattern
 * set.
 *
 * @return              End-Iterator for the pattern set.
 */
rasbhari::iterator rasbhari::end(){
    return _PatternList.end();
}


/**
 * Some Debug information, this is not for your concern.
 */
void rasbhari::_debug(){
    std::cout << "RasbPat Score Idx\tListPat Score Idx\t\tpairwise idx identic?" << std::endl;
    for(unsigned i = 0; i < _RasbhariPattern.size(); i++){
        std::cout << _RasbhariPattern[i].to_string() << " " << _RasbhariPattern[i].score() << " " << _RasbhariPattern[i].idx() << " \t";
        std::cout << _PatternList[i].to_string() << " " << _PatternList[i].score() << " " << _PatternList[i].idx() << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Coef Mat " << _RasbhariPattern.size() << "x" << _RasbhariPattern.size() << "\tSymmetric?" << std::endl;
    for(unsigned i = 0; i < _PatternList.size(); i++){
        for(unsigned j = 0; j < _PatternList.size(); j++){
            std::cout << _CoefMat[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "_RasbhariScore: " << _RasbhariScore/((rasb_opt::Size*(rasb_opt::Size+1))/2) << std::endl;
    std::cout << "_RasbhariSensitivity: " << _RasbhariSensitivity << std::endl << std::endl << std::endl;
}

/**
 * Calculates for minimal and maximal DC positions with its weight and pattern
 * set number the maximal possible number of unique patterns. If the number is
 * to high, the values will be adjusted, a warning is thrown.
 */
void rasbhari::_check_pattern_number(unsigned & Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    double MaxPatNo = _max_pat_no(Weight, MinDontCare, MaxDontCare);
    if((double) Size > MaxPatNo){
        printf("%c[1;33m", 27);
        std::cerr << "The chosen number of " << Size << " unique patterns is to high for your set configuration!" << std::endl;
        std::cerr << "The number of patterns will be adjusted to " << MaxPatNo << " patterns!" << std::endl << std::endl;
        printf("%c[0m", 27);
        Size = (unsigned)MaxPatNo;
        rasb_opt::Size = (unsigned)MaxPatNo;
    }
}
/**
 * Calculates the maximal number for a pattern set condition; uses the binomial
 * coeficient.
 *
 * @param Weight        The number of match positions in the pattern set
 *
 * @param MinDontCare   The minimal number of don't care positions
 *
 * @param MinDontCare   The maximal number of don't care positions
 */
double rasbhari::_max_pat_no(unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    if(MinDontCare > MaxDontCare){
        std::swap(MinDontCare,MaxDontCare);
    }
    double MaxPatCount = 0;
    for(unsigned d = MinDontCare; d <= MaxDontCare; d++){
        MaxPatCount += _binom_coef(d+Weight-2,Weight-2);
    }
    return MaxPatCount;
}
/**
 * To overcome problems with overflow for double and the binomial coefficient, 
 * an alternate form for it is used, which is calculated more efficient.
 * Formula: 'n over k' or 'n binom k'
 *
 * @param n             The size of a set with entities to be chosen.
 *
 * @param k             The number entities that should be chosen from a set
 *                          with n elements.
 */
double rasbhari::_binom_coef(unsigned n, unsigned k){
    if(n < k){
        std::swap(n,k);
    }
    double Coeff = 1;
    for(unsigned i = 1; i <= k; i++){
        Coeff *= (n-k+i)/((double)i);
    }
    return Coeff;
}