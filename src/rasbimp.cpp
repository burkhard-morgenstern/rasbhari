/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight. It is possible to improve your patternset and
 * read patterns from a file.
 *
 * rasbhari implementation namespace file
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
#include "rasbimp.hpp"

/**
 * Creates for patternset-parameters Size, Weight, Min-/MaxDontCare a randomly
 * initialised patternset, optimises it with rasb_opt::Limit many of random
 * permutations on patterns contributing to the variance at most.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_var(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = false;
    return _hillclimb(Size, Weight, MinDontCare, MaxDontCare);
}
/**
 * Creates for patternset-parameters Size, Weight, Min-/MaxDontCare a randomly
 * initialised patternset, optimises it with rasb_opt::Limit many of random
 * permutations on patterns contributing to the overlap-complexity at most.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_oc(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = true;
    return _hillclimb(Size, Weight, MinDontCare, MaxDontCare);
}

/**
 * Creates for patternset-parameters Size, Weight, Min-/MaxDontCare 
 * rasb_opt::OptOC randomly initialised patternset, optimises it with
 * rasb_opt::Limit many of random permutations on patterns contributing to the
 * overlap-complexity at most.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_oc_iterative(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = true;
    return _hillclimb_iterative(Size,Weight,MinDontCare,MaxDontCare);
}
/**
 * Creates for patternset-parameters Size, Weight, Min-/MaxDontCare 
 * rasb_opt::OptOc many randomly initialised patternset, optimises it with
 * rasb_opt::Limit many of random permutations on patterns contributing to the
 * variance at most.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_var_iterative(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = false;
    return _hillclimb_iterative(Size,Weight,MinDontCare,MaxDontCare);
}

/**
 * Chooses for patternset-parameters Size, Weight, Min-/MaxDontCare
 * rasb_opt::OptSens many patternsets, that are the best of rasb_opt::OptOc
 * many randomly initialised patternsets, which were OC optimised with
 * rasb_opt::Limit many of random permutations, the patternset with the highest
 * sensitivity.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_sens_oc(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = true;
    return _hillclimb_sens(Size,Weight,MinDontCare,MaxDontCare);
}
/**
 * Chooses for patternset-parameters Size, Weight, Min-/MaxDontCare
 * rasb_opt::OptSens many patternsets, that are the best of rasb_opt::OptOc
 * many randomly initialised patternsets, which were variance optimised with
 * rasb_opt::Limit many of random permutations, the patternset with the highest
 * sensitivity.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 */
rasbhari rasb_implement::hillclimb_sens_var(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::Oc = false;
    return _hillclimb_sens(Size,Weight,MinDontCare,MaxDontCare);
}

/**
 * Wrapper for the simple hill-climbing process for the OC/Var; depending 
 * on boolean output variables, some information will be printed to std::cout.
 * The final rasbhari set is returned.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 *
 * @return              Returns the final rasbhari patternset output.
 */
rasbhari rasb_implement::_hillclimb(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::ImproveMode = 1;
    rasb_opt::Sens = false;
    rasbhari RasbSet;
    if(rasb_opt::InFile.size() != 0){
        RasbSet = rasbhari(rasb_opt::InFile);
    }
    else{
        RasbSet = rasbhari(Size,Weight,MinDontCare,MaxDontCare);
    }
    rasb_opt::Size = RasbSet.pattern_set().size();
    rasb_opt::Weight = RasbSet.pattern_set().max_weight();
    rasb_opt::MinDontcare = RasbSet.pattern_set().min_dontcare();
    rasb_opt::MaxDontcare = RasbSet.pattern_set().max_dontcare();
    unsigned FillSizeBegin = 0, FillSizeEnd = 0;
    if(MaxDontCare + Weight > 13){
        FillSizeBegin = (MaxDontCare+Weight-13)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-13)%2 != 0){
            FillSizeEnd++;
        }
    }
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Initial Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    RasbSet.hill_climbing();
    if(MaxDontCare + Weight > 15){
        FillSizeBegin = (MaxDontCare+Weight-15)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-15)%2 != 0){
            FillSizeEnd++;
        }
    }
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Optimised Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    return RasbSet;
}
/**
 * Wrapper for the iterative hill-climbing process for the OC/Var; depending 
 * on boolean output variables, some information will be printed to std::cout.
 * The final rasbhari set is returned.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 *
 * @return              Returns the final rasbhari patternset output.
 */
rasbhari rasb_implement::_hillclimb_iterative(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::ImproveMode = 2;
    rasb_opt::Sens = false;
    rasbhari RasbSet;
    if(rasb_opt::InFile.size() != 0){
        RasbSet = rasbhari(rasb_opt::InFile);
    }
    else{
        RasbSet = rasbhari(Size,Weight,MinDontCare,MaxDontCare);
    }
    rasb_opt::Size = RasbSet.pattern_set().size();
    rasb_opt::Weight = RasbSet.pattern_set().max_weight();
    rasb_opt::MinDontcare = RasbSet.pattern_set().min_dontcare();
    rasb_opt::MaxDontcare = RasbSet.pattern_set().max_dontcare();
    unsigned FillSizeBegin = 0, FillSizeEnd = 0;
    if(MaxDontCare + Weight > 13){
        FillSizeBegin = (MaxDontCare+Weight-13)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-13)%2 != 0){
            FillSizeEnd++;
        }
    }
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Initial Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    RasbSet.iterate_hill_climbing();
    if(MaxDontCare + Weight > 15){
        FillSizeBegin = (MaxDontCare+Weight-15)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-15)%2 != 0){
            FillSizeEnd++;
        }
    }
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Optimised Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    return RasbSet;
}
/**
 * Wrapper for the iterative hill-climbing process for the sensitivity;
 * depending on boolean output variables, some information will be printed to
 * std::cout. The final rasbhari set is returned.
 *
 * @param Size          The number of pattern.
 *
 * @param Weight        The number of match positions in the patterns, i.e.
 *                          the weight.
 *
 * @param MinDontCare   The minimal number of don't-care positions in the
 *                          patterns.
 *
 * @param MaxDontCare   The maximal number of don't-care positions in the
 *                          patterns.
 *
 * @return              Returns the final rasbhari patternset output.
 */
rasbhari rasb_implement::_hillclimb_sens(unsigned Size, unsigned Weight, unsigned MinDontCare, unsigned MaxDontCare){
    rasb_opt::ImproveMode = 3;
    rasb_opt::Sens = true;
    rasbhari RasbSet;
    if(rasb_opt::InFile.size() != 0){
        RasbSet = rasbhari(rasb_opt::InFile);
    }
    else{
        RasbSet = rasbhari(Size,Weight,MinDontCare,MaxDontCare);
    }
    rasb_opt::Size = RasbSet.pattern_set().size();
    rasb_opt::Weight = RasbSet.pattern_set().max_weight();
    rasb_opt::MinDontcare = RasbSet.pattern_set().min_dontcare();
    rasb_opt::MaxDontcare = RasbSet.pattern_set().max_dontcare();
    unsigned FillSizeBegin = 0, FillSizeEnd = 0;
    if(MaxDontCare + Weight > 13){
        FillSizeBegin = (MaxDontCare+Weight-13)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-13)%2 != 0){
            FillSizeEnd++;
        }
    }
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Initial Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    if(MaxDontCare + Weight > 15){
        FillSizeBegin = (MaxDontCare+Weight-15)/2;
        FillSizeEnd = FillSizeBegin;
        if((MaxDontCare+Weight-15)%2 != 0){
            FillSizeEnd++;
        }
    }
    RasbSet.hill_climbing_sensitivity();
    if(!rasb_opt::Silent){
        std::cout << " #" << std::string(FillSizeBegin,'=') << " Optimised Set " << std::string(FillSizeEnd,'=') << "#" << std::endl;
        RasbSet.print();
    }
    return RasbSet;
}