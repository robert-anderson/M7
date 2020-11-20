//
// Created by RJA on 20/11/2020.
//

#ifndef M7_WEIGHTEDDRAWER_H
#define M7_WEIGHTEDDRAWER_H

#include "PRNG.h"

class WeightedDrawer {

    void update_cumprobs();

    void set_one(size_t iprob, defs::prob_t prob);

    template<typename ...Args>
    void set_one(size_t iprob, defs::prob_t first, Args... rest){
        m_probs[iprob] = first;
        set_one(iprob+1, rest...);
    }

public:
    std::vector<defs::prob_t> m_probs;
    std::vector<defs::prob_t> m_cumprobs;
    PRNG& m_prng;

    WeightedDrawer(size_t nprob, PRNG& prng);

    template<typename ...Args>
    void set(Args... probs){
        ASSERT(sizeof...(Args)+1>=m_probs.size());
        set_one(0, probs...);
        update_cumprobs();
    }

    size_t draw();

    const defs::prob_t& prob(size_t iprob);

};


#endif //M7_WEIGHTEDDRAWER_H
