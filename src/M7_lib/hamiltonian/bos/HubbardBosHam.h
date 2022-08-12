//
// Created by anderson on 19/07/2022.
//

#ifndef M7_HUBBARDBOSHAM_H
#define M7_HUBBARDBOSHAM_H

#include "BosHam.h"

class HubbardBosHam : public BosHam {
    /**
     * on-site repulsion scalar in units of the hopping
     */
    const ham_t m_u;

public:
    HubbardBosHam(ham_t u, const std::shared_ptr<lattice::Lattice>& lattice);

    HubbardBosHam(opt_pair_t opts);

    ham_t get_coeff_0011(uint_t a, uint_t i) const override;

    ham_t get_coeff_0022(uint_t a, uint_t b, uint_t i, uint_t j) const override;

    ham_t get_element_0000(const field::BosOnv& onv) const override;

    ham_t get_element_0011(const field::BosOnv& onv, const conn::BosOnv& conn) const override;

    ham_t get_element_0022(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const override {
        return 0;
    }

    void log_data() const override;

    uint_t default_nboson() const override;

    excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& opts) const override;

    conn_foreach_list_t make_foreach_iters() const override;
};


#endif //M7_HUBBARDBOSHAM_H
