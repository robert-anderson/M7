//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "ExcitGen.h"

struct UniformSingles : public FrmExcitGen {
    using FrmExcitGen::draw;

    UniformSingles(const Hamiltonian& ham, PRNG& prng);

    bool draw(const size_t& exsig, const field::FrmOnv &src, CachedOrbs& orbs,
               defs::prob_t &prob, conn::FrmOnv &conn) override;

    bool draw(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob, defs::ham_t &helem,
              conn::FrmOnv &conn) override;

    size_t approx_nconn() const override;

    std::string description() const override;

};


#endif //M7_UNIFORMSINGLES_H
