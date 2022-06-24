//
// Created by Robert J. Anderson on 25/11/2021.
//

#ifndef M7_BOSONSUMCONSERVINGDOUBLES_H
#define M7_BOSONSUMCONSERVINGDOUBLES_H

#include "BosExcitGen.h"

struct BosonSumConservingDoubles : public BosExcitGen {
    BosonSumConservingDoubles(const BosHam &h, PRNG &prng);

private:
    /**
     * @param min
     *  minimum allowed value for the draw of a
     * @param max
     *  maximum allowed value for the draw of a
     * @param nexclude
     *  index a should not be assigned the same value as i or j, step over these
     */
    void set_a_range(uint_t i, uint_t j, uint_t& min, uint_t &max, uint_t& nexclude) const;

    /**
     * @return
     *  number of possible values for index a
     */
    uint_t na(uint_t i, uint_t j) const;

public:

    bool draw_bos(uint_t /*exsig*/, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) override;

    defs::prob_t prob_bos(const field::BosOnv &src, const conn::BosOnv &conn) const override;

};


#endif //M7_BOSONSUMCONSERVINGDOUBLES_H
