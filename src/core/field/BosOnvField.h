//
// Created by rja on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include "NumberField.h"
#include "src/core/basis/BasisDims.h"

struct BosOnvField : NdNumberField<defs::bos_occ_t, 1> {

    BosOnvField(Row *row, size_t nmode, std::string name = "");

    BosOnvField(Row *row, BasisDims bd, std::string name = "");

    BosOnvField(const BosOnvField &other) :
            BosOnvField(other.row_of_copy(), {0ul, other.m_format.m_shape[0]}, other.m_name) {}

    BosOnvField &operator=(const BosOnvField &other) {
        NdNumberField<defs::bos_occ_t, 1>::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const defs::inds &inds) {
        DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
        for (size_t i = 0ul; i < inds.size(); ++i) (*this)[i] = inds[i];
        return *this;
    }

    /**
     * try to set using a value supplied by the user. raise an error if the input literal is nonempty and invalid
     * @param v
     *  user-defined literal
     * @return
     *  *this
     */
    BosOnvField &attempt_set_from_input(const defs::inds &v) {
        if (!v.empty()) {
            REQUIRE_EQ(v.size(), m_nelement, "Incorrectly sized boson ONV value given");
            *this = v;
        }
        return *this;
    }

};


#endif //M7_BOSONVFIELD_H
