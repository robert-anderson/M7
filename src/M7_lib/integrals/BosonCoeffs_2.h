//
// Created by Robert J. Anderson on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_2_H
#define M7_BOSONCOEFFS_2_H

#include <M7_lib/parallel/SharedArray.h>
#include <M7_lib/util/Integer.h>

class BosonCoeffs_2 {

    size_t index(const size_t &i, const size_t &j, const size_t& k, const size_t& l) const {
        using namespace utils::integer;
        auto ij = i*m_nmode+j;
        auto kl = k*m_nmode+l;
        return ij>=kl ? trigmap(ij, kl) : trigmap(kl, ij);
    }

    const size_t m_nmode;
public:
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_2(size_t nmode);

    void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;

    const defs::ham_t& phys_element(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;
};


#endif //M7_BOSONCOEFFS_2_H
