//
// Created by anderson on 27/07/2022.
//

#ifndef M7_DISTRIBUTION_H
#define M7_DISTRIBUTION_H

#include "Redistributor.h"

class Distribution {

    uintv_t m_block_iranks;

    /**
     * number of blocks stored on this MPI rank
     */
    uint_t m_nblock_local = 0ul;

public:
    uint_t nblock() const {
        return m_block_iranks.size();
    }

    uint_t nblock_() const {
        return m_nblock_local;
    }

    const uintv_t& block_iranks() const {
        return m_block_iranks;
    }

    explicit Distribution(size_t nblock);

    void update(const Redistributor& redist);

    template<typename field_t>
    uint_t iblock(const field_t& field) const {
        return field.hash()%m_block_iranks.size();
    }

    template<typename field_t>
    uint_t irank(const field_t& field) const {
        return m_block_iranks[iblock(field)];
    }
};

#endif //M7_DISTRIBUTION_H
