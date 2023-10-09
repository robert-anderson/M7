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

    template<typename field_t>
    static uintv_t one_irank_in_each_shmem_region(const field_t& field) {
        const auto hash = field.hash();
        uintv_t out(mpi::nshmem());
        uint_t ishmem = 0ul;
        for (auto& elem: out) {
            const auto i = hash % mpi::g_nrank_in_shmem_realms[ishmem];
            elem = mpi::g_iranks_world_in_shmem_realms[ishmem][i];
            ++ishmem;
        }
        return out;
    }
};

#endif //M7_DISTRIBUTION_H
