//
// Created by rja on 08/12/2020.
//

#ifndef M7_BUFFEREDTABLEVECTOR_H
#define M7_BUFFEREDTABLEVECTOR_H

#include "BufferedTable.h"

/**
 * Simply a wrapper for a std::vector of BufferedTables. This differs from the
 * BufferedTableArray in that the tables in the structure do not share an underlying
 * data buffer, and thus they are not contiguously communicatable e.g. by MPI_Alltoallv
 *
 * This class is however superior to the BufferedTableArray for point-to-point comms,
 * since adding to one rank's corresponding element of the BufferedTableVector only
 * requires the reallocation of that Buffer, not the Buffer on which the entire structure
 * is stored
 */
template<typename table_t>
class BufferedTableVector {
    std::vector<BufferedTable<table_t>> m_bts;

    template<typename ...Args>
    BufferedTableVector(size_t nelement, Args&&... args){
        m_bts.reserve(nelement);
        for (size_t i=0ul; i<nelement; ++i) m_bts.template emplace_back(args);
    }

    BufferedTable<table_t>& operator[](const size_t& i) {
        return m_bts[i];
    }

    const BufferedTable<table_t>& operator[](const size_t& i) const {
        return m_bts[i];
    }

    void gatherv(const table_t& local){
        if (m_bts.size()!=mpi::nrank())
            mpi::abort("attempting to gather to a BufferedTableVector of incorrect size");
        mpi::gatherv()
    }
};


#endif //M7_BUFFEREDTABLEVECTOR_H
