//
// Created by rja on 08/12/2020.
//

#ifndef M7_QUICKSORTER_H
#define M7_PARALLELEXTREMALVALUES_H

#include "TableExtremalValues.h"
#include "src/core/table/BufferedTableArray.h"

/**
 * Finds the extremal values over all processes, using distributed instances
 * of the binary heap-based ExtremalIndices class. The strategy for finding N
 * extremal-valued row indices is to find N/nrank on
 * each rank. This assumes a roughly equal distribution of extremal values
 * among MPI ranks.
 *
 * The values are then gathered, and the root process "hands out" indices to
 * the processes by finding the next extreme value
 *
 * as soon as all N/nrank locally-extremal values on irank have
 * been added to the globally-extremal selection, another find operation is
 * required on irank before the procedure can continue
 */
template<typename viewable_t>
class ParallelExtremalValues {
    TableExtremalValues<viewable_t> m_local_xvs;
    /**
     * only used for the initial gathering of inds,
     * subsequent communication of local extremal values is done point-to-point
     */
    defs::inds m_gathered_inds;


    defs::inds m_gathered_values;
    bool m_max, m_abs_val;
public:
    ParallelExtremalValues(const viewable_t &viewable, bool max = true, bool abs_val = false) :
            m_local_xvs(viewable, max, abs_val), m_max(max), m_abs_val(abs_val) {}

    void reset(const TableBase& table) {
        m_local_xvs.reset(table.m_hwm);
    }

    /**
     * @param[in] nfind - the number of extremal values to find globally
     */
    void find(size_t nfind) {
        auto nfind_local = std::max(nfind / mpi::nrank(), 1ul);
        m_local_xvs.find(nfind_local);
        defs::inds counts, displs;

        if (mpi::i_am_root()) {
            counts.resize(mpi::nrank());
            displs.resize(mpi::nrank());
        }
        mpi::gather(m_local_xvs.nfound(), counts, 0);
        if (mpi::i_am_root()) {
            mpi::counts_to_displs_consec(counts, displs);
            m_gathered_inds.resize(displs.back()+counts.back());
        }

//        if (mpi::i_am_root()) nfound.resize(mpi::nrank());
//        mpi::gather(m_local_xvs.nfound(), nfound, 0);
//        if (mpi::i_am_root()) m_gathered_inds.resize(nfind);

        mpi::gatherv(m_local_xvs.begin(), m_local_xvs.nfound(), m_gathered_inds.data(), counts, displs, 0);
        if (mpi::i_am_root()){
            utils::print(m_gathered_inds);
        }
    }

//    const size_t &nfound() const;
//
//    const size_t &operator[](const size_t &ifound) const;
//
//    void reset(size_t hwm);
//
//    void reset(const Table &table);
//
//    void find(size_t nfind);

};


#endif //M7_QUICKSORTER_H
