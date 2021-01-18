//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/RankAllocator.h"

struct WalkerTable : public MappedTable<fields::Onv<>> {
    fields::Onv<> m_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    fields::Number<defs::ham_comp_t> m_hdiag;
private:
    struct WalkerTableFlagSet : FlagSet {
        Flags<defs::ndim_wf> m_initiator;
        Flag m_reference_connection;
        Flag m_deterministic;

        WalkerTableFlagSet(fields::Bitset *bitset, size_t nroot, size_t nreplica) :
                FlagSet(bitset),
                m_initiator(this, "is initiator", nroot, nreplica),
                m_reference_connection(this, "is connected to reference ONV"),
                m_deterministic(this, "is in deterministic subspace") {}
    };
public:
    fields::Flags<WalkerTableFlagSet> m_flags;

    WalkerTable(size_t nbucket, size_t nsite, size_t nroot, size_t nreplica) :
            MappedTable<fields::Onv<>>(m_onv, nbucket),
            m_onv(this, nsite, "occupation number vectors"),
            m_weight(this, "weights", nroot, nreplica),
            m_hdiag(this, "hamiltonian diagonal element"),
            m_flags(this, "flags", nroot, nreplica),{}

    WalkerTable(const Options &opts, size_t nsite, ra_t& ra):
    WalkerTable(nbucket_guess(opts.nwalker_target, 4), nsite, opts.nroot, opts.nreplica, ra){}

#if 0

    std::list<size_t> highest_weighted_row_inds_local(size_t n) {
        return List::top_row_inds_local(m_weight, n);
    }

    void report_top_weighted() {
        const size_t n = 15;
        std::cout << "Top weighted configurations:" << std::endl;
        auto top = highest_weighted_row_inds_local(n);
        size_t i = 0ul;
        for (auto iter = top.begin(); iter != top.end(); iter++) {
            if (i++ == n) break;
            std::cout << m_determinant(*iter).to_string() << " "
                      << utils::num_to_string(*m_weight(*iter)) << " "
                      << (m_flags.m_initiator(*iter) ? "*" : "") << std::endl;
        }
    }

    defs::wf_comp_t l1_norm(const size_t &ielement = 0) {
        Reducible<defs::wf_comp_t> tot;
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            tot += std::abs(*m_weight(irow));
        }
        tot.mpi_sum();
        ASSERT(tot.reduced() == tot.reduced()); // else NaN
        ASSERT(tot.reduced() > 0);
        return tot.reduced();
    }

    defs::wf_comp_t square_norm(const size_t &ielement = 0) {
        Reducible<defs::wf_comp_t> tot;
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            tot += std::pow(std::abs(*m_weight(irow)), 2.0);
        }
        tot.mpi_sum();
        ASSERT(tot.reduced() == tot.reduced()); // else NaN
        ASSERT(tot.reduced() > 0);
        return tot.reduced();
    }

    void normalize(const size_t &ielement = 0, const defs::wf_t &norm = 1.0) {
        auto factor = std::sqrt(square_norm(ielement));
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            m_weight(irow) *= norm / factor;
        }
    }

    using PerforableMappedList<DeterminantElement>::push;

    size_t add(const DeterminantElement &key, const defs::wf_t &weight, const defs::ham_comp_t &hdiag,
               bool initiator = false, bool reference_connection = false, bool deterministic = false) {
        auto irow = PerforableMappedList::push(key);
        m_weight(irow) = weight;
        m_hdiag(irow) = hdiag;
        m_flags.m_initiator(irow) = initiator;
        ASSERT(m_flags.m_initiator(irow) == true);
        m_flags.m_reference_connection(irow) = reference_connection;
        m_flags.m_deterministic(irow) = deterministic;
        return irow;
    }

    //debugging only
    size_t verify_ninitiator(const double &nadd) {
        size_t ninitiator = 0ul;
        for (size_t i = 0ul; i < m_nrow_per_segment; ++i) {
            //auto abs_weight = std::abs(*m_weight(i, 0));
            bool is_initiator = m_flags.m_initiator(i, 0);
            //ninitiator += (abs_weight > nadd);
            ninitiator += is_initiator;
            //if((abs_weight > nadd) != is_initiator) return ~0ul;
        }
        return ninitiator;
    }
#endif //M7_WALKERTABLE_H

};

#endif //M7_WALKERTABLE_H
