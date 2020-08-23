//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include <src/core/list/PerforableMappedList.h>
#include <src/core/table/DeterminantField.h>
#include <src/core/table/FlagField.h>
#include <src/core/parallel/Reducible.h>
#include <list>

struct WalkerList : public PerforableMappedList<DeterminantElement> {
    DeterminantField m_determinant;
    NumericField<defs::wf_t> m_weight;
    NumericField<defs::ham_comp_t> m_hdiag;

private:
    struct Flags : public FlagField {
        Flag m_initiator, m_reference_connection, m_deterministic;

        Flags(Table *table, size_t nelement = 1, const std::string &description = "") :
            FlagField(table, nelement, description),
            m_initiator(this, nelement, "Initiator status"),
            m_reference_connection(this, nelement, "Connected to reference determinant"),
            m_deterministic(this, nelement, "In deterministic subspace") {}
    };

public:
    Flags m_flags;

    WalkerList(std::string name, size_t nsite, size_t nbucket) :
        PerforableMappedList<DeterminantElement>(name, m_determinant, nbucket),
        m_determinant(this, 1, nsite, "Determinant"),
        m_weight(this, 1, "Weight"),
        m_hdiag(this, 1, "Diagonal Hamiltonian matrix element"),
        m_flags(this, 1, "Flags") {
        ASSERT(m_determinant.element_dsize()==
            integer_utils::divceil(2*nsite, CHAR_BIT * sizeof(defs::data_t)))
    }

    ~WalkerList() {

    }

    void report_top_weighted(){
        const size_t n = 15;
        std::cout << "Top weighted configurations:" << std::endl;
        auto top = highest_weighted_row_inds_local(n);
        size_t i = 0ul;
        for (auto iter=top.begin(); iter!=top.end(); iter++) {
            if (i++==n) break;
            std::cout << m_determinant(*iter).to_string() << " "
            << utils::num_to_string(*m_weight(*iter)) << " "
            << (m_flags.m_initiator(*iter) ? "*":"") << std::endl;
        }
    }

    defs::wf_comp_t l1_norm(const size_t &ielement=0) {
        Reducible<defs::wf_comp_t> tot;
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            tot += std::abs(*m_weight(irow));
        }
        tot.mpi_sum();
        ASSERT(tot.reduced()==tot.reduced()); // else NaN
        ASSERT(tot.reduced() > 0);
        return tot.reduced();
    }

    defs::wf_comp_t square_norm(const size_t &ielement=0) {
        Reducible<defs::wf_comp_t> tot;
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            tot += std::pow(std::abs(*m_weight(irow)), 2.0);
        }
        tot.mpi_sum();
        ASSERT(tot.reduced()==tot.reduced()); // else NaN
        ASSERT(tot.reduced() > 0);
        return tot.reduced();
    }

    void normalize(const size_t &ielement=0, const defs::wf_t &norm=1.0) {
        auto factor = std::sqrt(square_norm(ielement));
        for (size_t irow = 0ul; irow < high_water_mark(0); irow++) {
            m_weight(irow)*=norm/factor;
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

    /**
     * @param n
     *      maximum number of rows to find
     * @return
     *      ordered list no longer than n elements of row indices with the
     *      largest weight magnitudes
     */
    std::list<size_t> highest_weighted_row_inds_local(size_t n) {
        std::list<size_t> result;
        for (size_t irow = 0ul; irow < high_water_mark(0); ++irow) {
            const auto abs_weight = std::abs(*m_weight(irow));
            auto iter = result.begin();
            size_t iiter=0ul;
            while (iter != result.end() && iiter < n){
                if (abs_weight>std::abs(*m_weight(*iter))) break;
                iter++; iiter++;
            }
            result.insert(iter, irow);
        }
        return result;
    }

    /**
     * @param row_inds
     *      will contain indices of the walker list local to
     * @param ndet_total
     *      the number of determinants sought in total
     * @param ndet_local
     *      the number of determinants sought locally before gatherv-ing,
     *      this is a memory-efficiency measure, since ndet_total may be
     *      large. Setting ndet_local to ndet_total will succeed in one
     *      call, but setting a small value of ndet_local may entail a
     *      subsequent recursive call with a larger value
     */
    void highest_weighted_row_inds(defs::inds& row_inds, size_t ndet_total, size_t ndet_local){
        auto local_irow_list = highest_weighted_row_inds_local(ndet_local);
        std::vector<size_t> local_irows(local_irow_list.begin(), local_irow_list.end());
        std::vector<defs::wf_comp_t> local_weights;
        local_weights.reserve(local_irows.size());
        for (auto & local_irow : local_irows){
            local_weights.push_back(*m_weight(local_irow));
        }
        std::vector<size_t> global_irows;
        mpi::all_gatherv(local_irows, global_irows);
        mpi::all_gatherv(local_irows, global_irows);
    }

};


#endif //M7_WALKERLIST_H
