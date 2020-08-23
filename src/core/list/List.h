//
// Created by Robert John Anderson on 2020-03-31.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include <src/core/parallel/MPIWrapper.h>
#include <src/core/sort/QuickSorter.h>
#include <list>
#include "src/core/table/Table.h"
#include "src/core/table/NumericField.h"

class List : public Table {
    List *m_recv = nullptr;
    defs::inds m_high_water_mark;
public:
    explicit List(std::string name, size_t nsegment = 1);

    void recv(List *list);

    void expand(size_t delta_nrow) override;

    const defs::inds &high_water_mark() const;

    const size_t &high_water_mark(const size_t isegment) const;

    size_t push(const size_t &isegment = 0, const size_t &nrow = 1);

    size_t expand_push(const size_t &isegment = 0, const size_t &nrow = 1, double factor = 1.5);

    void zero() override;

    virtual std::string to_string() const;

    virtual std::string to_string(const defs::inds &nrows) const;

    void communicate();

    void all_gather(List &local);

    /**
     * @param n
     *      maximum number of rows to find
     * @return
     *      ordered list no longer than n elements of row indices with the
     *      largest weight magnitudes
     */
    template<typename T>
    std::list<size_t> top_row_inds_local(const NumericField<T> &sort_field, size_t n, bool abs_mag=true) {
        if(!owns_field(&sort_field)) throw std::runtime_error("sort field must belong to List being selected from");
        std::list<size_t> result;
        /*

for (size_t irow = 0ul; irow < high_water_mark(0); ++irow) {
    const auto abs_weight = std::abs(*m_weight(irow));
    auto iter = result.begin();
    size_t iiter = 0ul;
    while (iter != result.end() && iiter < n) {
        if (abs_weight > std::abs(*m_weight(*iter))) break;
        iter++;
        iiter++;
    }
    result.insert(iter, irow);
} */

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
     *      subsequent call with a larger value
     * @returns
     *      true if successful
     */
    bool highest_weighted_row_inds(defs::inds &row_inds, size_t ndet_global, size_t ndet_local) {
        return true;
//        auto local_irow_list = highest_weighted_row_inds_local(ndet_local);
//        std::vector<size_t> local_irows(local_irow_list.begin(), local_irow_list.end());
//        std::vector<defs::wf_comp_t> local_weights, global_weights;
//        local_weights.reserve(local_irows.size());
//        for (auto &local_irow : local_irows) {
//            local_weights.push_back(*m_weight(local_irow));
//        }
//        std::vector<size_t> global_irows;
//        defs::inds recvcounts, recvdispls;
//        mpi::all_gatherv(local_irows, global_irows, recvcounts, recvdispls);
//        mpi::all_gatherv(local_weights, global_weights, recvcounts, recvdispls);
//        /*
//         * now we have top weights from all ranks, we need to find the ndet_global top weighted
//         * determinants. If the number of dets selected from irank is equal to recvcounts[irank],
//         * we need to try again, since there could have been weights that did not make the cut on
//         * the node, but would be included in the global selection.
//         */
//        defs::inds nrank_select(mpi::nrank()); // number of selections per rank
//        row_inds.reserve(ndet_global);
//        for (auto &it : row_inds) {
//            defs::wf_comp_t top_weight = 0.0;
//            size_t top_irank = ~0ul;
//            for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
//                auto ielement = recvdispls[irank] + nrank_select[irank];
//                auto irow = global_irows[ielement];
//                auto weight = global_weights[ielement];
//                if (weight > top_weight) {
//                    top_weight = weight;
//                    top_irank = irank;
//                }
//            }
//            ASSERT(top_irank != ~0ul)
//            nrank_select[top_irank]++;
//            if (nrank_select[top_irank] == ndet_local) {
//                /*
//                 * unsuccessful with specified ndet_local
//                 */
//                return false;
//            }
//            it =
//        }
//        return true;
    }


};

#endif //M7_LIST_H
