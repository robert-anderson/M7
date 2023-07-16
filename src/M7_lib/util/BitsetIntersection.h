//
// Created by rja on 16/07/23.
//

#ifndef M7_BITSETINTERSECTION_H
#define M7_BITSETINTERSECTION_H

#include "M7_lib/parallel/MPIAssert.h"
#include "M7_lib/util/Bit.h"
#include "M7_lib/util/Hash.h"
#include "Integer.h"

/**
 * a namespace of functions to efficiently compute and iterate over bitwise AND of vector bitsets
 */
namespace bitset_isect {

    /**
     * @param isetbits
     *  unordered indices to be converted into a bitset
     * @param nbit
     *  number of bits in the resulting vector bitset
     * @return
     *  vector bitset with 1 at positions corresponding to values present in isetbits, and 0 otherwise
     */
    uintv_t make_bitset(const uintv_t& isetbits, uint_t nbit);

    /**
     * as above, but infer bitset length
     */
    uintv_t make_bitset(const uintv_t& isetbits);

    /**
     * compute the intersection on the word range [iword_begin, iword_end) of two bitsets as a sparse intersection
     * vector "isect" whose elements are of the form:
     *  {non-zero word index, bitwise AND of the words at that index in bitset1 and bitset2}
     * @param bitset1
     *  a vector bitset
     * @param bitset2
     *  another vector bitset
     * @param iword_begin
     *  index of the first word to include in the intersection
     * @param iword_end
     *  index of the first word after iword_begin to exclude in the intersection
     * @param isect
     *  sparse intersection vector
     */
    void make_isect(const uintv_t& bitset1, const uintv_t& bitset2, uint_t iword_begin, uint_t iword_end, v_t<uintp_t>& isect);

    /**
     * as above, but do not restrict the word range
     */
    void make_isect(const uintv_t& bitset1, const uintv_t& bitset2, v_t<uintp_t>& isect);

    /**
     * as above, but return the intersection vector by value
     */
    v_t<uintp_t> make_isect(const uintv_t& bitset1, const uintv_t& bitset2);

    /**
     * as above, but take intersection of the bitset with itself (i.e. convert to intersection vector)
     */
    v_t<uintp_t> make_isect(const uintv_t& bitset);

    /**
     * as above, but place result into referenced intersection vector
     */
    void make_isect(const uintv_t& bitset, v_t<uintp_t>& isect);

    /**
     * update isect in-place by taking its intersection with the given bitset
     */
    void isect(v_t<uintp_t>& isect, const uintv_t& bitset);

    /**
     * compute intersection between isect_in and bitset and store the result in isect_out
     */
    void make_isect(const v_t<uintp_t>& isect_in, const uintv_t& bitset, v_t<uintp_t>& isect_out);

    /**
     * @tparam fn_t
     *  callable type which accepts a uint_t
     * @param isect
     *  intersection vector
     * @param fn
     *  function to call each set bit in the intersection vector
     */
    template<typename fn_t>
    void foreach_in_isect(const v_t<uintp_t>& isect, const fn_t& fn) {
        functor::assert_prototype<void(uint_t)>(fn);
        for (auto& pair: isect) {
            auto work = pair.second;
            while (work) {
                const auto ibit = bit::next_setbit(work);
                fn(pair.first * sizeof(uint_t) * CHAR_BIT + ibit);
            }
        }
    }

    /**
     * get intersection vector as std::set instance for testing / debugging
     */
    std::set<uint_t> to_std_set(const v_t<uintp_t>& isect);

    /**
     * get intersection vector as std::vector instance for testing / debugging
     */
    uintv_t to_std_vector(const v_t<uintp_t>& isect);

    /**
     * recurse through sets for intersection until maximum number of sets intersected is reached with the ordering that
     * each selected set has lower index than the previously intersected set
     * @tparam fn_t
     *  callable type which accepts a vector of set indices (in strict ascending order) and a sparse intersection vector
     * @param isets
     *  current set selection in strict ascending order whose size is the current number of sets intersected and whose
     *  capacity is the maximum number of sets to compute intersections of
     * @param isect
     *  result of taking the intersection of all sets indexed in isets as a sparse intersection vector
     * @param bitsets
     *  all sets in bitset form
     * @param fn
     *  function to be called each time a non-empty set is formed with n sets with n in [1, isets.capacity()]
     */
    template<typename fn_t>
    void foreach_unique_one_level(uintv_t& isets, const v_t<uintp_t>& isect, const v_t<uintv_t>& bitsets, const fn_t& fn) {
        functor::assert_prototype<void(const uintv_t& /*isets*/, const v_t<uintp_t>& /*isect*/)>(fn);
        // handle the current set combination which have intersection isect
        fn(isets, isect);
        // capacity of the vector gives the maximum number of sets to consider intersecting
        if (isets.size() == isets.capacity()) return;
        // set index of this level is constrained to be less than that of the next most senior level
        const auto iset_max = isets.back();
        v_t<uintp_t> next_isect;
        // add another level to the vector or set indices
        isets.push_back(~0ul);
        for (isets.back()=0ul; isets.back() < iset_max; ++isets.back()) {
            // compute intersection and store result in next_isect
            bitset_isect::make_isect(isect, bitsets[isets.back()], next_isect);
            // if there are any elements in the result, go another level deeper
            if (!next_isect.empty()) foreach_unique_one_level(isets, next_isect, bitsets, fn);
        }
        // revert to initial level
        isets.resize(isets.size()-1);
    }

    /**
     * loop over top-level isects (simply the bitsets converted to sparse intersection vector form) and dispatch full
     * recursive enumeration of non-empty intersections
     * @tparam fn_t
     *  callable type which accepts a vector of set indices (in strict ascending order) and a sparse intersection vector
     * @param top_isects
     *  top-level intersections. usually this will correspond to the full bitsets converted to intersection vector form
     *  but not necessarily - this could hold a subset of bits unique to each MPI rank as an efficient means of
     *  achieving parallel partitioning
     * @param nmax
     *  maximum number of sets between which to form intersection
     * @param fn
     *  function to be called each time a non-empty set is formed with n sets with n in [1, isets.capacity()]
     */
    template<typename fn_t>
    void foreach_unique(const v_t<uintv_t>& bitsets, const v_t<v_t<uintp_t>>& top_isects, uint_t nmax, const fn_t& fn) {
        functor::assert_prototype<void(const uintv_t& /*isets*/, const v_t<uintp_t>& /*isect*/)>(fn);
        const auto nset = bitsets.size();
        DEBUG_ASSERT_EQ(top_isects.size(), nset, "incompatible number of top-level intersections");
        uintv_t isets;
        isets.reserve(nmax);
        for (uint_t iset=0ul; iset < nset; ++iset) {
            isets.clear();
            isets.push_back(iset);
            const auto& isect = top_isects[iset];
            foreach_unique_one_level(isets, isect, bitsets, fn);
        }
    }
}


#endif //M7_BITSETINTERSECTION_H
