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

    typedef v_t<uintp_t> siv_t;
    typedef v_t<v_t<uintp_t>> vsiv_t;
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
     * compute the intersection on the word range [iword_begin, iword_end) of two bitsets as a sparse intersection
     * vector "siv" whose elements are of the form:
     *  {non-zero word index, bitwise AND of the words at that index in bitset1 and bitset2}
     * @param bitset1
     *  a vector bitset
     * @param bitset2
     *  another vector bitset
     * @param iword_begin
     *  index of the first word to include in the intersection
     * @param iword_end
     *  index of the first word after iword_begin to exclude in the intersection
     * @param siv
     *  sparse intersection vector
     */
    void isect(const uintv_t& bitset1, const uintv_t& bitset2, uint_t iword_begin, uint_t iword_end, siv_t& siv);

    /**
     * as above, but do not restrict the word range
     */
    void isect(const uintv_t& bitset1, const uintv_t& bitset2, siv_t& siv);

    /**
     * as above, but return the sparse intersection vector by value
     */
    siv_t isect(const uintv_t& bitset1, const uintv_t& bitset2);

    /**
     * as above, but take intersection of the bitset with itself (i.e. convert to sparse intersection vector)
     */
    siv_t bitset_to_siv(const uintv_t& bitset);

    /**
     * as above, but place result into referenced sparse intersection vector
     */
    void bitset_to_siv(const uintv_t& bitset, siv_t& siv);

    /**
     * repeat above process for many pairs of bitsets
     */
    vsiv_t isect_many(const v_t<uintv_t>& bitsets1, const v_t<uintv_t>& bitsets2, uint_t iword_begin, uint_t iword_end);
    /**
     * as above but take intersections of the bitsets with themselves (i.e. convert to intersection vectors)
     */
    vsiv_t bitset_to_siv_many(const v_t<uintv_t>& bitsets, uint_t iword_begin, uint_t iword_end);
    /**
     * as above, but infer the begin and end words from the size in words of the first bitset in the given vector
     */
    vsiv_t bitset_to_siv_many(const v_t<uintv_t>& bitsets);


    /**
     * update siv in-place by taking its intersection with the given bitset
     */
    void isect(siv_t& siv, const uintv_t& bitset);

    /**
     * compute intersection between siv_in and bitset and store the result in siv_out
     */
    void isect(const siv_t& siv_in, const uintv_t& bitset, siv_t& siv_out);

    /**
     * @tparam fn_t
     *  callable type which accepts a uint_t
     * @param siv
     *  intersection vector
     * @param fn
     *  function to call each set bit in the intersection vector
     */
    template<typename fn_t>
    void foreach_in_isect(const siv_t& siv, const fn_t& fn) {
        functor::assert_prototype<void(uint_t)>(fn);
        for (auto& pair: siv) {
            auto work = pair.second;
            while (work) {
                const auto ibit = bit::next_setbit(work);
                fn(pair.first * sizeof(uint_t) * CHAR_BIT + ibit);
            }
        }
    }

    /**
     * get sparse intersection vector as std::set instance for testing / debugging
     */
    std::set<uint_t> to_std_set(const siv_t& siv);

    /**
     * get sparse intersection vector as std::vector instance for testing / debugging
     */
    uintv_t to_std_vector(const siv_t& siv);

    /**
     * recurse through sets for intersection until maximum number of sets intersected is reached with the ordering that
     * each selected set has lower index than the previously intersected set
     * @tparam fn_t
     *  callable type which accepts a vector of set indices (in strict ascending order) and a sparse intersection vector
     * @param isets
     *  current set selection in strict ascending order whose size is the current number of sets intersected and whose
     *  capacity is the maximum number of sets to compute intersections of
     * @param siv
     *  result of taking the intersection of all sets indexed in isets as a sparse intersection vector
     * @param bitsets
     *  all sets in bitset form
     * @param max_only
     *  if true, only call the fn on a intersection when the size and capacity of the isets vector are identical
     * @param fn
     *  function to be called each time a non-empty set is formed with n sets with n in [1, isets.capacity()]
     */
    template<typename fn_t>
    void foreach_unique_one_level(uintv_t& isets, const siv_t& siv, const v_t<uintv_t>& bitsets, bool max_only, const fn_t& fn) {
        functor::assert_prototype<void(const uintv_t& /*isets*/, const siv_t& /*siv*/)>(fn);
        // handle the current set combination which have intersection siv
        if (!max_only || (isets.size() == isets.capacity())) fn(isets, siv);
        // capacity of the vector gives the maximum number of sets to consider intersecting
        if (isets.size() == isets.capacity()) return;
        const auto size = isets.size();
        // set index of this level is constrained to be less than that of the next most senior level
        const auto iset_max = isets.back();
        // total number of sets from which the unique intersections are being prepared
        const auto nset = bitsets.size();
        siv_t next_siv;
        // add another level to the vector or set indices
        isets.push_back(iset_max+1);
        for (; isets.back() < nset; ++isets.back()) {
            // compute intersection and store result in next_siv
            bitset_isect::isect(siv, bitsets[isets.back()], next_siv);
            // if there are any elements in the result, go another level deeper
            if (!next_siv.empty()) foreach_unique_one_level(isets, next_siv, bitsets, max_only, fn);
        }
        // revert to initial level
        isets.resize(size);
    }

    /**
     * loop over top-level sivs (simply the bitsets converted to sparse intersection vector form) and dispatch full
     * recursive enumeration of non-empty intersections
     * @tparam fn_t
     *  callable type which accepts a vector of set indices (in strict ascending order) and a sparse intersection vector
     * @param top_sivs
     *  top-level intersections. usually this will correspond to the full bitsets converted to intersection vector form
     *  but not necessarily - this could hold a subset of bits unique to each MPI rank as an efficient means of
     *  achieving parallel partitioning
     * @param nmax
     *  maximum number of sets between which to form intersection
     * @param max_only
     *  if true, only call the fn on a intersection of exactly nmax sets
     * @param fn
     *  function to be called each time a non-empty set is formed with n sets with n in [1, isets.capacity()]
     */
    template<typename fn_t>
    void foreach_unique(const v_t<uintv_t>& bitsets, const vsiv_t& top_sivs, uint_t nmax, bool max_only, const fn_t& fn) {
        functor::assert_prototype<void(const uintv_t& /*isets*/, const siv_t& /*siv*/)>(fn);
        const auto nset = bitsets.size();
        DEBUG_ASSERT_EQ(top_sivs.size(), nset, "incompatible number of top-level intersections");
        uintv_t isets;
        isets.reserve(nmax);
        for (uint_t iset=0ul; iset < nset; ++iset) {
            isets.clear();
            isets.push_back(iset);
            const auto& siv = top_sivs[iset];
            foreach_unique_one_level(isets, siv, bitsets, max_only, fn);
        }
    }
}


#endif //M7_BITSETINTERSECTION_H
