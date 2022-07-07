//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_BIT_H
#define M7_UTIL_BIT_H

#include <x86intrin.h>
#include "M7_lib/defs.h"
#include "Tag.h"


/**
 * utilities relating to the manipulation of integral types as bitsets
 */
namespace bit {
    /**
     * masks whose ith elements retain i bits when bitwise and-ed with another value of the same size (1, 2, 4 or 8 bytes)
     */
    template<typename T>
    using trunc_mask_array_t = std::array<T, sizeof(T)*8+1>;

    static constexpr trunc_mask_array_t<uint8_t> c_trunc_mask_8 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff};
    static constexpr trunc_mask_array_t<uint16_t> c_trunc_mask_16 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff};
    static constexpr trunc_mask_array_t<uint32_t> c_trunc_mask_32 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff,
             0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff,
             0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff,
             0x3fffffff, 0x7fffffff, 0xffffffff};
    static constexpr trunc_mask_array_t<uint64_t> c_trunc_mask_64 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff,
             0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff,
             0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff,
             0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff, 0x1ffffffff,
             0x3ffffffff, 0x7ffffffff, 0xfffffffff, 0x1fffffffff,
             0x3fffffffff, 0x7fffffffff, 0xffffffffff, 0x1ffffffffff,
             0x3ffffffffff, 0x7ffffffffff, 0xfffffffffff, 0x1fffffffffff,
             0x3fffffffffff, 0x7fffffffffff, 0xffffffffffff,
             0x1ffffffffffff, 0x3ffffffffffff, 0x7ffffffffffff,
             0xfffffffffffff, 0x1fffffffffffff, 0x3fffffffffffff,
             0x7fffffffffffff, 0xffffffffffffff, 0x1ffffffffffffff,
             0x3ffffffffffffff, 0x7ffffffffffffff, 0xfffffffffffffff,
             0x1fffffffffffffff, 0x3fffffffffffffff, 0x7fffffffffffffff,
             0xffffffffffffffff};


    static uint8_t get_trunc_element(uint_t n, tag::Int<1> /*size*/) {
        return c_trunc_mask_8[n];
    }
    static uint16_t get_trunc_element(uint_t n, tag::Int<2> /*size*/) {
        return c_trunc_mask_16[n];
    }
    static uint32_t get_trunc_element(uint_t n, tag::Int<4> /*size*/) {
        return c_trunc_mask_32[n];
    }
    static uint64_t get_trunc_element(uint_t n, tag::Int<8> /*size*/) {
        return c_trunc_mask_64[n];
    }

    /**
     * clear the bit in x at position i
     */
    template<typename T>
    static inline void clr(T &x, uint_t i) {
        x &= ~(T(1) << T(i));
    }
    /**
     * set the bit in x at position i
     */
    template<typename T>
    static inline void set(T &x, uint_t i) {
        x |= (T(1) << T(i));
    }
    /**
     * return true if the bit in x at position i is set, else false
     */
    template<typename T>
    static inline bool get(const T &x, uint_t i) {
        return x & (T(1) << T(i));
    }
    /**
     * use architecture dependent intrinsic to go to the next set bit position without looping
     * TODO: generalize for the absence of x86 instruction set with the BMI, SSE4.2, ABM etc (the exact set to which
     *  these instructions are taken to belong is vendor specific) extension e.g. ARM
     */
    inline uint_t next_setbit(unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        uint_t result = __tzcnt_u64(work);
        bit::clr(work, result);
        return result;
    }

    inline uint_t next_setbit(unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        uint_t result = __tzcnt_u64(work);
        bit::clr(work, result);
        return result;
    }

    inline uint_t next_setbit(unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        uint_t result = __tzcnt_u32(work);
        bit::clr(work, result);
        return result;
    }
    /**
     * use architecture dependent intrinsic to count the number of set bits in the word
     * TODO: see next_setbit
     */
    inline uint_t nsetbit(const unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    inline uint_t nsetbit(const unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    inline uint_t nsetbit(const unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        return _popcnt32(work);
    }

    /**
     * return v with the first n bits cleared
     */
    template<typename T>
    static T truncate(const T &v, uint_t n) {
        static_assert(std::is_integral<T>(), "can only bit-truncate integers");
        return v & get_trunc_element(n, tag::Int<sizeof(T)>());
    }

    /**
     * make a mask which is set (1) in the range [ibegin, iend), and clear (0) elsewhere
     */
    template<typename T>
    static void make_range_mask(T &v, uint_t ibegin, uint_t iend) {
        static_assert(std::is_integral<T>(), "can only make an integer-typed mask");
        v = get_trunc_element(iend, tag::Int<sizeof(T)>()) & ~get_trunc_element(ibegin, tag::Int<sizeof(T)>());
    }

    template<typename T>
    static T make_range_mask(uint_t ibegin, uint_t iend) {
        T v;
        make_range_mask(v, ibegin, iend);
        return v;
    }

    /**
     * apply an arbitrary bit mask to a single buffer word
     * e.g. for a range mask:
     *  F:    10101010101010110
     *  M:    00000111111000000
     *  F|M:  10101010101010110
     *  F&~M: 10101000000010110
     *
     * @tparam T
     *  buffer type of the bit field (integer of any size)
     * @param v
     *  buffer word
     * @param mask
     *  masking word
     * @param set
     *  true if bits that are set in the mask are to be set, else they are cleared
     */
    template<typename T>
    void apply_mask(T &v, const T &mask, bool set) {
        if (set) v |= mask;
        else v &= ~mask;
    }

    /**
     * make the range mask word then apply it
     */
    template<typename T>
    void apply_mask(T &v, uint_t ibegin, uint_t iend, bool set) {
        T mask;
        make_range_mask(mask, ibegin, iend);
        apply_mask(v, mask, set);
    }

    template<typename T>
    static inline uint_t nsetbit_before(const T &v, uint_t n) {
        return nsetbit(truncate(v, n));
    }

    template<typename T>
    static str_t to_string(const T &v) {
        str_t tmp;
        for (uint_t i = 0ul; i < sizeof(T) * CHAR_BIT; ++i) tmp += get(v, i) ? '1' : '0';
        return tmp;
    }

}

#endif //M7_UTIL_BIT_H
