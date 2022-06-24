//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_BIT_H
#define M7_UTIL_BIT_H

#include <x86intrin.h>
#include "M7_lib/defs.h"

using namespace defs;

/**
 * utilities relating to the manipulation of integral types as bitsets
 */
namespace bit {
    /**
     * masks whose ith elements retain i bits when bitwise and-ed with another value of the same size (4 or 8 bytes)
     */
    static constexpr std::array<uint8_t, 9> c_trunc_mask_8 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff};
    static constexpr std::array<uint16_t, 17> c_trunc_mask_16 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff};
    static constexpr std::array<uint32_t, 33> c_trunc_mask_32 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff,
             0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff,
             0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff,
             0x3fffffff, 0x7fffffff, 0xffffffff};
    static constexpr std::array<uint64_t, 65> c_trunc_mask_64 =
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
    static uint8_t truncate(const uint8_t &v, uint_t n) {
        return v & c_trunc_mask_8[n];
    }

    static uint16_t truncate(const uint16_t &v, uint_t n) {
        return v & c_trunc_mask_16[n];
    }

    static uint32_t truncate(const uint32_t &v, uint_t n) {
        return v & c_trunc_mask_32[n];
    }

    static uint64_t truncate(const uint64_t &v, uint_t n) {
        return v & c_trunc_mask_64[n];
    }

    static int8_t truncate(const int8_t &v, uint_t n) {
        return truncate(uint8_t(v), n);
    }

    static int16_t truncate(const int16_t &v, uint_t n) {
        return truncate(uint16_t(v), n);
    }

    static int32_t truncate(const int32_t &v, uint_t n) {
        return truncate(uint32_t(v), n);
    }

    static int64_t truncate(const int64_t &v, uint_t n) {
        return truncate(uint64_t(v), n);
    }

    /**
     * make a mask which is set (1) in the range [ibegin, iend), and clear (0) elsewhere
     */
    static void make_range_mask(uint8_t &v, uint_t ibegin, uint_t iend) {
        v = c_trunc_mask_8[iend] & ~c_trunc_mask_8[ibegin];
    }

    static void make_range_mask(uint16_t &v, uint_t ibegin, uint_t iend) {
        v = c_trunc_mask_16[iend] & ~c_trunc_mask_16[ibegin];
    }

    static void make_range_mask(uint32_t &v, uint_t ibegin, uint_t iend) {
        v = c_trunc_mask_32[iend] & ~c_trunc_mask_32[ibegin];
    }

    static void make_range_mask(uint64_t &v, uint_t ibegin, uint_t iend) {
        v = c_trunc_mask_64[iend] & ~c_trunc_mask_64[ibegin];
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
    static std::string to_string(const T &v) {
        std::string tmp;
        for (uint_t i = 0ul; i < sizeof(T) * CHAR_BIT; ++i) tmp += get(v, i) ? '1' : '0';
        return tmp;
    }

}

#endif //M7_UTIL_BIT_H
