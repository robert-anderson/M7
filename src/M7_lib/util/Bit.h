//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_BIT_H
#define M7_UTIL_BIT_H

#if ENABLE_ARM_NEON
#include <arm_neon.h>
#endif
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


    static constexpr uint8_t trailz_c_table[256] = {
        8,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,  4,  0,  1,  0,  2,  0,  1,  0,
        3,  0,  1,  0,  2,  0,  1,  0,  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
        4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,  6,  0,  1,  0,  2,  0,  1,  0,
        3,  0,  1,  0,  2,  0,  1,  0,  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
        5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,  4,  0,  1,  0,  2,  0,  1,  0,
        3,  0,  1,  0,  2,  0,  1,  0,  7,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
        4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,  5,  0,  1,  0,  2,  0,  1,  0,
        3,  0,  1,  0,  2,  0,  1,  0,  4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
        6,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,  4,  0,  1,  0,  2,  0,  1,  0,
        3,  0,  1,  0,  2,  0,  1,  0,  5,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0,
        4,  0,  1,  0,  2,  0,  1,  0,  3,  0,  1,  0,  2,  0,  1,  0
    };

    static constexpr uint8_t popcnt_c_table[256] = {
         0,  1,  1,  2,  1,  2,  2,  3,  1,  2,  2,  3,  2,  3,  3,  4,  1,  2,  2,  3,  2,  3,  3,  4,
         2,  3,  3,  4,  3,  4,  4,  5,  1,  2,  2,  3,  2,  3,  3,  4,  2,  3,  3,  4,  3,  4,  4,  5,
         2,  3,  3,  4,  3,  4,  4,  5,  3,  4,  4,  5,  4,  5,  5,  6,  1,  2,  2,  3,  2,  3,  3,  4,
         2,  3,  3,  4,  3,  4,  4,  5,  2,  3,  3,  4,  3,  4,  4,  5,  3,  4,  4,  5,  4,  5,  5,  6,
         2,  3,  3,  4,  3,  4,  4,  5,  3,  4,  4,  5,  4,  5,  5,  6,  3,  4,  4,  5,  4,  5,  5,  6,
         4,  5,  5,  6,  5,  6,  6,  7,  1,  2,  2,  3,  2,  3,  3,  4,  2,  3,  3,  4,  3,  4,  4,  5,
         2,  3,  3,  4,  3,  4,  4,  5,  3,  4,  4,  5,  4,  5,  5,  6,  2,  3,  3,  4,  3,  4,  4,  5,
         3,  4,  4,  5,  4,  5,  5,  6,  3,  4,  4,  5,  4,  5,  5,  6,  4,  5,  5,  6,  5,  6,  6,  7,
         2,  3,  3,  4,  3,  4,  4,  5,  3,  4,  4,  5,  4,  5,  5,  6,  3,  4,  4,  5,  4,  5,  5,  6,
         4,  5,  5,  6,  5,  6,  6,  7,  3,  4,  4,  5,  4,  5,  5,  6,  4,  5,  5,  6,  5,  6,  6,  7,
         4,  5,  5,  6,  5,  6,  6,  7,  5,  6,  6,  7,  6,  7,  7,  8,
    };

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
    static void clr(T &x, uint_t i) {
        x &= ~(T(1) << T(i));
    }
    /**
     * set the bit in x at position i
     */
    template<typename T>
    static void set(T &x, uint_t i) {
        x |= (T(1) << T(i));
    }
    /**
     * return true if the bit in x at position i is set, else false
     */
    template<typename T>
    static bool get(const T &x, uint_t i) {
        return x & (T(1) << T(i));
    }

    template<typename T>
    static void put_byte(T &x, uint_t ibyte, char c) {
        T mask;
        // first clear the byte in x
        mask = 0xff;
        mask <<= ibyte * CHAR_BIT;
        x &= ~mask;
        // then insert the char
        mask = c;
        mask <<= ibyte * CHAR_BIT;
        x |= mask;
    }

    template<typename T>
    static void clr_byte(T &x, uint_t ibyte) {
        put_byte(x, ibyte, 0);
    }

    template<typename T>
    static char get_byte(const T &x, uint_t ibyte) {
        return (x >> ibyte * CHAR_BIT) & 0xff;
    }

    /*
     * "count trailing zeros" implementations for 32-bit and 64-bit unsigned integers
     * _tzcnt: inline the x86 instruction TZCNT
     * _clz: inline the ARM instruction CLZ
     * _c: carry out the operation in software (slow), only used when neither of the above are supported
     */

    static uint_t trailz64_tzcnt(const uint64_t &n) {
#ifdef ENABLE_TZCNT
        uint_t res;
        asm("tzcntq %1, %0;": "=r" (res): "r" (n));
        return res;
#else
        (void) n;
        return ~0ul;
#endif
    }

    static uint_t trailz32_tzcnt(const uint32_t &n) {
#ifdef ENABLE_TZCNT
        if (!n) return 32;
        uint32_t res;
        asm("tzcnt %1, %0;": "=r" (res): "r" (n));
        return res;
#else
        (void) n;
        return ~0ul;
#endif
    }

    static uint_t trailz64_clz(const uint64_t &n) {
#ifdef ENABLE_CLZ
        if (!n) return 64;
        uint_t result;
        asm("clz %x0, %x1": "=r" (result): "r" (n));
        return 63 - result;
#else
        (void) n;
        return ~0ul;
#endif
    }

    static uint_t trailz32_clz(const uint32_t &n) {
#ifdef ENABLE_CLZ
        if (!n) return 32;
        uint32_t result;
        asm("clz %w0, %w1": "=r" (result): "r" (n));
        return 31 - result;
#else
        (void) n;
        return ~0ul;
#endif
    }

    static uint_t trailz64_c(const uint64_t &n) {
        if (!n) return 64;
        uint8_t index;
        for (uint_t shift = 0ul; shift != 64; shift += 8) {
            index = (n >> shift) & 0xff;
            if (index) return bit::trailz_c_table[index] + shift;
        }
        return 64;
    }

    static uint_t trailz32_c(const uint32_t &n) {
        if (!n) return 32;
        uint8_t index;
        for (uint_t shift = 0ul; shift != 32; shift += 8) {
            index = (n >> shift) & 0xff;
            if (index) return bit::trailz_c_table[index] + shift;
        }
        return 32;
    }

    static uint_t trailz64(const uint64_t &n) {
#if defined(ENABLE_TZCNT)
        return trailz64_tzcnt(n);
#elif defined(ENABLE_CLZ)
        return trailz64_clz(n);
#else
        // resort to software implementation
        return trailz64_c(n);
#endif
    }

    static uint_t trailz32(const uint32_t &n) {
#if defined(ENABLE_TZCNT)
        return trailz32_tzcnt(n);
#elif defined(ENABLE_CLZ)
        return trailz32_clz(n);
#else
        // resort to software implementation
        return trailz32_c(n);
#endif
    }


    /*
     * "count number of set bits" implementations for 32-bit and 64-bit unsigned integers
     * _popcnt: inline the x86 instructions POPCNT (32-bit) and POPCNTQ (64-bit)
     * _arm_neon: use the ARM Neon vector suite
     * _c: carry out the operation in software (slow), only used when neither of the above are supported
     */

    static uint_t nsetbit64_popcnt(const uint64_t &n) {
        uint_t res;
        asm("popcntq %1, %0;": "=r" (res): "r" (n));
        return res;
    }

    static uint_t nsetbit32_popcnt(const uint32_t &n) {
        uint32_t res;
        asm("popcnt %1, %0;": "=r" (res): "r" (n));
        return res;
    }

    static uint_t nsetbit64_arm_neon(const uint64_t &n) {
#ifdef ENABLE_ARM_NEON
        uint8x8_t vec = vld1_u8(reinterpret_cast<const unsigned char*>(&n));
        uint64x1_t popcount = vpaddl_u32(vpaddl_u16(vpaddl_u8(vcnt_u8(vec))));
        return vget_lane_u64(popcount, 0);
#else
        (void) n;
        return ~0ul;
#endif
    }

    static uint_t nsetbit32_arm_neon(const uint32_t &n) {
        return nsetbit64_arm_neon(n);
    }

    static uint_t nsetbit64_c(const uint64_t &n) {
        if (!n) return 0;
        uint_t tot = 0ul;
        for (uint_t shift = 0ul; shift != 64; shift += 8) tot += popcnt_c_table[(n >> shift) & 0xff];
        return tot;
    }

    static uint_t nsetbit32_c(const uint32_t &n) {
        if (!n) return 0;
        uint_t tot = 0ul;
        for (uint_t shift = 0ul; shift != 32; shift += 8) tot += popcnt_c_table[(n >> shift) & 0xff];
        return tot;
    }

    static uint_t nsetbit64(const uint64_t &n) {
#if defined(ENABLE_POPCNT)
        return nsetbit64_popcnt(n);
#elif defined(ENABLE_ARM_NEON)
        return nsetbit64_arm_neon(n);
#else
        // resort to software implementation
        return nsetbit64_c(n);
#endif
    }

    static uint_t nsetbit32(const uint32_t &n) {
#if defined(ENABLE_POPCNT)
        return nsetbit32_popcnt(n);
#elif defined(ENABLE_ARM_NEON)
        return nsetbit32_arm_neon(n);
#else
        // resort to software implementation
        return nsetbit32_c(n);
#endif
    }

    static uint_t next_setbit(uint64_t &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        const auto result = trailz64(work);
        bit::clr(work, result);
        return result;
    }

    static uint_t next_setbit(uint32_t &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        const auto result = trailz32(work);
        bit::clr(work, result);
        return result;
    }

    template<typename T>
    static uint_t next_setbyte(T &work) {
        static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value, "invalid type for bit operations");
        const uint_t ibyte = next_setbit(work) / CHAR_BIT;
        clr_byte(work, ibyte);
        return ibyte;
    }

    static uint_t nsetbit(const uint64_t &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return nsetbit64(work);
    }

    static uint_t nsetbit(const uint32_t &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        return nsetbit32(work);
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
    static uint_t nsetbit_before(const T &v, uint_t n) {
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
