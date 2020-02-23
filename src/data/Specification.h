//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_SPECIFICATION_H
#define M7_SPECIFICATION_H

#include <cstddef>
#include <complex>
#include <vector>
#include <array>

namespace numtypes {
    template<typename T>
    constexpr size_t itype() {return ~0ul;}

    // complex numbers
    template<>
    constexpr size_t itype<std::complex<float>>() {return 0;}
    template<>
    constexpr size_t itype<std::complex<double>>() {return 1;}
    template<>
    constexpr size_t itype<std::complex<long double>>() {return 2;}


    // real numbers
    template<>
    constexpr size_t itype<float>() {return 3;}
    template<>
    constexpr size_t itype<double>() {return 4;}
    template<>
    constexpr size_t itype<long double>() {return 5;}

    // signed integers
    template<>
    constexpr size_t itype<char>() {return 6;}
    template<>
    constexpr size_t itype<short int>() {return 7;}
    template<>
    constexpr size_t itype<int>() {return 8;}
    template<>
    constexpr size_t itype<long int>() {return 9;}
    template<>
    constexpr size_t itype<long long int>() {return 10;}

    // unsigned integers
    template<>
    constexpr size_t itype<unsigned char>() {return 11;}
    template<>
    constexpr size_t itype<unsigned short int>() {return 12;}
    template<>
    constexpr size_t itype<unsigned int>() {return 13;}
    template<>
    constexpr size_t itype<unsigned long int>() {return 14;}
    template<>
    constexpr size_t itype<unsigned long long int>() {return 15;}

    // boolean
    template<>
    constexpr size_t itype<bool>() {return 16;}

    constexpr bool is_complex(size_t itype){
        if (itype<3) return true;
        else return false;
    }

    constexpr size_t ntype = 17;

    constexpr std::array<size_t, ntype> sizes = {
            sizeof(std::complex<float>),
            sizeof(std::complex<double>),
            sizeof(std::complex<long double>),
            sizeof(float),
            sizeof(double),
            sizeof(long double),
            sizeof(char),
            sizeof(short int),
            sizeof(int),
            sizeof(long int),
            sizeof(long long int),
            sizeof(unsigned char),
            sizeof(unsigned short int),
            sizeof(unsigned int),
            sizeof(unsigned long int),
            sizeof(unsigned long long int),
            sizeof(bool)
    };
}
using namespace numtypes;

class Specification {
public:
    std::array<size_t, ntype> m_numeric_lengths;
    std::array<size_t, ntype> m_numeric_datawords_used;
    std::array<size_t, ntype> m_numeric_offsets;
    std::vector<size_t> m_bitfield_lengths;
    std::vector<size_t> m_bitfield_datawords_used;
    std::vector<size_t> m_bitfield_offsets;
    size_t m_total_numeric_datawords_used;
    size_t m_total_bitfield_datawords_used;
    size_t m_total_datawords_used;

    explicit Specification(
            const std::array<size_t, ntype> &numeric_lengths = {},
            const std::vector<size_t> &bitfield_lengths = {});

    explicit Specification(const std::vector<size_t> &bitfield_lengths);

    template<typename T>
    size_t add(size_t n) {
        assert(itype<T>()!=~0ul);
        m_numeric_lengths[itype<T>()] += n;
        compile();
        return m_numeric_lengths[itype<T>()]-n;
    }

    void compile();

    bool operator==(const Specification &rhs) const;

    bool operator!=(const Specification &rhs) const;

private:
    static std::array<size_t, ntype> numeric_datawords_used(const std::array<size_t, ntype> &lengths);

    static std::array<size_t, ntype> numeric_offsets(
            const std::array<size_t, ntype> &lengths,
            const std::array<size_t, ntype> &datawords_used
    );

    static std::vector<size_t> bitfield_datawords_used(const std::vector<size_t> &lengths);

    static std::vector<size_t> bitfield_offsets(
            const std::vector<size_t> &lengths,
            const std::vector<size_t> &datawords_used
    );
};


#endif //M7_SPECIFICATION_H
