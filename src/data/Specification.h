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
    constexpr size_t itype = ~0ul;

    // complex numbers
    template<>
    constexpr size_t itype<std::complex<float>> = 0;
    template<>
    constexpr size_t itype<std::complex<double>> = 1;
    template<>
    constexpr size_t itype<std::complex<long double>> = 2;


    // real numbers
    template<>
    constexpr size_t itype<float> = 3;
    template<>
    constexpr size_t itype<double> = 4;
    template<>
    constexpr size_t itype<long double> = 5;

    // signed integers
    template<>
    constexpr size_t itype<char> = 6;
    template<>
    constexpr size_t itype<short int> = 7;
    template<>
    constexpr size_t itype<int> = 8;
    template<>
    constexpr size_t itype<long int> = 9;
    template<>
    constexpr size_t itype<long long int> = 10;

    // unsigned integers
    template<>
    constexpr size_t itype<unsigned char> = 11;
    template<>
    constexpr size_t itype<unsigned short int> = 12;
    template<>
    constexpr size_t itype<unsigned int> = 13;
    template<>
    constexpr size_t itype<unsigned long int> = 14;
    template<>
    constexpr size_t itype<unsigned long long int> = 15;

    // boolean
    template<>
    constexpr size_t itype<bool> = 16;

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
    void create(size_t n) {
        m_numeric_lengths[itype<T>] = n;
    }

    /*
     * fruitful method to be called at DataTable instantiation, whereupon the
     * Specification is copied to a const owned by the DataTable instance.
     * This allows the Specification of a datatable to be created by multiple
     * calls to the "create" methods, rather than having to pass ready-formed length
     * arrays to a DataTable constructor
     */
    Specification &commit();

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
