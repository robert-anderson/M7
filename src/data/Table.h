//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include "src/defs.h"
#include <omp.h>
#include <cstring>
#include <assert.h>
#include "Specification.h"
#include "BitfieldNew.h"
#include "NumericView.h"
#include "src/fermion/Determinant.h"

class Table {
    const Specification m_spec;
protected:
    size_t m_nrow = 0ul;
    std::vector<defs::data_t> m_data_internal{};
    defs::data_t *m_data = nullptr;

public:

    const Specification &spec() const;

    const size_t nrow() const;

    const size_t row_length() const;

    virtual const size_t high_water_mark() const;

    Table(Specification spec, size_t nrow);

    Table(Specification spec, size_t nrow, defs::data_t *data_external);

    virtual void grow(const size_t &nrow);

    virtual void grow(defs::data_t *const new_ptr, size_t nrow);

    void zero(size_t irow = ~0ul);

    defs::data_t * data() const{
        return m_data;
    }

    template<typename T>
    using view_t = typename std::conditional<
        std::is_same<T, BitfieldNew>::value || std::is_same<T, Determinant>::value,
        T, NumericView<T> >::type;

    template<typename T>
    view_t<T> view(const size_t &irow, const size_t &ientry) const {
        if (ientry >= m_spec.m_numeric_lengths[itype<T>()]) {
            return NumericView<T>(nullptr, 0);
        }
        return NumericView<T>(
            (char *) (m_data + irow * row_length() +
                      m_spec.m_numeric_offsets[itype<T>()]) + ientry * sizeof(T),
            m_spec.m_numeric_lengths[itype<T>()] - ientry);
    }

    template<typename T>
    view_t<T> view(const size_t &irow) const {
        return view<T>(irow, 0);
    }

    void print(size_t nrow = ~0ul) const;

};

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &irow, const size_t &ientry) const;

template<>
Determinant Table::view<Determinant>(const size_t &irow, const size_t &ientry) const;

#endif //M7_TABLE_H
