//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include "src/defs.h"
#include <omp.h>
#include <cstring>
#include <memory>
#include <assert.h>
#include "Specification.h"
#include "BitfieldNew.h"
#include "NumericView.h"
#include "src/fermion/Determinant.h"

template<typename T>
using view_t = typename std::conditional<
        std::is_same<T, BitfieldNew>::value || std::is_same<T, Determinant>::value,
        T, NumericView<T> >::type;

class Table {
public:

    template<typename T>
    struct Field {
        Table *m_table;
        const size_t m_ientry;
        const size_t m_n;

        Field(Table *table, size_t n=1) :
                m_table(table), m_ientry(table->m_spec.add<T>(n)), m_n(n) {}

        view_t<T> get(const size_t &irow, const size_t &ientry=0) {
            assert(ientry<m_n);
            assert(m_table->data());
            return m_table->view<T>(irow, m_ientry+ientry);
        }
    };

protected:
    Specification m_spec;
    size_t m_nrow = 0ul;
    std::vector<defs::data_t> m_data_internal{};
    defs::data_t *m_data = nullptr;

public:

    size_t nrow() const;

    size_t row_length() const;

    virtual size_t high_water_mark() const;

    explicit Table(defs::data_t *data_external = nullptr);

    virtual void grow(const size_t &nrow);

    virtual void grow(defs::data_t *const new_ptr, size_t nrow);

    virtual void zero(size_t irow);

    virtual void zero();

    defs::data_t *data() const {
        return m_data;
    }

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

    void print(size_t nrow, size_t irank) const;
    virtual void print(size_t irank) const;
    void print() const;

};

template<>
BitfieldNew Table::view<BitfieldNew>(const size_t &irow, const size_t &ientry) const;

template<>
Determinant Table::view<Determinant>(const size_t &irow, const size_t &ientry) const;


#endif //M7_TABLE_H
