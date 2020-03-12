//
// Created by Robert John Anderson on 2020-03-09.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include <src/multidim/ArrayIndexer.h>
#include "src/utils.h"
#include "TableNew.h"
#include <cstring>
#include <algorithm>


template<typename T, size_t nind = 1>
struct Field : public TableNew::FieldBase {
    const ArrayIndexer<nind> m_indexer;
private:
    static std::array<size_t, nind> data_shape(size_t nelement) {
        std::array<size_t, nind> result;
        result.fill(nelement);
        return result;
    }

public:
    Field(TableNew *table, const std::array<size_t, nind> &shape) :
            FieldBase(table, TableNew::nbit_per_element<T>(), ArrayIndexer<nind>(shape).nelement(), typeid(T)),
            m_indexer(ArrayIndexer<nind>(shape)) {}

    Field(TableNew *table, const size_t &nelement) : Field(table, data_shape(nelement)) {}

    Field(TableNew *table) : Field(table, 1) {}

    inline T *flat_get_ptr(const size_t &irow, const size_t &flat) const {
        assert(flat < m_nelement);
        return (T *) (m_table->m_data.data() + irow * m_table->m_ndatawords_padded +
                      m_offset.first) + m_offset.second + flat;
    }

    inline T &flat_get(const size_t &irow, const size_t &flat) const {
        assert(flat < m_nelement);
        return *((T *) (m_table->m_data.data() + irow * m_table->m_ndatawords_padded +
                        m_offset.first) + m_offset.second + flat);
    }

    struct Row {
        Field<T, nind> &m_field;
        const size_t m_i;

        Row(Field<T, nind> &field, const size_t &i) : m_field(field), m_i(i) {}

        template<typename U=T>
        typename std::enable_if<!std::is_same<U, bool>::value, T &>::type
        operator()(const size_t &flat = 0) {
            return m_field.flat_get(m_i, flat);
        }

        template<typename U=T>
        typename std::enable_if<!std::is_same<U, bool>::value, T &>::type
        operator()(const std::array<size_t, nind> &inds) {
            return m_field.flat_get(m_i, m_field.m_indexer.get(inds));
        }

        void operator=(const Row &rhs) {
            assert(m_field.m_indexer == rhs.m_field.m_indexer);
            std::copy(
                    rhs.m_field.flat_get_ptr(m_i, 0),
                    rhs.m_field.flat_get_ptr(m_i, m_field.m_nelement),
                    m_field.flat_get_ptr(m_i, 0));
        }

        void operator=(const T &rhs) {
            /*
             * set all elements to the rhs
             */
            std::memset(m_field.flat_get_ptr(m_i, 0), rhs, m_field.m_nelement * sizeof(T));
        }

        int compare(const Row &rhs) const {
            assert(m_field.m_indexer == rhs.m_field.m_indexer);
            return std::memcmp(m_field.flat_get_ptr(m_i, 0), rhs.m_field.flat_get_ptr(m_i, 0), m_field.m_nelement);
        }

        bool operator==(const Row &rhs) const {
            return !compare(rhs);
        }

        bool operator==(const T &rhs) const {
            /*
             * all elements equal to the rhs?
             */
            return std::all_of(m_field.flat_get_ptr(m_i, 0), m_field.flat_get_ptr(m_i, m_field.m_nelement - 1),
                               [&rhs](const T &v) { return v == rhs; });
        }
    };

    template<typename U=T>
    typename std::enable_if<!std::is_same<U, void>::value, void>::type
    zero(const size_t &irow) {
        std::memset(flat_get(irow, 0), 0, m_nelement * sizeof(U));
    }

    template<typename U=T>
    typename std::enable_if<!std::is_same<U, void>::value, void>::type
    zero(const defs::pair &pair) {
        zero(m_table->pair_to_irow(pair));
    }

    Row row(const size_t &irow) {
        assert(irow < m_table->nrow());
        return Row(*this, irow);
    }

    Row row(const defs::pair &pair) {
        return Row(*this, m_table->pair_to_irow(pair));
    }

    T &operator()(const size_t &irow, const size_t &iflat = 0) {
        return flat_get(irow, iflat);
    }

    T &operator()(const defs::pair &pair, const size_t &iflat = 0) {
        return flat_get(m_table->pair_to_irow(pair), iflat);
    }

    T &operator()(const size_t &irow, const std::array<size_t, nind> &inds) {
        return flat_get(irow, m_indexer.get(inds));
    }

    T &operator()(const defs::pair &pair, const std::array<size_t, nind> &inds) {
        return flat_get(m_table->pair_to_irow(pair), m_indexer.get(inds));
    }

    std::string to_string(size_t irow, const std::array<size_t, nind> &inds) {
        return std::to_string(element(irow, inds));
    }

    std::string to_string(size_t irow) override {
        std::string out = "";
        for (size_t i = 0ul; i < m_nelement; ++i) out += std::to_string(flat_get(irow, i)) + " ";
        return out;
    }
};


template<size_t nind = 1>
struct Flag : public Field<bool, nind> {

    Flag(TableNew *table, const std::array<size_t, nind> &shape) : Field<bool, nind>(table, shape) {}

    using Field<bool, nind>::m_table;
    using Field<bool, nind>::m_offset;
    using Field<bool, nind>::m_indexer;
    using Field<bool, nind>::m_nelement;

    void set(const size_t &irow, const size_t &flat = 0) {
        assert(flat < m_nelement);
        assert(irow < m_table->nrow());
        auto bit_offset = m_offset.second + flat;
        auto dataword = m_table->m_data.data() + irow * m_table->m_ndatawords_padded
                        + m_offset.first + bit_offset / TableNew::nbit_per_element<defs::data_t>();
        bit_utils::set(*dataword, bit_offset % TableNew::nbit_per_element<defs::data_t>());
    }

    void set(const size_t &irow, const std::array<size_t, nind> &inds) {
        set(irow, m_indexer.get(inds));
    }

    void set(const defs::pair &pair, const std::array<size_t, nind> &inds) {
        set(m_table->pair_to_irow(pair), inds);
    }

    void clr(const size_t &irow, const size_t &flat = 0) {
        assert(flat < m_nelement);
        assert(irow < m_table->nrow());
        auto bit_offset = m_offset.second + flat;
        auto dataword = m_table->m_data.data() + irow * m_table->m_ndatawords_padded
                        + m_offset.first + bit_offset / TableNew::nbit_per_element<defs::data_t>();
        bit_utils::clr(*dataword, bit_offset % TableNew::nbit_per_element<defs::data_t>());
    }

    void clr(const size_t &irow, const std::array<size_t, nind> &inds) {
        clr(irow, m_indexer.get(inds));
    }

    void clr(const defs::pair &pair, const std::array<size_t, nind> &inds) {
        clr(m_table->pair_to_irow(pair), inds);
    }

    bool get(const size_t &irow, const size_t &flat = 0) {
        assert(flat < m_nelement);
        assert(irow < m_table->nrow());
        auto bit_offset = m_offset.second + flat;
        auto dataword = m_table->m_data.data() + irow * m_table->m_ndatawords_padded
                        + m_offset.first + bit_offset / TableNew::nbit_per_element<defs::data_t>();
        return bit_utils::get(*dataword, bit_offset % TableNew::nbit_per_element<defs::data_t>());
    }

    bool get(const size_t &irow, const std::array<size_t, nind> &inds) {
        return get(irow, m_indexer.get(inds));
    }

    bool get(const defs::pair &pair, const std::array<size_t, nind> &inds) {
        return get(m_table->pair_to_irow(pair), inds);
    }

};

#endif //M7_FIELD_H
