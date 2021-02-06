//
// Created by rja on 06/02/2021.
//

#ifndef M7_FIELDX_H
#define M7_FIELDX_H

#include <map>
#include "src/core/nd/NdSequence.h"


struct RowX;

struct FieldBaseX {
    RowX *m_row;
    const size_t m_element_size;
    const size_t m_nelement;
    const size_t m_size;
    const std::type_info &m_type_info;
    const size_t m_offset;
    mutable size_t m_element_offset = 0;
    std::map<std::string, std::string> m_details;

    FieldBaseX(RowX *row, size_t nelement, const std::type_info &type_info, size_t element_size);

    FieldBaseX(const FieldBaseX &other);

protected:
    FieldBaseX &copy(const FieldBaseX &other) {
        std::memcpy(raw_view(), other.raw_view(), m_element_size);
        return *this;
    }

    FieldBaseX &copy_all(const FieldBaseX &other) {
        std::memcpy(begin(), other.begin(), m_size);
        return *this;
    }

    FieldBaseX &operator=(const FieldBaseX &other) {
        copy(other);
        return *this;
    }

public:

    char *begin();

    const char *begin() const;

    char *raw_view();

    const char *raw_view() const;

    void zero();

    void zero_all();

    bool is_same_type_as(const FieldBaseX &other) const {
        return m_type_info == other.m_type_info;
    }

    virtual std::string to_string() const = 0;

    std::string to_string_all() const {
        std::string tmp;
        for (m_element_offset = 0ul; m_element_offset < m_size; m_element_offset += m_element_size) {
            tmp += to_string() + " ";
        }
        return tmp;
    }

};


template<typename T>
struct NumericFieldX : FieldBaseX {
    NumericFieldX(RowX *row, size_t nelement) :
            FieldBaseX(row, nelement, typeid(T), sizeof(T)) {}

    NumericFieldX &operator=(const T &v) {
        (T &) *this = v;
        return *this;
    }

    operator const T &() const {
        return *(T *) raw_view();
    }

    operator T &() {
        return *(T *) raw_view();
    }

    std::string to_string() const override {
        return std::to_string((const T &) *this);
    }
};

template<typename T, size_t nind>
struct NumericArrayFieldX : FieldBaseX {
    NdFormat<nind> m_format;

    NumericArrayFieldX(RowX *row, size_t nelement, const NdFormat<nind> &format) :
            FieldBaseX(row, nelement, typeid(NumericFieldX<T>), sizeof(T) * format.nelement()) {}

    NumericArrayFieldX &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == m_format.nelement());
        std::memcpy(raw_view(), v.data(), m_element_size);
        return *this;
    }

    const T &operator[](const size_t &i) const {
        ASSERT(i < m_format.nelement());
        return ((T *) raw_view())[i];
    }

    T &operator[](const size_t &i) {
        ASSERT(i < m_format.nelement());
        return ((T *) raw_view())[i];
    }

    template<typename ...Args>
    T &operator()(Args... inds) {
        return ((T *) raw_view())[m_format.flatten(inds...)];
    }

    template<typename ...Args>
    const T &operator()(Args... inds) const {
        return ((T *) raw_view())[m_format.flatten(inds...)];
    }

    std::string to_string() const override {
        std::string tmp;
        for ()
    }
};


struct BitsetFieldX : FieldBaseX {
    const size_t m_nbit;
    const size_t m_dsize;

public:
    BitsetFieldX(RowX *row, size_t nelement, size_t nbit) :
            FieldBaseX(row, nelement, typeid(BitsetFieldX),
                       integer_utils::divceil(nbit, defs::nbit_data) * defs::nbyte_data),
            m_nbit(nbit), m_dsize(m_size * defs::nbyte_data) {}

    const size_t &nbit() const {
        return m_nbit;
    }

    bool getbit(const size_t &ibit) const {
        ASSERT(ibit < nbit());
        return bit_utils::get(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void setbit(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::set(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void setbit(const defs::inds &setinds) {
        for (auto ibit: setinds) setbit(ibit);
    }

    void clrbit(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::clr(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void clrbit(const defs::inds &setinds) {
        for (auto ibit: setinds) clrbit(ibit);
    }

    defs::data_t get_dataword(const size_t &idataword) const {
        auto tmp = ((defs::data_t *) raw_view())[idataword];
        if (idataword + 1 == m_dsize) {
            auto n = nbit() - defs::nbit_data * idataword;
            tmp = bit_utils::truncate(tmp, n);
        }
        return tmp;
    }

    defs::data_t get_antidataword(const size_t &idataword) const {
        auto tmp = ~(((defs::data_t *) raw_view())[idataword]);
        if (idataword + 1 == m_dsize) {
            auto n = nbit() - defs::nbit_data * idataword;
            tmp = bit_utils::truncate(tmp, n);
        }
        return tmp;
    }

    size_t nsetbit() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            result += bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }


    std::string to_string() const override {
        std::string res;
        res.reserve(nbit());
        for (size_t i = 0ul; i < nbit(); ++i) res += getbit(i) ? "1" : "0";
        return res;
    }

    BitsetFieldX &operator=(const defs::inds &ibits) {
        zero();
        for (auto ibit: ibits) setbit(ibit);
        return *this;
    }
};


template<typename field_t, size_t nind>
struct NdFieldX : field_t {
    static_assert(std::is_base_of<FieldBaseX, field_t>::value, "Template arg must be derived from FieldBase");
    static_assert(!std::is_same<FieldBaseX, field_t>::value, "Template arg must not be FieldBase directly");
    typename NdSequence<nind>::Cursor m_inds;
    using FieldBaseX::m_element_offset;
    using FieldBaseX::m_element_size;
    using field_t::operator=;

    template<typename ...Args>
    NdFieldX(RowX *row, const NdSequence<nind> &sequence, Args... args) :
            field_t(row, sequence.m_format.nelement(), args...),
            m_inds(sequence.cursor()) {}

    NdFieldX &select_first() {
        m_inds.to_front();
        m_element_offset = 0;
        return *this;
    }

    NdFieldX &select_next() {
        m_inds++;
        m_element_offset += m_element_size;
        return *this;
    }

    template<typename ...Args>
    NdFieldX &select(Args... inds) {
        m_inds.to(inds...);
        m_element_offset = m_element_size * m_inds.ielement();
        return *this;
    }

};


struct TableBaseX {
    std::vector<defs::data_t> m_data;

    TableBaseX() : m_data(100, 0ul) {}
};


struct RowX {
    TableBaseX &m_table;
    defs::data_t *m_dbegin;
    std::vector<const FieldBaseX *> m_fields;
    size_t m_size;
    size_t m_dsize;
    size_t m_current_offset = 0ul;
    mutable RowX *m_last_copied = nullptr;

    RowX(TableBaseX &table) : m_table(table), m_dbegin(m_table.m_data.data()) {

    }

    RowX(const RowX &other) : m_table(other.m_table), m_dbegin(m_table.m_data.data()) {
        other.m_last_copied = this;
    }


    std::string to_string() const {
        std::string tmp;
        for (auto field: m_fields) tmp += " " + field->to_string_all();
        return tmp;
    }

    size_t add_field(const FieldBaseX *field) {
        // returns the offset in bytes for the column being added
        auto offset = 0ul;
        if (!m_fields.empty()) {
            offset = m_fields.back()->m_offset + m_fields.back()->m_size;
            if (!m_fields.back()->is_same_type_as(*field)) {
                // go to next whole dataword
                offset = integer_utils::divceil(offset, defs::nbyte_data) * defs::nbyte_data;
            }
        }

        m_current_offset = offset + field->m_size;
        m_dsize = integer_utils::divceil(m_current_offset, defs::nbyte_data);
        m_size = m_dsize * defs::nbyte_data;

        m_fields.push_back(field);
        return offset;
    }

};


#endif //M7_FIELDX_H
