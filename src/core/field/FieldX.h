//
// Created by rja on 06/02/2021.
//

#ifndef M7_FIELDX_H
#define M7_FIELDX_H

#include <map>
#include "src/core/nd/NdSequence.h"
#include "src/core/hash/Hashing.h"


struct RowX;


struct FieldBaseX {
    RowX *m_row;
    const bool m_is_key;
    const size_t m_element_size;
    const size_t m_nelement;
    const size_t m_size;
    const std::type_info &m_type_info;
    const size_t m_offset;
    mutable size_t m_element_offset = 0;
    std::map<std::string, std::string> m_details;

    FieldBaseX(RowX *row, bool is_key, size_t nelement, const std::type_info &type_info, size_t element_size);

    FieldBaseX(const FieldBaseX &other);

private:

    FieldBaseX &select_first() {
        m_element_offset = 0;
        return *this;
    }

    FieldBaseX &select_next() {
        m_element_offset += m_element_size;
        return *this;
    }

    template<typename ...Args>
    FieldBaseX &select(const size_t iflat) {
        m_element_offset = m_element_size * iflat;
        return *this;
    }

public:
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

    char *begin();

    const char *begin() const;

    char *raw_view();

    const char *raw_view() const;

    void zero();

    void zero_all();

    bool equals(const FieldBaseX& other) const{
        return std::memcmp(raw_view(), other.raw_view(), m_size)==0;
    }

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
struct NumberFieldX : FieldBaseX {
    NumberFieldX(RowX *row, bool is_key, size_t nelement) :
            FieldBaseX(row, is_key, nelement, typeid(T), sizeof(T)) {}

    const T & operator()() const {
        return *(T *) raw_view();
    }

    T & operator()() {
        return *(T *) raw_view();
    }

    /*
    NumberFieldX &operator=(const T &v) {
        (T &) *this = v;
        return *this;
    }

    operator const T &() const {
        return *(T *) raw_view();
    }

    operator T &() {
        return *(T *) raw_view();
    }
     */

    std::string to_string() const override {
        return std::to_string((const T &) *this);
    }
};

template<typename T, size_t nind>
struct NumberArrayFieldX : FieldBaseX {
    NdFormat<nind> m_format;

    NumberArrayFieldX(RowX *row, bool is_key, size_t nelement, const NdFormat<nind> &format) :
            FieldBaseX(row, is_key, nelement, typeid(NumberArrayFieldX<T, nind>), sizeof(T) * format.nelement()),
            m_format(format) {}

    NumberArrayFieldX &operator=(const std::vector<T> &v) {
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
        std::string tmp = "[";
        const auto nelement = m_format.nelement();
        for (size_t i = 0ul; i < nelement; ++i)
            tmp += std::to_string(this->operator[](i)) + " ";
        return tmp + "]";
    }
};


struct BitsetFieldX : FieldBaseX {
    const size_t m_nbit;
    const size_t m_dsize;

public:
    BitsetFieldX(RowX *row, bool is_key, size_t nelement, size_t nbit) :
            FieldBaseX(row, is_key, nelement, typeid(BitsetFieldX),
                       integer_utils::divceil(nbit, defs::nbit_data) * defs::nbyte_data),
            m_nbit(nbit), m_dsize(m_size * defs::nbyte_data) {}

    const size_t &nbit() const {
        return m_nbit;
    }

    bool get(const size_t &ibit) const {
        ASSERT(ibit < nbit());
        return bit_utils::get(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void set(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::set(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void set(const defs::inds &setinds) {
        for (auto ibit: setinds) set(ibit);
    }

    void clr(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::clr(((defs::data_t *) raw_view())[ibit / defs::nbit_data], ibit % defs::nbit_data);
    }

    void clr(const defs::inds &setinds) {
        for (auto ibit: setinds) clr(ibit);
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

    size_t nset() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            result += bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }


    std::string to_string() const override {
        std::string res;
        res.reserve(nbit());
        for (size_t i = 0ul; i < nbit(); ++i) res += get(i) ? "1" : "0";
        return res;
    }

    BitsetFieldX &operator=(const defs::inds &ibits) {
        zero();
        for (auto ibit: ibits) set(ibit);
        return *this;
    }
};


struct FermionOnvFieldX : BitsetFieldX {
    const size_t m_nsite;

    FermionOnvFieldX(RowX *row, bool is_key, size_t nelement, size_t nsite) :
            BitsetFieldX(row, is_key, nelement, 2 * nsite), m_nsite(nsite) {}

    using BitsetFieldX::set;
    using BitsetFieldX::get;
    using BitsetFieldX::clr;
    using BitsetFieldX::operator=;

    void set(const size_t &ispin, const size_t &iorb) {
        set(ispin * m_nsite + iorb);
    }

    void set(const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) set(0, i);
        for (auto i: beta) set(1, i);
    }

    void set(const std::string &s) {
        zero();
        size_t i = 0ul;
        for (auto c: s) {
            // divider
            if (c != ',') {
                if (c != '0' && c != '1')
                    throw std::runtime_error(
                            R"(FermionOnv-defining string must contain only "0", "1", or ",")");
                if (c == '1') set(i);
                ++i;
            } else {
                if (i != m_nsite)
                    throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
            }
        }
        MPI_REQUIRE(i > m_nbit, "FermionOnv-defining string not long enough");
        MPI_REQUIRE(i < m_nbit, "FermionOnv-defining string too long");
    }

    void clr(const size_t &ispin, const size_t &iorb) {
        clr(ispin * m_nsite + iorb);
    }

    bool get(const size_t &ispin, const size_t &iorb) const {
        return get(ispin * m_nsite + iorb);
    }

    void excite(const size_t &i, const size_t &j) {
        /*
         * single excitation i->j
         */
        clr(i);
        set(j);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        /*
         * double excitation i,j->k,l
         */
        clr(i);
        clr(j);
        set(k);
        set(l);
    }

    int spin() const {
        int spin = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            work = get_dataword(idataword);
            while (work) {
                size_t ibit = idataword * defs::nbit_data + bit_utils::next_setbit(work);
                if (ibit < m_nsite) ++spin;
                else if (ibit >= 2 * m_nsite) return spin;
                else --spin;
            }
        }
        return spin;
    }

    int nalpha() const {
        // number of electrons occupying spinors in the alpha spin channel
        int nalpha = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            work = get_dataword(idataword);
            while (work) {
                size_t ibit = idataword * defs::nbit_data + bit_utils::next_setbit(work);
                if (ibit >= m_dsize) return nalpha;
                nalpha++;
            }
        }
        return nalpha;
    }
};

struct BosonOnvFieldX : NumberArrayFieldX<uint8_t, 1> {
    size_t m_nmode;
    using NumberArrayFieldX<uint8_t, 1>::operator=;
    using NumberArrayFieldX<uint8_t, 1>::operator[];
    using NumberArrayFieldX<uint8_t, 1>::operator();

    BosonOnvFieldX(RowX *row, bool is_key, size_t nelement, size_t nmode) :
            NumberArrayFieldX<uint8_t, 1>(row, is_key, nelement, {nmode}), m_nmode(nmode) {}
};

template<size_t nind>
struct NdFieldBaseX {
    NdSequence<nind> m_sequence;
    typename NdSequence<nind>::Cursor m_inds;

    template<typename ...Args>
    NdFieldBaseX(std::array<size_t, nind> shape) :
            m_sequence(shape),
            m_inds(m_sequence.cursor()) {}
};


template<typename field_t, size_t nind>
struct NdFieldX : field_t, NdFieldBaseX<nind> {
    static_assert(std::is_base_of<FieldBaseX, field_t>::value, "Template arg must be derived from FieldBase");
    static_assert(!std::is_same<FieldBaseX, field_t>::value, "Template arg must not be FieldBase directly");
    using NdFieldBaseX<nind>::m_inds;
    using FieldBaseX::begin;
    using FieldBaseX::m_size;
    using FieldBaseX::m_element_offset;
    using FieldBaseX::m_element_size;
    using field_t::operator=;

    template<typename ...Args>
    NdFieldX(RowX *row, bool is_key, std::array<size_t, nind> shape, Args... args) :
            field_t(row, is_key, NdFormat<nind>(shape).nelement(), args...),
            NdFieldBaseX<nind>(shape) {
            }

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

/*
template<size_t nind>
struct FermiBosOnvNdFieldX : NdFieldBaseX<nind> {
    using NdFieldBaseX<nind>::m_inds;
    using Hashable::m_key_fields;

    FermionOnvFieldX m_fonv;
    BosonOnvFieldX m_bonv;

    FermiBosOnvNdFieldX(RowX *row, std::array<size_t, nind> shape, size_t nsite) :
            NdFieldBaseX<nind>(shape),
            m_fonv(row, NdFormat<nind>(shape).nelement(), nsite),
            m_bonv(row, NdFormat<nind>(shape).nelement(), nsite) {
        m_key_fields = {&m_fonv, &m_bonv};
    }

    FermiBosOnvNdFieldX &select_first() {
        m_inds.to_front();
        m_fonv.m_element_offset = 0l;
        m_bonv.m_element_offset = 0l;
        return *this;
    }

    FermiBosOnvNdFieldX &select_next() {
        m_inds++;
        m_fonv.m_element_offset += m_fonv.m_element_size;
        m_bonv.m_element_offset += m_bonv.m_element_size;
        return *this;
    }

    template<typename ...Args>
    FermiBosOnvNdFieldX &select(Args... inds) {
        m_inds.to(inds...);
        m_fonv.m_element_offset += m_fonv.m_element_size * m_inds.ielement();
        m_bonv.m_element_offset += m_bonv.m_element_size * m_inds.ielement();
        return *this;
    }
};
 */

struct TableBaseX;

struct RowX {
    TableBaseX *m_table = nullptr;
    mutable defs::data_t *m_dbegin = nullptr;
    std::vector<FieldBaseX *> m_fields;
    defs::inds m_key_field_inds;
    size_t m_size;
    size_t m_dsize;
    size_t m_current_offset = 0ul;
    mutable RowX *m_last_copied = nullptr;

    void set_table(TableBaseX* table){
        m_table = table;
        select_first();
    }

    const RowX& select_first() const;
    const RowX& select_next() const;
    const RowX& select(const size_t& irow) const;

    defs::inds field_format() const {
        defs::inds tmp;
        tmp.reserve(m_fields.size());
        for (auto field: m_fields) tmp.push_back(field->m_size*(size_t)&field->m_type_info);
        return tmp;
    }

    RowX(){}
    RowX(TableBaseX* table){
        set_table(table);
    }
    RowX(const RowX &other);

    std::string to_string() const {
        std::string tmp;
        for (auto field: m_fields) tmp += " " + field->to_string_all();
        return tmp;
    }

    size_t add_field(FieldBaseX *field) {
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

        if (field->m_is_key) m_key_field_inds.push_back(m_fields.size());
        m_fields.push_back(field);
        return offset;
    }

    void copy_keys(const RowX& other){
        ASSERT(other.m_key_field_inds.size()==m_key_field_inds.size());
        for (const auto& ifield : m_key_field_inds) {
            m_fields[ifield]->copy(*other.m_fields[ifield]);
        }
    }

    defs::hash_t hash_keys() const {
        defs::hash_t tmp = 0ul;
        for (const auto& ifield: m_key_field_inds) {
            const auto field = m_fields[ifield];
            tmp^=hashing::fnv_hash(field->begin(), field->m_size);
        }
        return tmp;
    }

    bool equal_keys(const RowX& other) const {
        ASSERT(other.m_key_field_inds.size()==m_key_field_inds.size());
        for (const auto& ifield : m_key_field_inds) {
            if (!m_fields[ifield]->equals(*other.m_fields[ifield])) return false;
        }
        return true;
    }

};

namespace fieldxs {

    template<typename T, size_t nind>
    struct Numbers : NdFieldX<NumberFieldX<T>, nind> {
        Numbers(RowX *row, std::array<size_t, nind> shape, bool is_key=false) :
        NdFieldX<NumberFieldX<T>, nind>(row, is_key, shape) {}
    };

    template<typename T>
    struct Number : Numbers<T, 0> {
        Number(RowX *row, bool is_key=false) : Numbers<T, 0>(row, is_key, {}) {}
    };

    template<typename T, size_t nind, size_t nind_element>
    struct NumberArrays : NdFieldX<NumberArrayFieldX<T, nind_element>, nind> {
        NumberArrays(RowX *row, std::array<size_t, nind> shape, const NdFormat<nind_element> &element_format, bool is_key=false) :
                NdFieldX<NumberArrayFieldX<T, nind_element>, nind>(row, is_key, shape, element_format) {}
    };

    template<typename T, size_t nind_element>
    struct NumberArray : NumberArrays<T, 0, nind_element> {
        NumberArray(RowX *row, const NdFormat<nind_element> &element_format, bool is_key=false) :
                NumberArrays<T, 0, nind_element>(row, is_key, {}, element_format) {}
    };

    template<size_t nind>
    struct Bitsets : NdFieldX<BitsetFieldX, nind> {
        Bitsets(RowX *row, std::array<size_t, nind> shape, size_t nbit, bool is_key=false) :
                NdFieldX<BitsetFieldX, nind>(row, is_key, shape, nbit) {}
    };

    struct Bitset : Bitsets<0> {
        Bitset(RowX *row, size_t nbit, bool is_key=false) : Bitsets<0>(row, {}, nbit, is_key) {}
    };


    template<size_t nind>
    struct FermionOnvs : NdFieldX<FermionOnvFieldX, nind> {
        FermionOnvs(RowX *row, std::array<size_t, nind> shape, size_t nsite, bool is_key=false) :
                NdFieldX<FermionOnvFieldX, nind>(row, is_key, shape, nsite) {}
    };

    struct FermionOnv : FermionOnvs<0> {
        FermionOnv(RowX *row, size_t nsite, bool is_key=false) : FermionOnvs<0>(row, {}, nsite, is_key) {}
    };

    template<size_t nind>
    struct BosonOnvs : NumberArrays<uint8_t, nind, 1> {
        BosonOnvs(RowX *row, std::array<size_t, nind> shape, size_t nmode, bool is_key=false):
        NumberArrays<uint8_t, nind, 1>(row, shape, {nmode}, is_key){}
    };

    struct BosonOnv : BosonOnvs<0>{
        BosonOnv(RowX *row, size_t nmode, bool is_key=false):
        BosonOnvs<0>(row, {}, nmode, is_key){}
    };


//    template<size_t nind>
//    using FermiBosOnvs = FermiBosOnvNdFieldX<nind>;
//
//    struct FermiBosOnv : FermiBosOnvs<0> {
//        FermiBosOnv(RowX *row, size_t nbit) : FermiBosOnvs<0>(row, {}, nbit) {}
//    };

}


#endif //M7_FIELDX_H
