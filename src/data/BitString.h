//
// Created by rja on 12/03/2020.
//

#ifndef M7_BITSTRING_H
#define M7_BITSTRING_H

#include "Field.h"

template<size_t nind = 1>
struct BitStringField : Field<defs::data_t, nind + 1> {
    const size_t m_nbit;
private:
    static std::array<size_t, nind + 1> data_shape(const std::array<size_t, nind> &shape, size_t nbit) {
        std::array<size_t, nind + 1> result;
        std::copy(shape.begin(), shape.end(), result.begin());
        result[nind] = integer_utils::divceil(nbit, TableNew::nbit_per_element<defs::data_t>());
        return result;
    }

public:

    BitStringField(TableNew *table, size_t nbit, const std::array<size_t, nind> &shape) :
            Field<defs::data_t, nind + 1>(
                    table, data_shape(shape, nbit)), m_nbit(nbit) {}

    BitStringField(TableNew *table, size_t nbit) :
            BitStringField(table, nbit, std::array<size_t, 1>{1}) {}


    using Field<defs::data_t, nind + 1>::m_table;
    using Field<defs::data_t, nind + 1>::m_indexer;

    struct Element {
        BitStringField &m_field;
        const size_t m_irow;
        const size_t m_flat;

    private:
        defs::data_t *begin() const{
            return m_field.flat_get_ptr(m_irow, m_flat);
        }

    public:

        Element(BitStringField<nind> &field, const size_t &irow, const size_t &flat) :
                m_field(field), m_irow(irow), m_flat(flat) {}

        void set(const size_t &ibit) {
            bit_utils::set(*(begin() + ibit / TableNew::nbit_per_element<defs::data_t>()),
                           ibit % TableNew::nbit_per_element<defs::data_t>());
        }

        void clr(const size_t &ibit) {
            bit_utils::clr(*(begin() + ibit / TableNew::nbit_per_element<defs::data_t>()),
                           ibit % TableNew::nbit_per_element<defs::data_t>());
        }

        bool get(const size_t &ibit) const {
            return bit_utils::get(*(begin() + ibit / TableNew::nbit_per_element<defs::data_t>()),
                                  ibit % TableNew::nbit_per_element<defs::data_t>());
        }


        void operator=(const Element &rhs) {
            std::memset(rhs.begin(), begin(), m_indexer.shape().back());
        }

        std::string to_string() const{
            std::string result = "";
            for (size_t ibit=0ul; ibit<m_field.m_nbit; ++ibit) result+= get(ibit)?"1":"0";
            return result;
        }

        void print() const{
            std::cout << to_string() << std::endl;
        }
    };

    Element operator()(const size_t &irow, const size_t &iflat=0){
        return Element(*this, irow, iflat);
    }
    Element operator()(const defs::pair &pair, const size_t &iflat=0){
        return Element(*this, m_table->pair_to_irow(pair), iflat);
    }
    Element operator()(const size_t &irow, const std::array<size_t, nind> &inds){
        return Element(*this, irow, m_indexer.get(inds));
    }
    Element operator()(const defs::pair &pair, const std::array<size_t, nind> &inds){
        return Element(*this, m_table->pair_to_irow(pair), m_indexer.get(inds));
    }

    std::string to_string(size_t irow, size_t flat){
        return (*this)(irow, flat).to_string();
    }

};


#endif //M7_BITSTRING_H
