//
// Created by rja on 12/03/2020.
//

#ifndef M7_DETERMINANT_H
#define M7_DETERMINANT_H

#include "BitString.h"

template<size_t nind = 1>
struct DeterminantField : public BitStringField<nind + 1> {
private:
    static std::array<size_t, nind + 1> data_shape(const std::array<size_t, nind> &shape) {
        std::array<size_t, nind + 1> result;
        std::copy(shape.begin(), shape.end(), result.begin());
        result[nind] = 2;
        return result;
    }

public:
    DeterminantField(TableNew *table, size_t nbit, const std::array<size_t, nind> &shape) : BitStringField<nind + 1>(
            table, nbit, data_shape(shape)) {}

    DeterminantField(TableNew *table, size_t nbit) :
            DeterminantField(table, nbit, std::array<size_t, 1>{1}) {}

    using Base = BitStringField<nind + 1>;
    using Base::m_table;
    using Base::m_indexer;

    struct Element {
        DeterminantField &m_field;
        const size_t m_irow;
        const size_t m_flat;

        std::array<typename Base::Element, 2> m_bitfields;

    public:

        Element(DeterminantField<nind> &field, const size_t &irow, const size_t &flat) :
                m_field(field), m_irow(irow), m_flat(flat),
                m_bitfields{
                        typename Base::Element(field, irow, flat),
                        typename Base::Element(field, irow, flat + m_field.m_indexer.shape().back())
                } {}


        typename Base::Element& operator[](const size_t &ispin){
            assert(ispin==0 || ispin==1);
            return m_bitfields[ispin];
        }

        void operator=(const Element &rhs) {
            m_bitfields[0] = rhs[0];
            m_bitfields[1] = rhs[1];
        }

        void set(const size_t &ibit) {
            if (ibit<m_field.m_nbit) m_bitfields[0].set(ibit);
            else m_bitfields[1].set(ibit-m_field.m_nbit);
        }

        void clr(const size_t &ibit) {
            if (ibit<m_field.m_nbit) m_bitfields[0].clr(ibit);
            else m_bitfields[1].clr(ibit-m_field.m_nbit);
        }

        bool get(const size_t &ibit) const {
            if (ibit<m_field.m_nbit) return m_bitfields[0].get(ibit);
            else return m_bitfields[1].get(ibit-m_field.m_nbit);
        }

        std::string to_string() const {
            return m_bitfields[0].to_string()+" "+m_bitfields[1].to_string();
        }

        void print() const{
            std::cout << to_string() << std::endl;
        }
    };

    Element operator()(const size_t &irow, const size_t &iflat = 0) {
        return Element(*this, irow, iflat);
    }

    Element operator()(const defs::pair &pair, const size_t &iflat = 0) {
        return Element(*this, m_table->pair_to_irow(pair), iflat);
    }

    Element operator()(const size_t &irow, const std::array<size_t, nind> &inds) {
        return Element(*this, irow, m_indexer.get(inds));
    }

    Element operator()(const defs::pair &pair, const std::array<size_t, nind> &inds) {
        return Element(*this, m_table->pair_to_irow(pair), m_indexer.get(inds));
    }

    std::string to_string(size_t irow, size_t flat) {
        return (*this)(irow, flat).to_string();
    }

};


#endif //M7_DETERMINANT_H
