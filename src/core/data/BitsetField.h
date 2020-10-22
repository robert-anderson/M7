//
// Created by rja on 21/10/2020.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include "NdField.h"
#include "src/core/util/utils.h"

template<size_t nind>
struct BitsetFieldX : NdFieldX<nind> {
    const size_t m_nbit;
    BitsetFieldX(TableX *table, std::array<size_t, nind> shape, size_t nbit, std::string description) :
            NdFieldX<nind>(table, shape, defs::nbyte_data*integer_utils::divceil(nbit, defs::nbit_data), description),
            m_nbit(nbit){}

    struct View : FieldX::View {
        View(const BitsetFieldX& field, char* ptr): FieldX::View(field, ptr){}

        const size_t& nbit() const {
            return static_cast<const BitsetFieldX<nind>&>(m_field).m_nbit;
        }

        bool get(const size_t& ibit) const {
            ASSERT(ibit<nbit());
            return bit_utils::get(*dptr(ibit/defs::nbit_data), ibit%defs::nbit_data);
        }
        void set(const size_t& ibit) {
            ASSERT(ibit<nbit());
            bit_utils::set(*dptr(ibit/defs::nbit_data), ibit%defs::nbit_data);
        }
        void set(const defs::inds& setinds){
            for (auto& i: setinds) set(i);
        }
        void clr(const size_t& ibit) {
            ASSERT(ibit<nbit());
            bit_utils::clr(*dptr(ibit/defs::nbit_data), ibit%defs::nbit_data);
        }
        std::string to_string() const override{
            std::string res;
            res.reserve(nbit());
            for (size_t i=0ul; i<nbit(); ++i) res+=get(i)?"1":"0";
            return res;
        }
    };

    template<typename ...Args>
    View operator()(const size_t& irow, Args...inds){
        return View(*this, NdFieldX<nind>::raw_ptr(irow, inds...));
    }

    std::string element_string(size_t irow, size_t ielement) const override {
        return View(*this, FieldX::raw_ptr(irow, ielement)).to_string();
    }

};


#endif //M7_BITSETFIELD_H
