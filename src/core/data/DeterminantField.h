//
// Created by rja on 21/10/2020.
//

#ifndef M7_DETERMINANTFIELD_H
#define M7_DETERMINANTFIELD_H

#include "BitsetField.h"

template <size_t nind>
struct DeterminantFieldX : BitsetFieldX<nind>{
    const size_t m_nsite;
    DeterminantFieldX(TableX *table, std::array<size_t, nind> shape, size_t nsite, std::string description) :
            BitsetFieldX<nind>(table, shape, 2*nsite, description), m_nsite(nsite){}

    struct View : BitsetFieldX<nind>::View {
        View(const DeterminantFieldX& field, char* ptr):
        BitsetFieldX<nind>::View(field, ptr){}

        using BitsetFieldX<nind>::View::m_field;
        const size_t& nsite() const {
            return static_cast<const DeterminantFieldX<nind>&>(m_field).m_nsite;
        }

        using BitsetFieldX<nind>::View::nbit;
        using BitsetFieldX<nind>::View::get;
        std::string to_string() const override {
            std::string res;
            res.reserve(nbit()+3);
            res+="(";
            for (size_t i=0ul; i<nbit(); ++i) {
                if (i==nsite()) res+=",";
                res+=get(i)?"1":"0";// spin channel delimiter
            }
            res+=")";
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



#endif //M7_DETERMINANTFIELD_H
