//
// Created by rja on 21/10/2020.
//

#ifndef M7_DETERMINANTFIELD_H
#define M7_DETERMINANTFIELD_H

#include "BitsetField.h"

#if 0

template <size_t nind>
struct DeterminantFieldX : BitsetFieldX<nind>{
    const size_t m_nsite;
    DeterminantFieldX(TableX *table, std::array<size_t, nind> shape, size_t nsite, std::string description) :
            BitsetFieldX<nind>(table, shape, 2*nsite, description), m_nsite(nsite){}

    struct View : BitsetFieldX<nind>::View {
        View(const DeterminantFieldX& field, const size_t& irow, const size_t& iflat):
        BitsetFieldX<nind>::View(field, irow, iflat){}

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
        return View(*this, irow, BitsetFieldX<nind>::m_format.flatten(inds...));
    }

    std::string element_string(size_t irow, size_t ielement) const override {
        return View(*this, irow, ielement).to_string();
    }

    std::map<std::string, std::string> details() const override {
        auto map = BitsetFieldX<nind>::details();
        map["field type"] = "Determinant";
        map["number of sites"] = std::to_string(m_nsite);
        return map;
    }

};



#endif //M7_DETERMINANTFIELD_H
#endif //M7_DETERMINANTFIELD_H
