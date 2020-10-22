//
// Created by rja on 12/10/2020.
//

#ifndef M7_FLAGFIELD_H
#define M7_FLAGFIELD_H

#include "BitsetField.h"

template <size_t nind>
class FlagField : Field_NEW<nind>{

    std::string element_to_string(size_t irow, size_t ielement) const override {
        return std::__cxx11::string();
    }

    using Field_NEW<nind>::m_format;
    using FieldBase::m_bit_offset;
    std::map<std::string, std::string> details() const override {
        std::map<std::string, std::string> map;
        map["rank"] = std::to_string(nind);
        if (nind) map["multidimensional shape"] = utils::to_string(m_format.shape());
        map["offset (bytes)"] = std::to_string(m_offset);
        map["bit offset"] = std::to_string(m_bit_offset);
        map["number of elements"] = std::to_string(m_nelement);
        return map;
    }

public:
    using FieldBase::m_offset;
    using FieldBase::m_nelement;
    using FieldBase::back_offset;
    template<typename ...Args>
    FlagField(Table_NEW* table, std::string description, Args... shape):
    Field_NEW<nind>(table, 
                    integer_utils::divceil(NdFormat<nind>(shape...).nelement(), (size_t)CHAR_BIT),
                    typeid(bool), description, shape...){
        if (!table->m_fields.empty()) {
            const auto &last_field = *table->m_fields.back();
            if (FieldBase::is_same_type_as(last_field)){
                auto last_nbit = last_field.m_nelement+last_field.m_bit_offset;
                m_offset = last_field.m_offset+last_nbit/CHAR_BIT;
                m_bit_offset = last_nbit%CHAR_BIT;
            }
            else m_offset = integer_utils::round_up(last_field.back_offset(), sizeof(defs::data_t));
        }
        table->m_tight_row_size = back_offset();
        table->add_field(this);

    }
    //BitsetField<nind>::View
};


#endif //M7_FLAGFIELD_H
