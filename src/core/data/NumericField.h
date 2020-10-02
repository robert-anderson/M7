//
// Created by rja on 02/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include <climits>
#include <src/core/util/utils.h>
#include "Field.h"
#include "Table.h"

template<typename T, size_t nind>
struct NumericField : public Field<nind> {
    struct View : Field<nind>::View {
        using Field<nind>::View::m_ptr;
        View(const NumericField<T, nind>& field, const size_t& irow, const size_t& ielement):
                Field<nind>::View(field.begin(irow)+ielement*sizeof(T)){
        }
        operator T&() {
            return *((T*)m_ptr);
        }
        operator const T&() const {
            return *((T*)m_ptr);
        }
        View& operator =(const T& v){
            *((T*)m_ptr) = v;
        }
        std::string to_string() const {
            return utils::num_to_string((T)*this);
        }
        std::string to_bit_string() const {
            std::string res;
            size_t& tmp = *(size_t*)this;
            for (size_t ibit=0ul; ibit<CHAR_BIT*sizeof(T); ++ibit)
                res+=bit_utils::get(tmp, ibit) ? '1':'0';
            return res;
        }
    };

    using FieldBase::m_table;
    using FieldBase::m_offset;
    using FieldBase::back_offset;
    using FieldBase::is_same_type_as;
    void set_offsets() {
        if (!m_table->m_fields.empty()) {
            const auto &last_field = *m_table->m_fields.back();
            if (is_same_type_as(last_field)) m_offset = last_field.back_offset();
            else m_offset = integer_utils::round_up(last_field.back_offset(), sizeof(defs::data_t));
        }
        m_table->m_tight_row_size = back_offset();
        m_table->add_field(this);
    }

    using FieldBase::m_nelement;
    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=View(*this, irow, i).to_string()+" ";
        return res;
    }

    std::string to_bit_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=View(*this, irow, i).to_bit_string();
        return res;
    }

    template<typename ...Args>
    NumericField(Table* table, Args&& ...shape) :
    Field<nind>(table, sizeof(T), typeid(T), shape...){
        set_offsets();
    }

    using Field<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds){
        return View(*this, irow, m_format.flat(inds...));
    }
};

#endif //M7_NUMERICFIELD_H
