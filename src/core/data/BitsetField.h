//
// Created by rja on 02/10/2020.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include <climits>
#include <src/core/util/utils.h>
#include "Field.h"
#include "Table.h"

template<size_t nind>
struct BitsetField : public Field<nind> {

    const size_t m_nbit;
    // size in data_t
    const size_t m_dsize;

    struct View : Field<nind>::View {
        const size_t& m_nbit;
        const size_t m_dsize;
        using Field<nind>::View::m_ptr;
        View(const BitsetField<nind>& field, const size_t& irow, const size_t& ielement):
                Field<nind>::View(field.begin(irow)+ielement*field.m_element_size),
                m_nbit(field.m_nbit), m_dsize(field.m_dsize){
        }

        inline defs::data_t *dptr() const {
            return (defs::data_t *) m_ptr;
        }

        struct BitView {
            View& m_view;
            size_t m_ibit;
            BitView(View& view, const size_t& ibit):m_view(view), m_ibit(ibit){}
            BitView& operator =(bool v){
                if (v) m_view.set(m_ibit);
                else m_view.clr(m_ibit);
            }
            operator bool() {
                return m_view.get(m_ibit);
            }
        };

        BitView operator[](const size_t& ibit){
            return BitView(*this, ibit);
        }

        void set(const size_t& ibit){
            ASSERT(ibit < m_nbit);
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<m_dsize);
            bit_utils::set( dptr()[iword], ibit-iword*defs::nbit_data);
        }

        void clr(const size_t& ibit){
            ASSERT(ibit < m_nbit);
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<m_dsize);
            bit_utils::clr(dptr()[iword], ibit-iword*defs::nbit_data);
        }

        bool get(const size_t& ibit) const{
            ASSERT(ibit < m_nbit);
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<m_dsize);
            return bit_utils::get(dptr()[iword], ibit-iword*defs::nbit_data);
        }

        size_t nsetbit() const {
            size_t result = 0;
            for (size_t idataword = 0ul; idataword<m_dsize; ++idataword){
                result+=bit_utils::nsetbit(dptr()[idataword]);
            }
            return result;
        }

        View& operator =(const View& v){
            ASSERT(v.m_nbit==m_nbit);
            std::copy(v.m_ptr, v.m_ptr+v.m_size, m_ptr);
        }

        std::string to_string() const {
            std::string res;
            res.reserve(m_nbit);
            for (size_t i=0ul; i<m_nbit; ++i) res.append(get(i)?"1":"0");
            return res;
        }

        std::string to_bit_string() const {
            std::string res;
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
        return res;
    }

    template<typename ...Args>
    BitsetField(Table* table, size_t nbit, Args&& ...shape) :
            Field<nind>(table, integer_utils::divceil(nbit, (size_t)CHAR_BIT), typeid(std::vector<bool>), shape...),
            m_nbit(nbit), m_dsize(integer_utils::divceil(nbit, defs::nbit_data)){
        set_offsets();
    }

    using Field<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds) {
        return View(*this, irow, m_format.flat(inds...));
    }
};



#endif //M7_BITSETFIELD_H
