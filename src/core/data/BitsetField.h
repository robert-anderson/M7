//
// Created by rja on 02/10/2020.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include <climits>
#include <src/core/util/utils.h>
#include "Field.h"
#include "Table.h"

/*
 * forward declaration of the "buffered view" type for Bitsets
 */
template<size_t nind> class Bitset;

template<size_t nind>
struct BitsetField : public Field<nind> {
    const size_t m_nbit;

    class View : FieldBase::View {
        using FieldBase::View::m_ptr;
        View():FieldBase::View(){}
        friend class Bitset<nind>;
    public:
        View(const BitsetField<nind> *field, const size_t& irow, const size_t& ielement):
            FieldBase::View(field, irow, ielement){}

        /*
         * A bitset is byte-addressable, so overload the address-of operator
         * clangtidy will complain that this is dangerous in general, but for
         * View subtypes it makes sense, they are just wrappers around a byte string.
         */
        const char* operator& () const {return m_ptr;}
        char* operator&() {return m_ptr;}

        inline const size_t& nbit() const{
            return static_cast<const BitsetField<nind>*>(m_field)->m_nbit;
        }
        
        struct BitView {
            View& m_view;
            size_t m_ibit;
            BitView(View& view, const size_t& ibit):m_view(view), m_ibit(ibit){}
            BitView& operator =(bool v){
                if (v) m_view.set(m_ibit);
                else m_view.clr(m_ibit);
                return *this;
            }
            operator bool() {
                return m_view.get(m_ibit);
            }
        };

        BitView operator[](const size_t& ibit){
            return BitView(*this, ibit);
        }

        void set(const size_t& ibit){
            ASSERT(ibit < nbit());
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<dsize());
            bit_utils::set( dptr()[iword], ibit-iword*defs::nbit_data);
        }

        void clr(const size_t& ibit){
            ASSERT(ibit < nbit());
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<dsize());
            bit_utils::clr(dptr()[iword], ibit-iword*defs::nbit_data);
        }

        bool get(const size_t& ibit) const{
            ASSERT(ibit < nbit());
            const size_t iword = ibit/defs::nbit_data;
            ASSERT(iword<dsize());
            return bit_utils::get(dptr()[iword], ibit-iword*defs::nbit_data);
        }

        size_t nsetbit() const {
            size_t result = 0;
            for (size_t idataword = 0ul; idataword<dsize(); ++idataword){
                result+=bit_utils::nsetbit(dptr()[idataword]);
            }
            return result;
        }

        View& operator =(const View& v){
            ASSERT(v.nbit()==nbit());
            std::copy(v.m_ptr, v.m_ptr+v.m_size, m_ptr);
        }

        std::string to_string() const {
            std::string res;
            res.reserve(nbit());
            for (size_t i=0ul; i<nbit(); ++i) res.append(get(i)?"1":"0");
            return res;
        }
    };

    using FieldBase::m_nelement;
    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=View(this, irow, i).to_string()+" ";
        return res;
    }

    template<typename ...Args>
    BitsetField(Table* table, size_t nbit, Args&& ...shape) :
            Field<nind>(table, integer_utils::divceil(nbit, (size_t)CHAR_BIT), typeid(std::vector<bool>), shape...),
            m_nbit(nbit){
        FieldBase::set_offsets();
    }

    using Field<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds) {
        return View(this, irow, m_format.flat(inds...));
    }
};

template<size_t nind>
class Bitset : public BitsetField<nind>::View {
    struct InternalTable : Table {
        BitsetField<nind> field;
        template<typename ...Args>
        InternalTable(size_t nbit, Args... shape):Table(), field(this, nbit, shape...){}
    };
    BufferedTable<InternalTable> m_table;
public:
    template<typename ...Args>
    Bitset(size_t nbit, Args... args):
    BitsetField<nind>::View(),
    m_table(nbit, args...){
        m_table.expand(1);
        m_table.push_back();
        Field<nind>::View::init(&m_table.field, 0, 0);
    }
};


#endif //M7_BITSETFIELD_H
