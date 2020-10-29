//
// Created by rja on 21/10/2020.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include "NdField.h"
#include "src/core/util/utils.h"

struct BitsetFieldX : FieldBaseX {
    const size_t m_nbit;

    BitsetFieldX(size_t nbit):
    FieldBaseX(defs::nbyte_data * integer_utils::divceil(nbit, defs::nbit_data),
               typeid(BitsetFieldX)), m_nbit(nbit){
        m_details["type"] = "Bitset";
        m_details["number of bits"] = std::to_string(m_nbit);
    }

    struct View : FieldBaseX::View {
        View(const BitsetFieldX &field, char* ptr) : FieldBaseX::View(field, ptr){}

        struct BitView {
            std::unique_ptr<View> m_view;
            const size_t m_ibit;
            BitView(const View& view, const size_t &ibit):
                    m_view(std::unique_ptr<View>(new View(view))), m_ibit(ibit){}
            BitView& operator=(bool t){
                if (t) m_view->set(m_ibit);
                else m_view->clr(m_ibit);
                return *this;
            }
            operator bool() {
                return m_view->get(m_ibit);
            }
        };
        BitView operator[](const size_t& ibit){
            return BitView(*this, ibit);
        }

        const size_t &nbit() const {
            return static_cast<const BitsetFieldX &>(m_field).m_nbit;
        }

        bool get(const size_t &ibit) const {
            ASSERT(ibit < nbit());
            return bit_utils::get(*dptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
        }

        void set(const size_t &ibit) {
            ASSERT(ibit < nbit());
            bit_utils::set(*dptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
        }

        void set(const defs::inds &setinds) {
            for (auto &i: setinds) set(i);
        }

        void clr(const size_t &ibit) {
            ASSERT(ibit < nbit());
            bit_utils::clr(*dptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
        }

        std::string to_string() const override {
            std::string res;
            res.reserve(nbit());
            for (size_t i = 0ul; i < nbit(); ++i) res += get(i) ? "1" : "0";
            return res;
        }
    };

    std::string element_string(char *ptr) const override {
        return View(*this, ptr).to_string();
    }

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const {
        return View(*this, ptr);
    }
};


#endif //M7_BITSETFIELD_H
