//
// Created by rja on 21/10/2020.
//

#ifndef M7_BITSETSPECIFIER_H
#define M7_BITSETSPECIFIER_H

#include "FieldSpecifier.h"
#include "src/core/util/utils.h"

struct BitsetSpecifier : FieldSpecifier {
    const size_t m_nbit;
    const size_t m_ndataword;

    struct View : FieldSpecifier::View {

        struct BitView {
            /*
             * This is a view on a view, if the BitView is to live longer than the View
             * from which it was derived, the BitView will need its own copy of that view,
             * which we store in a smart ptr.
             */
            std::unique_ptr<View> m_view;
            const size_t m_ibit;
            BitView(const View& view, const size_t &ibit);

            /*
             * This smart ptr is not trivially copyable, so we need to explicitly implement
             * the copy ctor
             */
            BitView(const BitView& bv);
            BitView& operator=(bool t);
            operator bool() const;
        };
        typedef BitView bitview_t;
        typedef const BitView const_bitview_t;

        View(const BitsetSpecifier &field, char* ptr);

        View(const View& other);

        bitview_t operator[](const size_t& ibit);
        const_bitview_t operator[](const size_t& ibit) const;

        const size_t &nbit() const;

        const size_t &ndataword() const;

        bool get(const size_t &ibit) const;

        defs::data_t get_dataword(const size_t &idataword) const;

        defs::data_t get_antidataword(const size_t &idataword) const;

        size_t nsetbit() const;

        std::string to_string() const override;

        void set(const size_t &ibit);

        void set(const defs::inds &setinds);

        void clr(const size_t &ibit);

        View& operator=(const defs::inds& ibits);

        View& operator=(const View& other);
    };

    typedef View view_t;
    typedef const View const_view_t;

    BitsetSpecifier(size_t nbit);

    std::string element_string(char *ptr) const override;

    const View operator()(char *ptr) const;
    View operator()(char *ptr);
};


#endif //M7_BITSETSPECIFIER_H
