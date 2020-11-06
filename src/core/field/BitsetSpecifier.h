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

    BitsetSpecifier(size_t nbit);
    struct params_t{size_t nbit;};
    BitsetSpecifier(params_t p): BitsetSpecifier(p.nbit){}

    struct View : FieldSpecifier::View {
        View(const BitsetSpecifier &field, char* ptr);

        struct BitView {
            std::unique_ptr<View> m_view;
            const size_t m_ibit;
            BitView(const View& view, const size_t &ibit);
            BitView& operator=(bool t);
            operator bool();
        };
        BitView operator[](const size_t& ibit);

        View& operator=(const defs::inds& ibits){
            zero();
            for (auto& ibit : ibits) set(ibit);
            return *this;
        }

        const size_t &nbit() const;

        const size_t &ndataword() const;

        bool get(const size_t &ibit) const;

        void set(const size_t &ibit);

        void set(const defs::inds &setinds);

        void clr(const size_t &ibit);

        defs::data_t get_dataword(const size_t &idataword) const;

        defs::data_t get_antidataword(const size_t &idataword) const;

        size_t nsetbit() const;

        std::string to_string() const override;
    };

    std::string element_string(char *ptr) const override;

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const;
};


#endif //M7_BITSETSPECIFIER_H
