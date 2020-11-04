//
// Created by rja on 21/10/2020.
//

#ifndef M7_DETERMINANTSPECIFIER_H
#define M7_DETERMINANTSPECIFIER_H

#include "BitsetSpecifier.h"
#include "src/core/enumerator/Enumerator.h"


struct DeterminantSpecifier : BitsetSpecifier {
    const size_t m_nsite;

    DeterminantSpecifier(const size_t &nsite);

    struct View : BitsetSpecifier::View {
        View(const DeterminantSpecifier &field, char *ptr);

        const DeterminantSpecifier& spec() const {
            return static_cast<const DeterminantSpecifier&>(m_spec);
        }

        std::string to_string() const override;

        const size_t& nsite()const;

        using BitsetSpecifier::View::set;
        using BitsetSpecifier::View::get;
        using BitsetSpecifier::View::clr;

        void set(const size_t &ispin, const size_t &iorb);

        void set(const defs::inds &ispinorbs);

        void set(const std::string& s);

        void clr(const size_t &ispin, const size_t &iorb){
            clr(ispin*nsite()+iorb);
        }

        bool get(const size_t &ispin, const size_t &iorb) const{
            return get(ispin*nsite()+iorb);
        }

        void excite(const size_t &i, const size_t &j) {
            /*
             * single excitation i->j
             */
            clr(i); set(j);
        }

        void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
            /*
             * double excitation i,j->k,l
             */
            clr(i); clr(j); set(k); set(l);
        }

        int spin() const;

        int nalpha() const;

    };


    std::string element_string(char *ptr) const override;

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const;
};


#endif //M7_DETERMINANTSPECIFIER_H