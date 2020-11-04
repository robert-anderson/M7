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

        const DeterminantSpecifier& field() const {
            return static_cast<const DeterminantSpecifier&>(m_field);
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

//        class DatawordEnumerator : public Enumerator<defs::data_t> {
//        protected:
//            size_t m_idataword = 0ul;
//            const View &m_view;
//
//            virtual size_t get_dataword(const size_t &idataword);
//
//            virtual size_t get_dataword(const size_t &idataword, const size_t &nbit);
//
//            const size_t& ndataword() const;
//
//            const size_t& nbit() const;
//
//            bool next_element(defs::data_t &result) override;
//
//        public:
//            explicit DatawordEnumerator(const View &view);
//        };
//
//        class AntiDatawordEnumerator : public Enumerator<defs::data_t> {
//        protected:
//            size_t m_idataword = 0ul;
//            const View &m_view;
//
//            virtual size_t get_dataword(const size_t &idataword) {
//                return m_view.get_antidataword(idataword);
//            }
//
//            virtual size_t get_dataword(const size_t &idataword, const size_t &nbit) {
//                return m_view.get_antidataword(idataword, nbit);
//            }
//
//            const size_t& ndataword() const {
//                return m_view.field().m_ndataword;
//            }
//
//            const size_t& nbit() const {
//                return m_view.field().m_nbit;
//            }
//
//            bool next_element(defs::data_t &result) override {
//                if (m_idataword == ndataword()) return false;
//                if (m_idataword + 1 == ndataword())
//                    result = get_dataword(m_idataword, nbit() - defs::nbit_data * m_idataword);
//                else result = get_dataword(m_idataword);
//                m_idataword++;
//                return true;
//            }
//
//        public:
//            explicit AntiDatawordEnumerator(const View &view) :
//                    Enumerator<defs::data_t>(), m_view(view) {}
//        };

        int spin() const;

        int nalpha() const;

    };


    std::string element_string(char *ptr) const override;

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const;
};


#if 0

template <size_t nind>
struct DeterminantSpecifier : BitsetSpecifier<nind>{
    const size_t m_nsite;
    DeterminantSpecifier(TableX *table, std::array<size_t, nind> shape, size_t nsite, std::string description) :
            BitsetSpecifier<nind>(table, shape, 2*nsite, description), m_nsite(nsite){}

    struct View : BitsetSpecifier<nind>::View {
        View(const DeterminantSpecifier& field, const size_t& irow, const size_t& iflat):
        BitsetSpecifier<nind>::View(field, irow, iflat){}

        using BitsetSpecifier<nind>::View::m_field;
        const size_t& nsite() const {
            return static_cast<const DeterminantSpecifier<nind>&>(m_field).m_nsite;
        }

        using BitsetSpecifier<nind>::View::nbit;
        using BitsetSpecifier<nind>::View::get;
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
        return View(*this, irow, BitsetSpecifier<nind>::m_format.flatten(inds...));
    }

    std::string element_string(size_t irow, size_t ielement) const override {
        return View(*this, irow, ielement).to_string();
    }

    std::map<std::string, std::string> details() const override {
        auto map = BitsetSpecifier<nind>::details();
        map["field type"] = "Determinant";
        map["number of sites"] = std::to_string(m_nsite);
        return map;
    }

};



#endif //M7_DETERMINANTSPECIFIER_H
#endif //M7_DETERMINANTSPECIFIER_H
