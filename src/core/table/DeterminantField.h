//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef SANDBOX2_DETERMINANTFIELD_H
#define SANDBOX2_DETERMINANTFIELD_H


#include "BitsetField.h"
#include "src/core/enumerator/Enumerator.h"

class DeterminantField;

class DeterminantElement : public BitsetElement {
public:
    typedef DeterminantField Field_T;

    DeterminantElement(DeterminantField *field, char *begin);

    std::string to_string() const override;

    void set(const size_t &ispin, const size_t &iorb);

    using BitsetElement::set;
    void set(const defs::inds &ispinorbs);

    using BitsetElement::clr;
    void clr(const size_t &ispin, const size_t &iorb);

    using BitsetElement::get;
    bool get(const size_t &ispin, const size_t &iorb) const;

    size_t nsite() const;

    void excite(const size_t &i, const size_t &j){
        /*
         * single excitation i<-j
         */
        BitsetElement::clr(j); BitsetElement::set(i);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l){
        /*
         * double excitation i,j<-k,l
         */
        BitsetElement::clr(k); BitsetElement::clr(l);
        BitsetElement::set(i); BitsetElement::set(j);
    }

    class DatawordEnumerator : public Enumerator<defs::data_t> {
        size_t m_idataword = 0ul;
        const DeterminantElement &m_data;

        typedef defs::data_t (DeterminantElement::*Getter)(const size_t &) const;

        typedef defs::data_t (DeterminantElement::*TruncGetter)(const size_t &, const size_t &) const;

        const Getter m_getter;
        const TruncGetter m_trunc_getter;

        Getter getter(bool anti) {
            if (anti) return &DeterminantElement::get_antidataword;
            else return &DeterminantElement::get_dataword;
        }

        TruncGetter trunc_getter(bool anti) {
            if (anti) return &DeterminantElement::get_antidataword;
            else return &DeterminantElement::get_dataword;
        }

        bool next_element(defs::data_t &result) override {
            if (m_idataword == m_data.m_field->element_dsize()) return false;
            if (m_idataword + 1 == m_data.m_field->element_dsize())
                result = (m_data.*m_trunc_getter)(m_idataword, m_data.nbit() - defs::nbit_data * m_idataword);
            else result = (m_data.*m_getter)(m_idataword);
            m_idataword++;
            return true;
        }

    public:
        DatawordEnumerator(const DeterminantElement &data, bool anti = false) :
            Enumerator<defs::data_t>(), m_data(data),
            m_getter(getter(anti)), m_trunc_getter(trunc_getter(anti)){}
    };

    int spin() const {
        int spin = 0;
        size_t idataword = ~0ul;
        defs::data_t work;
        DeterminantElement::DatawordEnumerator enumerator(*this);
        while(enumerator.next(work, idataword)){
            while (work) {
                size_t ibit = bit_utils::next_setbit(work);
                if (ibit<nsite()) ++spin;
                else --spin;
            }
        }
        return spin;
    }
};

class DeterminantField : public BitsetField {
public:
    const size_t m_nsite;

    DeterminantField(Table *table, size_t nelement, size_t nsite);

    DeterminantElement element(const size_t &irow, const size_t &isegment=0, const size_t &ielement=0);

    virtual std::string to_string(size_t irow, size_t isegment, size_t ibegin, size_t iend);
};


#endif //SANDBOX2_DETERMINANTFIELD_H
