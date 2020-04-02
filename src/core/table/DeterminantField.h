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

    void set(const defs::inds &ispinorbs);

    void clr(const size_t &ispin, const size_t &iorb);

    bool get(const size_t &ispin, const size_t &iorb) const;

    size_t nsite() const;

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
};

class DeterminantField : public BitsetField {
public:
    const size_t m_nsite;

    DeterminantField(Table *table, size_t nelement, size_t nsite) :
        BitsetField(table, nelement, nsite*2), m_nsite(nsite) {}

    DeterminantElement element(const size_t &irow, const size_t &isegment, const size_t &ielement) {
        return DeterminantElement(this, element_begin(irow, isegment, ielement));
    }

    virtual std::string to_string(size_t irow, size_t isegment, size_t ibegin, size_t iend) {
        std::string result = "";
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement) {
            result += element(irow, isegment, ielement).to_string() + " ";
        }
        return result;
    }
};


#endif //SANDBOX2_DETERMINANTFIELD_H
