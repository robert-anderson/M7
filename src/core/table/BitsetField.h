//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_BITSETFIELD_H
#define SANDBOX2_BITSETFIELD_H


#include "Field.h"
#include "Element.h"
#include "src/core/util/defs.h"
#include "src/core/util/utils.h"
#include "limits.h"

class BitsetField;

class BitsetElement : public Element {
public:
    typedef BitsetField Field_T;

    BitsetElement(BitsetField *field, char *begin);

    void set(const defs::pair &pair);

    virtual void set(const size_t &ibit);

    void set(const defs::inds &inds);

    void clr(const defs::pair &pair);

    virtual void clr(const size_t &ibit);

    void clr(const defs::inds &inds);

    bool get(const defs::pair &pair) const;

    bool get(const size_t &ibit) const;

    std::string to_string() override;

    size_t nsetbit() const;

    bool is_zero() const override;
};

class BitsetField : public Field {
protected:
    size_t m_nbit;
public:
    BitsetField(Table *table, size_t nelement, size_t nbit, const std::string &description = "") :
        Field(table, sizeof(defs::data_t), nelement, typeid(defs::data_t), description) {
        update_nbit(nbit);
    }

    size_t nbit() const override {
        return m_nbit;
    }

    BitsetElement operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0){
        ASSERT(ielement < m_nelement);
        return BitsetElement(this, element_begin(irow, isegment, ielement));
    }

    static defs::pair rectify_offset(const defs::pair &pair);

    std::string to_string(size_t irow, size_t isegment, size_t ielement) override;

protected:
    virtual void update_nbit(size_t nbit) {
        m_nbit = nbit;
        m_element_dsize = integer_utils::divceil(nbit, CHAR_BIT * sizeof(defs::data_t));
        m_element_size = sizeof(defs::data_t)*m_element_dsize;
    }

    void increment_nbit(size_t nbit) {
        update_nbit(m_nbit + nbit);
    }
};


#endif //SANDBOX2_BITSETFIELD_H
