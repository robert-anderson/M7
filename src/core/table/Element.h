//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_ELEMENT_H
#define SANDBOX2_ELEMENT_H


#include <cstddef>
#include <cstring>
#include <src/defs.h>
#include "Field.h"

class Field;

class Element {

protected:
    Field *m_field;
    char *m_begin;

public:
    typedef Field Field_T;
    Element(Field* field, char* begin);
    virtual size_t hash() const;
    virtual size_t size() const;
    virtual size_t dsize() const;
    virtual std::string to_string() const;
    bool compatible_with(const Element &rhs) const;
    defs::data_t &dataword(const size_t &idataword) const;
    defs::data_t get_dataword(const size_t &idataword) const;
    defs::data_t get_dataword(const size_t &idataword, const size_t &nbit) const;
    defs::data_t get_antidataword(const size_t &idataword) const;
    defs::data_t get_antidataword(const size_t &idataword, const size_t &nbit) const;
    int cmp(const Element& rhs) const;
    Element& operator=(const Element& rhs);
    bool operator==(const Element& rhs) const;
    bool operator!=(const Element& rhs) const;
    char *begin() const;
    size_t nbit() const;
    void zero();
    virtual bool is_zero() const;
};

#endif //SANDBOX2_ELEMENT_H
