//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DETERMINANT_H
#define M7_DETERMINANT_H

#include "src/core/Table/Table.h"
#include "src/core/Table/DeterminantField.h"

class Determinant : public DeterminantElement {
    struct DeterminantTable : public Table {
        DeterminantField field;

        DeterminantTable(Determinant *determinant, size_t nsite) :
            Table(),
            field(this, 1, nsite) {
            expand(1);
        }
    };

    DeterminantTable internal_table;
public:
    Determinant(size_t nsite) : DeterminantElement(nullptr, nullptr), internal_table(this, nsite) {
        m_field = &internal_table.field;
        auto element = m_field->element(0);
        m_begin = element.begin();
    }

    Determinant& operator=(const DeterminantElement& rhs){
        DeterminantElement::operator=(rhs);
        return *this;
    }

    Determinant& operator=(const Determinant& rhs){
        DeterminantElement::operator=(rhs);
        return *this;
    }

};

#endif //M7_DETERMINANT_H