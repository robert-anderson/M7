//
// Created by RJA on 22/09/2020.
//

#ifndef M7_PERMANENT_H
#define M7_PERMANENT_H

#include "src/core/table/Table.h"
#include "PermanentField.h"

class Permanent : public PermanentElement {
    struct PermanentTable : public Table {
        PermanentField field;

        PermanentTable(size_t nmode, size_t occ_cutoff) :
                Table("Working permanent"),
                field(this, 1, nmode, occ_cutoff) {
            expand(1);
        }
    };

    PermanentTable internal_table;
public:
    using NumericArrayElement<uint8_t>::operator=;
    Permanent(size_t nmode, size_t occ_cutoff) :
    PermanentElement(nullptr, nmode, 0, 0, 0),
    internal_table(nmode, occ_cutoff) {
        m_field = &internal_table.field;
    }
};

#endif //M7_PERMANENT_H
