//
// Created by RJA on 01/11/2020.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include "Table.h"

struct CompositeField: FieldGroup {

    TableX* m_table;
    CompositeField(TableX *table);

    size_t add_field(const NdFieldBaseX *field) override;

    const char *raw_ptr(const size_t& icomponent, const size_t &irow, const size_t &ielement) const;

    const size_t& size(const size_t& icomponent=0) const;

    const size_t& offset(const size_t& icomponent=0) const;

    struct View {
        virtual std::string to_string() const = 0;
    };
};


#endif //M7_COMPOSITEFIELD_H
