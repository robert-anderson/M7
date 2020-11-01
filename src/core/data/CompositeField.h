//
// Created by RJA on 01/11/2020.
//

#ifndef M7_COMPOSITEFIELD_H
#define M7_COMPOSITEFIELD_H

#include "Table.h"

struct CompositeField: FieldGroup {

    TableX* m_table;
    CompositeField(TableX *table): m_table(table) {}

    size_t add_field(const NdFieldBaseX *field) override {
        FieldGroup::add_field(field);
        auto offset = m_table->add_field(field);
        return offset;
    }

    const char *raw_ptr(const size_t& icomponent, const size_t &irow, const size_t &ielement) const {
        ASSERT(icomponent<nfield())
        return m_fields[icomponent]->raw_ptr(irow, ielement);
    }

//    std::pair<const char *, size_t> raw_view(const size_t& icomponent, const size_t &irow, const size_t &ielement) const {
//        ASSERT(icomponent<nfield())
//        return m_fields[icomponent]->raw_view(irow, ielement);
//    }

    const size_t& size(const size_t& icomponent=0) const {
        return m_fields[icomponent]->m_size;
    }

    const size_t& offset(const size_t& icomponent=0) const {
        return m_fields[icomponent]->m_offset;
    }

    struct View {
        virtual std::string to_string() const = 0;
    };
};


#endif //M7_COMPOSITEFIELD_H
