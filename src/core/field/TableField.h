//
// Created by RJA on 27/10/2020.
//

#ifndef M7_TABLEFIELD_H
#define M7_TABLEFIELD_H

#include "FieldSpecifier.h"
#include "src/core/hash/Hashing.h"
#include "src/core/nd/NdFormat.h"

struct TableX;

struct TableField {
    TableX *m_table;
    FieldData m_data;
    const std::string m_description;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_offset;

    char *begin(const size_t &irow) const;

    char *raw_ptr(const size_t &irow, const size_t &ielement) const;

    TableField(TableX *table, FieldData field_data,
               size_t nelement, std::string description);

    bool is_same_type_as(const TableField &other) const;

    virtual std::string to_string(size_t irow) const = 0;
};


template<typename spec_t>
struct Field : TableField {
    static_assert(std::is_base_of<FieldSpecifier, spec_t>::value, "Template arg must be derived from FieldSpecifier");
    const spec_t m_spec;

    Field(TableX *table, spec_t spec, size_t nelement, std::string description) :
            TableField(table, static_cast<const FieldSpecifier &>(spec).m_data, nelement, description), m_spec(spec) {}

    const FieldSpecifier &spec() const {
        return static_cast<const FieldSpecifier &>(m_spec);
    }

    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement)
            res += spec().element_string(raw_ptr(irow, ielement)) + " ";
        return res;
    }

    typename spec_t::view_t operator()(const size_t &irow, const size_t &ielement) {
        return typename spec_t::view_t(m_spec, raw_ptr(irow, ielement));
    }
};

template<typename spec_t, size_t nind>
struct NdFieldBase : Field<spec_t> {
    const NdFormat<nind> &m_format;

    NdFieldBase(TableX *table, spec_t spec, std::string description, const NdFormat<nind> &format) :
            Field<spec_t>(table, spec, format.nelement(), description), m_format(format) {}

    template<typename ...Args>
    typename spec_t::view_t operator()(const size_t &irow, Args... inds) {
        return Field<spec_t>::operator()(irow, m_format.flatten(inds...));
    }
};


template<size_t nind>
struct NdFieldGroup {
    NdFormat<nind> m_format;

    template<typename ...Args>
    NdFieldGroup(Args... shape): m_format(shape...) {}
};


template<typename spec_t, size_t nind>
struct NdField : NdFieldGroup<nind> {
    NdFieldBase<spec_t, nind> m_field;
    using NdFieldGroup<nind>::m_format;
    typedef typename spec_t::view_t view_t;
    typedef typename spec_t::const_view_t const_view_t;
    typedef typename spec_t::params_t params_t;

    template<typename ...Args>
    NdField(TableX *table, spec_t spec, std::string description, Args... shape):
            NdFieldGroup<nind>(shape...), m_field(table, spec, description, m_format) {}

    template<typename ...Args>
    typename spec_t::view_t operator()(const size_t &irow, Args... inds) {
        return m_field(irow, inds...);
    }

    struct hash_fn {
        defs::hash_t operator()(const view_t &view) const {
            return spec_t::hash(view);
        }
    };
};


#endif //M7_TABLEFIELD_H
