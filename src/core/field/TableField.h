/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 * The inheritance developed between the class definitions in this file ultimately provide
 * the NdField template, which has the ability to multi-dimensionally access formatted views
 * on a Table's BufferWindow. The NdField's parent class, NdFieldGroup, can be extended as in
 * the case of FermionBosonOnv to provide analogous functionality for composite fields.
 */


#ifndef M7_TABLEFIELD_H
#define M7_TABLEFIELD_H

#include "FieldSpecifier.h"
#include "src/core/hash/Hashing.h"
#include "src/core/nd/NdFormat.h"

struct Table;

/**
 * Primitive field type, identifies the starting location of a
 * string of bytes within a table's buffer given a row index. This string
 * of bytes is a *raw view* on an *element*, whose exact format within that
 * byte string is defined in the templated subclasses
 */
struct TableField {
    Table *m_table;
    FieldData m_data;
    const std::string m_description;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_offset;

    char *begin(const size_t &irow) const;

    char *raw_ptr(const size_t &irow, const size_t &ielement) const;

    TableField(Table *table, FieldData field_data,
               size_t nelement, std::string description);

    TableField(const TableField& other);

    bool is_same_type_as(const TableField &other) const;

    virtual std::string to_string(size_t irow) const = 0;
};

/**
 * Builds on TableField by introducing the field specification, which contains typedefs
 * view_t and const_view_t whose constructors accept a byte (aka raw view). The get_view
 * methods return a spec_t::view_t or spec_t::const_view_t instance by treating
 * the nelement elements of the field as a flat (1D) vector.
 *
 * Instances of Field should not directly be added to Tables.
 *
 * Multidimensionality is added in subclasses.
 */
template<typename spec_t>
struct Field : TableField {
    static_assert(std::is_base_of<FieldSpecifier, spec_t>::value, "Template arg must be derived from FieldSpecifier");
    const spec_t m_spec;

    Field(Table *table, spec_t spec, size_t nelement, std::string description) :
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

    typename spec_t::view_t get_view(const size_t &irow, const size_t &ielement) {
        return m_spec(raw_ptr(irow, ielement));
    }

    typename spec_t::const_view_t get_view(const size_t &irow, const size_t &ielement) const {
        return m_spec(raw_ptr(irow, ielement));
    }
};

/**
 * Extends Field by enabling multidimensional access. Instances should not be added directly to Tables
 */
template<typename spec_t, size_t nind>
struct NdFieldBase : Field<spec_t> {
    const NdFormat<nind> &m_format;

    NdFieldBase(Table *table, spec_t spec, std::string description, const NdFormat<nind> &format) :
            Field<spec_t>(table, spec, format.nelement(), description), m_format(format) {}

    using Field<spec_t>::get_view;
    template<typename ...Args>
    typename spec_t::view_t operator()(const size_t &irow, Args... inds) {
        return get_view(irow, m_format.flatten(inds...));
    }

    template<typename ...Args>
    typename spec_t::const_view_t operator()(const size_t &irow, Args... inds) const {
        return get_view(irow, m_format.flatten(inds...));
    }
};

/**
 * FieldGroups can be composited from TableField subclasses, e.g. FermiBosOnv,
 * which contains a
 * NdFieldBase<FermionOnvSpecifier, nind> for the fermion ONV, and a
 * NdFieldBase<BosonOnvSpecifier, nind> for the boson ONV
 */
template<size_t nind>
struct NdFieldGroup {
    NdFormat<nind> m_format;

    template<typename ...Args>
    NdFieldGroup(Args... shape): m_format(shape...) {}

};

/**
 * This class *should be used* to add non-composite fields to Tables, since it
 * provides a uniformity of interface with the composite subclasses of NdFieldGroup
 */
template<typename spec_t, size_t nind>
struct NdField : NdFieldGroup<nind> {
    NdFieldBase<spec_t, nind> m_field;
    using NdFieldGroup<nind>::m_format;
    typedef typename spec_t::view_t view_t;
    typedef typename spec_t::const_view_t const_view_t;

    template<typename ...Args>
    NdField(Table *table, spec_t spec, std::string description, Args... shape):
            NdFieldGroup<nind>(shape...), m_field(table, spec, description, m_format) {}


    typename spec_t::view_t get_view(const size_t &irow, const size_t &ielement) {
        return m_field.get_view(irow, ielement);
    }

    typename spec_t::const_view_t get_view(const size_t &irow, const size_t &ielement) const {
        return m_field.get_view(irow, ielement);
    }

    template<typename ...Args>
    typename spec_t::view_t operator()(const size_t &irow, Args... inds) {
        return m_field(irow, inds...);
    }

    template<typename ...Args>
    typename spec_t::const_view_t operator()(const size_t &irow, Args... inds) const {
        return m_field(irow, inds...);
    }

    struct hash_fn {
        defs::hash_t operator()(const view_t &view) const {
            return spec_t::hash(view);
        }
    };
};


#endif //M7_TABLEFIELD_H
