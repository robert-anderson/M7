//
// Created by rja on 25/01/2021.
//

#ifndef M7_COLUMN_H
#define M7_COLUMN_H

#include "ColumnSpecifier.h"
#include "src/core/nd/NdFormat.h"

struct Table;

/**
 * Describes the memory formatting of a table buffer in terms of the primitive
 * data types encoded (e.g. Numbers, Arrays of Numbers, Bitsets).
 * Identifies the starting location of a string of bytes within a table's buffer
 * given a row index. This string of bytes is a *raw view* on an *element*,
 * whose exact format within that byte string is defined in the templated subclasses
 */
struct ColumnBase {
    Table *m_table;
    ColumnData m_data;
    const std::string m_description;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_offset;

    char *begin(const size_t &irow) const;

    char *raw_ptr(const size_t &irow, const size_t &ielement) const;

    ColumnBase(Table *table, ColumnData column_data,
               size_t nelement, std::string description);

    ColumnBase(const ColumnBase &other);

    bool is_same_type_as(const ColumnBase &other) const;

    virtual std::string to_string(size_t irow) const = 0;
};

/**
 * Builds on ColumnBase by introducing the column specification, which contains typedef
 * view_t whose constructor accepts a byte (aka raw view). The get_view
 * methods return a (const or modifiable) spec_t::view_t instance by treating
 * the nelement elements of the column as a flat (1D) vector.
 *
 * Instances of Column should not directly be added to Tables. This is the purpose of
 * the NdFieldGroup (which may be composited out of many Columns), and NdField (which
 * is the commonly-used single-column derived class of NdFieldGroup.
 *
 * Multidimensionality is added in subclass.
 */
template<typename spec_t>
struct Column : ColumnBase {
    static_assert(std::is_base_of<ColumnSpecifier, spec_t>::value, "Template arg must be derived from ColumnSpecifier");
    typedef typename spec_t::view_t view_t;
    const spec_t m_spec;

    Column(Table *table, spec_t spec, size_t nelement, std::string description) :
            ColumnBase(table, static_cast<const ColumnSpecifier &>(spec).m_data, nelement, description), m_spec(spec) {}

    const ColumnSpecifier &spec() const {
        return static_cast<const ColumnSpecifier &>(m_spec);
    }

    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement)
            res += spec().element_string(raw_ptr(irow, ielement)) + " ";
        return res;
    }

    view_t get_view(const size_t &irow, const size_t &ielement) {
        return m_spec(raw_ptr(irow, ielement));
    }

    const view_t get_view(const size_t &irow, const size_t &ielement) const {
        return m_spec(raw_ptr(irow, ielement));
    }

    view_t get_view(char* ptr){
        return m_spec(ptr);
    }

    const view_t get_view(const char* ptr) const{
        return m_spec(ptr);
    }
};


/**
 * Extends Column by enabling multidimensional access. Instances should
 * not be added directly to Tables as explained above.
 */
template<typename spec_t, size_t nind>
struct NdColumn : Column<spec_t> {
    typedef typename spec_t::view_t view_t;
    const NdFormat<nind> &m_format;

    NdColumn(Table *table, spec_t spec, std::string description, const NdFormat<nind> &format) :
            Column<spec_t>(table, spec, format.nelement(), description), m_format(format) {}

    using Column<spec_t>::get_view;

    template<typename ...Args>
    view_t operator()(const size_t &irow, Args... inds) {
        return get_view(irow, m_format.flatten(inds...));
    }

    template<typename ...Args>
    const view_t operator()(const size_t &irow, Args... inds) const {
        return get_view(irow, m_format.flatten(inds...));
    }

    struct RowView {
        NdColumn &m_column;
        // pointer to beginning of column in
        char *m_ptr;

        RowView(NdColumn &column, const size_t &irow) :
                m_column(column), m_ptr(m_spec(column.begin(irow))) {}

        template<typename ...Args>
        view_t operator()(Args... inds) {
            return get_view(m_ptr+m_column.m_format.flatten(inds...)*m_column.m_data.m_element_size);
        }

        template<typename ...Args>
        const view_t operator()(Args... inds) const {
            return get_view(m_ptr+m_column.m_format.flatten(inds...)*m_column.m_data.m_element_size);
        }
    };

    RowView operator[](const size_t &irow) {
        return RowView(*this, irow);
    }

    const RowView operator[](const size_t &irow) const {
        return RowView(*this, irow);
    }
};


#endif //M7_COLUMN_H
