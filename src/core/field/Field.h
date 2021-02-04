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


#ifndef M7_FIELD_H
#define M7_FIELD_H

#include "Column.h"
#include "src/core/hash/Hashing.h"
#include "src/core/nd/NdFormat.h"
#include "src/core/parallel/MPIAssert.h"

/**
 * FieldGroups can be composited from ColumnBase subclasses, e.g. FermiBosOnv,
 * which contains a
 * NdFieldBase<FermionOnvSpecifier, nind> for the fermion ONV, and a
 * NdFieldBase<BosonOnvSpecifier, nind> for the boson ONV
 */
template<size_t nind>
struct NdFieldGroup {
    NdFormat<nind> m_format;

    /**
     * @return
     * The offset in bytes between the address of a table and this Field
     */
    size_t symbol_offset(Table* table) const {
        auto tmp = std::distance((const char*)table, (const char *)this);
        MPI_REQUIRE_ALL(tmp>0, "Invalid Table pointer");
        return tmp;
    }

    template<typename ...Args>
    NdFieldGroup(Args... shape): m_format(shape...) {}

};

/**
 * This class *should be used* to add non-composite fields to Tables, since it
 * provides a uniformity of interface with the composite subclasses of NdFieldGroup
 */
template<typename spec_t, size_t nind>
struct NdField : NdFieldGroup<nind> {
    static_assert(std::is_base_of<ColumnSpecifier, spec_t>::value, "Template arg should be derived from ColumnSpecifier");
    NdColumn<spec_t, nind> m_column;
    using NdFieldGroup<nind>::m_format;
    typedef typename spec_t::view_t view_t;
    typedef typename spec_t::cview_t cview_t;

    typedef typename NdColumn<spec_t, nind>::rview_t rview_t;
    typedef typename NdColumn<spec_t, nind>::crview_t crview_t;

    template<typename ...Args>
    NdField(Table *table, spec_t spec, std::string description, Args... shape):
            NdFieldGroup<nind>(shape...), m_column(table, spec, description, m_format) {}

    bool next(cview_t& view) const {
        auto tmp = static_cast<ColumnSpecifier::View&>(view).m_ptr + m_column.m_data.m_row_size;
        if (tmp < (char*)static_cast<const Table&>(m_column.m_table).m_bw.m_dend){
            static_cast<ColumnSpecifier::View&>(view).m_ptr = tmp;
            return true;
        }
        return false;
    }

    view_t get_view(const size_t &irow, const size_t &ielement){
        return m_column.get_view(irow, ielement);
    }

    cview_t get_view(const size_t &irow, const size_t &ielement) const {
        return m_column.get_view(irow, ielement);
    }

    template<typename ...Args>
    view_t operator()(const size_t &irow, Args... inds) {
        return m_column(irow, inds...);
    }

    template<typename ...Args>
    cview_t operator()(const size_t &irow, Args... inds) const {
        return m_column(irow, inds...);
    }

    rview_t operator[](const size_t &irow) {
        return m_column[irow];
    }

    crview_t operator[](const size_t &irow) const {
        return m_column[irow];
    }

    struct hash_fn {
        defs::hash_t operator()(const view_t &view) const {
            return spec_t::hash(view);
        }
    };
};


#endif //M7_FIELD_H
