//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumberField.h"
#include "MultiField.h"
#include "FermionOnvField.h"
#include "BosonOnvField.h"

namespace fields {
    template<typename T, size_t nind>
    using Numbers = NdNumberField<T, nind>;

    template<typename T>
    using Number = NumberField<T>;

    template<size_t nind>
    using Flags = BitsetField<uint8_t, nind>;
    using Flag = BitField<uint8_t>;

    using FermionOnv = FermionOnvField;

    using BosonOnv = BosonOnvField;

    struct FermiBosOnv : MultiField<FermionOnv, BosonOnv> {
        const std::string m_name;
        FermionOnv &m_frm;
        BosonOnv &m_bos;

        FermiBosOnv(Row *row, size_t nsite, std::string name = "") :
                MultiField<FermionOnv, BosonOnv>(row,
                                                 {nullptr, nsite, name.empty() ? "" : name + " (fermion)"},
                                                 {nullptr, nsite, name.empty() ? "" : name + " (boson)"}),
                m_name(name), m_frm(get<0>()), m_bos(get<1>()) {
        }

        FermiBosOnv(const FermiBosOnv &other) : FermiBosOnv(other.m_frm.row_of_copy(), other.m_frm.m_nsite,
                                                            other.m_name) {}

        FermiBosOnv &operator=(const FermiBosOnv &other) {
            m_frm = other.m_frm;
            m_bos = other.m_bos;
            return *this;
        }

        FermiBosOnv &operator=(const std::pair<defs::inds, defs::inds> &inds) {
            m_frm = inds.first;
            m_bos = inds.second;
            return *this;
        }
    };

    template<bool enable_bosons = defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    struct FermionMevInds : MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> {
        typedef MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> base_t;
        const std::string m_name;
        Numbers<defs::mev_ind_t, 1> &m_ann;
        Numbers<defs::mev_ind_t, 1> &m_cre;

        FermionMevInds(Row *row, size_t nann, size_t ncre, std::string name = "") :
                base_t(row, {nullptr, {nann}, name.empty() ? "" : name + " (annihilation)"},
                       {nullptr, {ncre}, name.empty() ? "" : name + " (creation)"}),
                m_name(name), m_ann(get<0>()), m_cre(get<1>()) {
        }

        FermionMevInds(const FermionMevInds &other) :
                FermionMevInds(other.m_ann.row_of_copy(),
                               other.m_ann.m_format.nelement(),
                               other.m_cre.m_format.nelement(), other.m_name) {}

        FermionMevInds &operator=(const FermionMevInds &other) {
            m_ann = other.m_ann;
            m_cre = other.m_cre;
            return *this;
        }

        FermionMevInds &operator=(const std::pair<defs::inds, defs::inds> &inds) {
            for (size_t i = 0ul; i < inds.first.size(); ++i) m_ann[i] = inds.first[i];
            for (size_t i = 0ul; i < inds.second.size(); ++i) m_cre[i] = inds.second[i];
            return *this;
        }
    };

    /**
     * Fields are defined as symbols, i.e. members within a Row-derived class.
     * They also retain a pointer to the Row with which they are associated, and the Row in turn keeps a vector of
     * corresponding pointers to FieldBase objects. In some situations, we may want to specify a field within a row,
     * and them use that field but in the context of another row. In order to do this, the offset in memory between
     * the specified row and field is computed, then that offset is applied to the pointer to the other row, resulting
     * in a reference to the corresponding field in the other row.
     * This function is subject to run time checks ensuring safety of the casting operations.
     *
     * @tparam row_t
     *  Row-derived class containing the field_t symbol
     * @tparam field_t
     *  FieldBase-derived class defined as a member of row_t
     * @param target
     *  Row object for which the field reference is desired
     * @param source
     *  Row object for which the field reference is specified
     * @param field
     *  FieldBase-derived object which is the referenced symbol within the source row
     * @return
     *  reference to the same field within "target" as is represented by "field" within "source"
     */
    template<typename row_t, typename field_t>
    static field_t& identify(row_t& target, row_t& source, field_t& field){
        static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
        static_assert(std::is_base_of<FieldBase, field_t>::value, "Template arg must be derived from FieldBase");
        MPI_REQUIRE(static_cast<FieldBase&>(field).belongs_to_row(source), "field arg must belong to source arg");
        auto row_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&field);
        long byte_offset = field_ptr-row_ptr;
        MPI_ASSERT(byte_offset>0, "field pointer is not positively offset from row!");
        return *reinterpret_cast<field_t*>(row_ptr+byte_offset);
    }

    template<typename row_t, typename ...Args>
    static MultiField<Args...>& identify(const row_t& target, row_t& source, MultiField<Args...>& multifield){
        static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
        MPI_REQUIRE(multifield.belongs_to_row(source), "multifield arg must belong to source arg");
        auto row_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&multifield);
        long byte_offset = field_ptr-row_ptr;
        MPI_ASSERT(byte_offset>0, "field pointer is not positively offset from row!");
        return *reinterpret_cast<MultiField<Args...>*>(row_ptr+byte_offset);
    }

}

#endif //M7_FIELDS_H