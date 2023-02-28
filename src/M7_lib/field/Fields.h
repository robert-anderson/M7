//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumberField.h"
#include "FrmBosOnvField.h"
#include "XonvField.h"
#include "SpecMomIndsField.h"

namespace field {
    template<typename T, uint_t nind>
    using Numbers = NdNumberField<T, nind>;

    using String = StringField;

    template<typename T>
    using Number = NumberField<T>;

    template<typename T, uint_t nind>
    using NdBitset = BitsetField<T, nind>;
    template<typename T>
    using Bitset = NdBitset<T, 1ul>;

    template<uint_t nind>
    using Flags = BitsetField<uint8_t, nind>;
    using Flag = BitField<uint8_t>;

    using FrmOnv = FrmOnvField;

    using BosOnv = BosOnvField;

    using FrmBosOnv = FrmBosOnvField;

    using FrmXonv = FrmXonvField;
    using BosXonv = BosXonvField;
    using FrmBosXonv = FrmBosXonvField;

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<uint_t mbf_type_ind>
    using mbf_t = typename std::tuple_element<mbf_type_ind, mbf_tup_t>::type;
    typedef mbf_t<c_mbf_type_ind> Mbf;

    typedef RdmIndsField RdmInds;
    typedef SpecMomIndsField SpecMomInds;

    /**
     * Fields are defined as symbols, i.e. members within a Row-derived class.
     * They also retain a pointer to the Row with which they are associated, and the Row in turn keeps a vector of
     * corresponding pointers to FieldBase objects. In some situations, we may want to specify a field within a row,
     * and them use that field but in the context of another row. In order to do this, the offset in memory between
     * the specified row and field is computed, then that offset is applied to the pointer to the other row, resulting
     * in a reference to the corresponding field in the other row.
     * This function is subject to run time checks ensuring safety of the casting operations.
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
        REQUIRE_TRUE(static_cast<FieldBase&>(field).belongs_to_row(source), "field arg must belong to source arg");
        auto target_ptr = reinterpret_cast<char*>(&target);
        auto source_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&field);
        long byte_offset = field_ptr - source_ptr;
        DEBUG_ASSERT_GT(byte_offset, 0l, "field pointer is not positively offset from row!");
        auto ptr = reinterpret_cast<field_t*>(target_ptr + byte_offset);
        DEBUG_ASSERT_TRUE(ptr->belongs_to_row(target), "field identification failed");
        return *ptr;
    }

    /**
     * same as the above function, but for references to CompositeFields instead
     * @tparam row_t
     *  Row-derived class containing the field_t symbol
     * @tparam Args
     *  FieldBase-derived arguments to the CompositeField class template
     * @param target
     *  Row object for which the CompositeField reference is desired
     * @param source
     *  Row object for which the CompositeField reference is specified
     * @param multifield
     *  CompositeField-derived object which is the referenced symbol within the source row
     * @return
     *  reference to the same CompositeField within "target" as is represented by "field" within "source"
     */
    template<typename row_t, typename ...Args>
    static CompositeField<Args...>& identify(row_t& target, row_t& source, CompositeField<Args...>& composite){
        static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
        auto target_ptr = reinterpret_cast<char*>(&target);
        auto source_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&composite);
        long byte_offset = field_ptr - source_ptr;
        DEBUG_ASSERT_GT(byte_offset, 0, "field pointer is not positively offset from row!");
        auto ptr = reinterpret_cast<CompositeField<Args...>*>(target_ptr + byte_offset);
        DEBUG_ASSERT_TRUE(ptr->belongs_to_row(target), "multifield identification failed");
        return *ptr;
    }

}

#endif //M7_FIELDS_H