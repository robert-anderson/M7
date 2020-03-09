//
// Created by rja on 05/03/2020.
//

#ifndef M7_FIELDSET_H
#define M7_FIELDSET_H

#include <src/defs.h>
#include <iostream>
#include <assert.h>


/*
 * FieldSets both define the layout of data stored in a defs::data_t buffer,
 * and provides a namespace in which the stored data can be accessed by
 * variables with the Field type.
 *
 * Specification of the data layout is done by building up the member fields
 * in the list initialization of the FieldSet struct. Consecutively added 
 * fields of the same data type are allowed to share a defs::data_t word,
 * with the exception of boolean fields of length>1 (aka bitfields). Bitfields
 * are always allotted a whole number of data words.
 *
 * e.g. A: 3 floats
 *      B: 5 chars
 *      C: 6 chars
 *      D: 28 bools
 *      E: 1 bool
 *      F: 1 bool
 *      G: 1 bool
 *
 * |A0A0A0A0|A0A0A0A0|A0A0A0A0|A0A0A0A0|A1A1A1A1|A1A1A1A1|A1A1A1A1|A1A1A1A1|
 * |A2A2A2A2|A2A2A2A2|A2A2A2A2|A2A2A2A2|________|________|________|________|
 * |C0C0C0C0|C1C1C1C1|C2C2C2C2|C3C3C3C3|C4C4C4C4|C5C5C5C5|________|________|
 * |DDDDDDDD|DDDDDDDD|DDDDDDDD|DDDD____|________|________|________|________|
 * |EFG_____|________|________|________|________|________|________|________|
 *
 */

struct FieldSet {
    defs::data_t *m_buffer = nullptr;
    size_t m_nbit = 0ul;
    size_t m_length = 0ul;
    size_t m_type_offset = 0ul;
    size_t m_last_type = 0ul;

    template<size_t nind=1>
    struct FieldBase {
        FieldSet *m_field_set;
        const std::string m_description;
        const Indexer<nind> m_indexer;
        const size_t m_offset=0ul;

        virtual std::string to_string(size_t irow, size_t ientry) = 0;

        FieldBase(FieldSet *field_set, Indexer<nind> indexer, std::string description = "") :
                m_field_set(field_set), m_description(std::move(description)), m_indexer(indexer) {}
    };

    template<typename T, size_t nind=1>
    struct Field : FieldBase<nind> {
        const size_t m_offset;

        Field(FieldSet *field_set, Indexer<nind> indexer, std::string description = "") :
                FieldBase<nind>(field_set, indexer, description), m_offset(field_set->add_field(this, length)) {}

        Field(FieldSet *field_set, size_t length, std::string description = "") :
                Field(field_set, length, description){}

        T *get(const size_t &irow, const size_t &ientry = 0) {
            assert(ientry < m_length);
            return (T *) (m_field_set->m_buffer + irow * m_field_set->m_length) + m_offset + ientry;
        }

        std::string to_string(size_t irow, size_t ientry) override {
            return std::to_string(*get(irow, ientry));
        }
    };

    template<typename T>
    static size_t size_in_bits() { return std::is_same<bool, T>::value ? 1 : 8*sizeof(T); }

    template<typename T>
    size_t add_field(Field<T> *field, size_t length) {
        auto advance = [&](){
            auto tmp = m_nbit%size_in_bits<defs::data_t>();
            if (tmp) m_nbit+=size_in_bits<defs::data_t>()-tmp;
        };

        size_t this_type = typeid(T).hash_code();
        if (m_nbit>0){
            if (this_type!=m_last_type) advance();
        }

        auto offset = m_nbit/size_in_bits<T>();
        m_nbit+=length*size_in_bits<T>();
        m_length = m_nbit/size_in_bits<defs::data_t>();
        if (m_nbit%size_in_bits<defs::data_t>()) m_length++;
        m_last_type = this_type;
        m_fields.push_back(field);
        return offset;
    }

    std::vector<FieldBase *> m_fields;

    FieldSet(defs::data_t *buffer) : m_buffer(buffer) {}

    defs::inds field_lengths(){
        defs::inds result{};
        for (auto field: m_fields) result.push_back(field->m_length);
        return result;
    }

    defs::inds field_offsets(){
        defs::inds result{};
        for (auto field: m_fields) result.push_back(field->m_offset);
        return result;
    }

    void print(size_t irow) {
        std::cout << irow << " |";
        for (auto field: m_fields) {
            for (size_t ientry = 0ul; ientry < field->m_length; ++ientry) {
                std::cout << field->to_string(irow, ientry);
                std::cout << "   ";
            }
        }
        std::cout << std::endl;
    }
};


#endif //M7_FIELDSET_H
