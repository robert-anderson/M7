//
// Created by rja on 05/03/2020.
//

#ifndef M7_TABLENEW_H
#define M7_TABLENEW_H

#include <src/defs.h>
#include <iostream>
#include <assert.h>
#include <typeindex>
#include "src/utils.h"

/*
 * Tables both define the layout of data stored in a defs::data_t buffer,
 * and provide a namespace in which the stored data can be accessed by
 * symbols with the Field type.
 *
 * Specification of the data layout is done by building up the member fields
 * in the list initialization of the Table struct. Consecutively added
 * fields of the same data type are allowed to share a defs::data_t word,
 * with the exception of boolean fields of length >1 (aka bitfields). Bitfields
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
 * Table rows are always padded to fill a whole cache line
 *
 */

struct TableNew {
    std::vector<defs::data_t> m_data{};
    // the number of defs::data_t words required to store the encoded data
    size_t m_ndatawords = 0ul;
    /*
     * the number of defs::data_t words required to store the encoded data padded
     * to the smallest number of whole cache lines to prevent false sharing
     */
    size_t m_ndatawords_padded = 0ul;
    const size_t m_nsegment;
    size_t m_nrow_per_segment = 0;

    template<typename T>
    static size_t nbit_per_element() {
        return std::is_same<T, bool>::value ? 1 : sizeof(T) * 8;
    }

    size_t pair_to_irow(const defs::pair &pair) const {
        return pair.first*m_nrow_per_segment+pair.second;
    }

    size_t nrow() const {
        return m_nrow_per_segment*m_nsegment;
    }

    size_t ndataword_per_segment(size_t nrow_per_segment) const {
        return nrow_per_segment*m_ndatawords_padded;
    }
    size_t ndataword_per_segment() const {
        return ndataword_per_segment(m_nrow_per_segment);
    }


    size_t extend(size_t delta_nrow_per_segment) {
        /*
         * add more rows to each segment
         */
        m_data.resize(m_data.size() + delta_nrow_per_segment * m_nsegment * m_ndatawords_padded, 0);
        /*
         * move segments backwards to help avoid overlap. std::move will handle overlap correctly if it occurs
         */
        for (size_t isegment=m_nsegment-1; isegment>0; --isegment){
            std::move(
                    m_data.begin()+isegment*ndataword_per_segment(),
                    m_data.begin()+(isegment+1)*ndataword_per_segment(),
                    m_data.begin()+isegment*ndataword_per_segment(m_nrow_per_segment + delta_nrow_per_segment)
            );
        }
        m_nrow_per_segment+=delta_nrow_per_segment;
    }

    struct FieldBase {
        TableNew *m_table;
        // length of the indexer
        const size_t m_nelement;
        const size_t m_nbit_per_element;
        const std::type_index m_type_index;
        // first value is dataword offset, second is offset in elements from start of dataword
        const std::pair<size_t, size_t> m_offset;

        virtual std::string to_string(size_t irow) = 0;

        FieldBase(TableNew *table, size_t nbit, size_t nelement, const std::type_info &type_info) :
                m_table(table), m_nbit_per_element(nbit), m_nelement(nelement),
                m_type_index(type_info), m_offset(table->add_field(this)) {}

    };

    std::vector<FieldBase *> m_fields;

    std::pair<size_t, size_t> add_field(FieldBase *field) {
        if (m_nrow_per_segment) throw std::runtime_error("Cannot add a field to an initialised Table");
        size_t dataword_offset = 0ul;
        size_t element_offset = 0ul;
        if (!m_fields.empty()) {
            if (field->m_type_index != m_fields.back()->m_type_index) {
                // different type to last field
                dataword_offset = m_ndatawords;
                element_offset = 0;
            } else {
                // same type as last field
                element_offset = m_fields.back()->m_offset.second + m_fields.back()->m_nelement;
                element_offset *= field->m_nbit_per_element;
                dataword_offset = m_fields.back()->m_offset.first;
                dataword_offset += element_offset / nbit_per_element<defs::data_t>();
                element_offset %= nbit_per_element<defs::data_t>();
            }
        }
        m_ndatawords = dataword_offset * nbit_per_element<defs::data_t>();
        m_ndatawords += (element_offset + field->m_nelement) * field->m_nbit_per_element;
        m_ndatawords = integer_utils::divceil(m_ndatawords, nbit_per_element<defs::data_t>());
        m_ndatawords_padded = sizeof(defs::data_t)*integer_utils::divceil(m_ndatawords*sizeof(defs::data_t), defs::cache_line_size);
        m_fields.push_back(field);
        return {dataword_offset, element_offset};
    }

    TableNew(size_t nsegment=1) : m_nsegment(nsegment) {}

    defs::inds field_nelements();

    std::vector<std::pair<size_t, size_t>> field_offsets();

    void print(size_t irow);

    virtual void zero(){
        m_data.assign(m_data.size(), 0);
    }
};


#endif //M7_TABLENEW_H
