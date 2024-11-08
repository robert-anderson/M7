//
// Created by rja on 02/05/23.
//

#ifndef M7_DISTTABLELOADER_H
#define M7_DISTTABLELOADER_H

#include "M7_lib/field/FieldBase.h"
#include "DatasetTransaction.h"

/**
 * helper class to load a distributed MappedTable from fields stored as HDF5 datasets
 */
struct DistTableLoader {

    struct NameFieldPair {
        const str_t m_name;
        FieldBase * const m_field;
        NameFieldPair(str_t name, FieldBase *field): m_name(std::move(name)), m_field(field){}
    };

private:
    struct FieldLoaderPair {
        FieldBase * const m_field;
        hdf5::DatasetLoader * const m_loader;
    };

    /**
     * allocate a dataset loader for each field
     */
    v_t<hdf5::DatasetLoader> m_loaders;
    /**
     * put field and loader pointer pairs together in a vector so they can be easily iterated over simultaneously
     */
    v_t<FieldLoaderPair> m_fields_loaders;

    uint_t max_item_size() const;
    uint_t nitem_next(uint_t max_nitem_per_op) const;

    static v_t<NameFieldPair> make_pairs(const v_t<FieldBase*>& fields);

public:

    uint_t nitem() const;

    uint_t nitem_local() const;

    DistTableLoader(const hdf5::NodeReader& nr, v_t<NameFieldPair> pairs);

    DistTableLoader(const hdf5::NodeReader& nr, const v_t<FieldBase *>& fields);

    DistTableLoader(const hdf5::NodeReader& nr, Row& row): DistTableLoader(nr, row.m_fields){}


    template<typename fn_t>
    void load(uint_t max_nitem_per_op, const fn_t& fn) {
        functor::assert_prototype<void(uint_t)>(fn);
        max_nitem_per_op = std::min(max_nitem_per_op, nitem_local());
        for (auto field_loader : m_fields_loaders) {
            auto& table = field_loader.m_field->m_row->m_table;
            if (table->nrow_in_use() < max_nitem_per_op) table->push_back(max_nitem_per_op);
        }

        // use same contiguous buffer for all fields, so allocate enough space for the field with the largest items
        v_t<buf_t> buf(max_nitem_per_op * max_item_size());

        bool all_done = false;
        while (!all_done) {
            buf.clear();
            const auto nitem_to_find = nitem_next(max_nitem_per_op);
            for (auto& field_loader: m_fields_loaders) {
                auto& field = field_loader.m_field;
                auto& loader = field_loader.m_loader;
                all_done = loader->read(nitem_to_find ? buf.data() : nullptr, nitem_to_find);
                auto src = buf.data();
                for (uint_t irow = 0ul; irow < nitem_to_find; ++irow) {
                    REQUIRE_TRUE(field->check_buffer(src),
                        logging::format("data in buffer is invalid for the {} field", field->m_name));
                    field->from_buffer(src, irow);
                    src += field->m_size;
                }
            }
            fn(nitem_to_find);
        }
    }
};

#endif //M7_DISTTABLELOADER_H
