//
// Created by anderson on 03/02/2023.
//

#include "Field.h"

void hdf5::field::save(const FieldBase &field, const hdf5::NodeWriter &nw, const str_t &name, hdf5::Type type,
                       bool is_complex, uintv_t item_shape, strv_t item_dim_names, bool this_rank,
                       uint_t max_nitem_per_op, std::list<Attr> attrs) {
    const dataset::ItemFormat item_format(type, std::move(item_shape), std::move(item_dim_names), is_complex);
    const dataset::PartDistListFormat format(item_format, this_rank ? field.m_row->m_table->nrecord() : 0ul, "row");
    REQUIRE_EQ(item_format.m_size, field.m_size, "item size implied by given shape inconsistent with field size");
    DatasetSaver ds(nw, name, format);

    // need not allocate a buffer longer than needed
    max_nitem_per_op = std::min(format.m_local.m_nitem, max_nitem_per_op);
    v_t<buf_t> buf(max_nitem_per_op * item_format.m_size);

    // items are given by records, which are rows below the high water mark that have not been freed
    uint_t nitem_found = 0ul;
    uint_t irow = 0ul;
    bool all_done = false;
    while (!all_done) {
        buf.clear();
        const auto next_nitem_found = std::min(nitem_found + max_nitem_per_op, format.m_local.m_nitem);
        const void* src = nullptr;
        if (next_nitem_found != nitem_found) {
            auto buf_ptr = buf.data();
            while (nitem_found != next_nitem_found) {
                if (!field.m_row->m_table->is_freed(irow)) {
                    field.to_buffer(buf_ptr, irow);
                    buf_ptr += field.m_size;
                    ++nitem_found;
                }
                ++irow;
            }
            src = buf.data();
        }
        all_done = ds.write(src, ds.nitem_next(max_nitem_per_op));
    }
    ds.save_attrs(attrs);
}
