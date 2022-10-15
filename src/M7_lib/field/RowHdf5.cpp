//
// Created by Robert J. Anderson on 17/03/2021.
//

#include "RowHdf5.h"

RowHdf5Base::RowHdf5Base(const Row &row, uint_t nitem, strv_t field_names) :
        m_nitem(nitem), m_nitem_max(mpi::all_max(m_nitem)),
        m_nitem_total(mpi::all_sum(m_nitem)), m_field_names(std::move(field_names)),
        m_selected_field_inds(make_selected_field_inds(row.m_fields)) {}


uintv_t RowHdf5Base::make_selected_field_inds(const v_t<FieldBase *> &fields) const {
    uintv_t inds;
    inds.reserve(m_field_names.size());
    for (const auto &field_name: m_field_names) {
        REQUIRE_TRUE(!field_name.empty(), "Selected field name must be non-zero in length");
        bool match_found = false;
        uint_t i = 0ul;
        for (FieldBase *field_ptr: fields) {
            if (field_ptr->m_name == field_name) {
                match_found = true;
                inds.emplace_back(i);
            }
            ++i;
        }
        REQUIRE_TRUE(match_found, "Invalid field name \"" + field_name + "\"");
    }
    return inds;
}

void RowHdf5WriterBase::write_all_fields() {
    for (uint_t ifield = 0ul; ifield < m_selected_field_inds.size(); ++ifield) {
        auto &writer = m_column_writers[ifield];
        auto field = m_row.m_fields[ifield];
        field->save(writer, m_iitem);
    }
}

void RowHdf5WriterBase::null_write_all_fields() {
    for (auto &writer: m_column_writers) writer.write_h5item_bytes(m_iitem, nullptr);
}

bool RowHdf5WriterBase::write(const uint_t &n) {
    logging::debug_("beginning HDF5 write loop over rows");
    if (!m_iitem) m_row.restart();
    auto limit = std::min(m_iitem + n, m_nitem_max);
    while (m_iitem < limit){
        if (m_iitem >= m_nitem || m_row.is_h5_write_exempt()) null_write_all_fields();
        else write_all_fields();
        if (m_iitem < m_nitem) ++m_row;
        ++m_iitem;
    }
    mpi::barrier();
    logging::debug_("ending HDF5 write loop over rows");
    return mpi::all_lor(m_iitem >= m_nitem);
}

void RowHdf5WriterBase::write() {
    write(m_nitem_max);
}

RowHdf5WriterBase::RowHdf5WriterBase(Row &row, const hdf5::NodeWriter &parent, str_t name, uint_t nitem,
                                     strv_t field_names) :
        RowHdf5Base(row, nitem, std::move(field_names)), m_row(row), m_group(parent, name){
    m_column_writers.reserve(m_selected_field_inds.size());
    for (auto &ifield: m_selected_field_inds) {
        auto field = m_row.m_fields[ifield];
        str_t field_name = field->m_name;
        m_column_writers.emplace_back(m_group, field_name, field->h5_shape(), nitem, field->h5_type(),
                                      field->h5_dim_names());
        /*
         * the writer creates the dataset, to which the field's attributes, if any, must be added
         */
        field->h5_write_attrs(m_column_writers.back().m_dataset_handle);
    }
}

strv_t RowHdf5WriterBase::get_all_field_names(const Row &row) {
    strv_t out;
    out.reserve(row.m_fields.size());
    for (FieldBase *field_ptr: row.m_fields) {
        REQUIRE_FALSE(field_ptr->m_name.empty(), "All fields selected for HDF5 write must be named");
        out.emplace_back(field_ptr->m_name);
    }
    return out;
}


RowHdf5ReaderBase::RowHdf5ReaderBase(Row &row, const hdf5::NodeReader &parent, str_t name,
                                     strv_t field_names) :
        RowHdf5Base(row, get_nitem(row, parent, name, field_names), field_names),
        m_row(row), m_group(parent, name) {
    m_column_readers.reserve(m_selected_field_inds.size());
    for (auto &ifield: m_selected_field_inds) {
        auto &field = m_row.m_fields[ifield];
        str_t field_name = field->m_name;
        m_column_readers.emplace_back(m_group, field_name, field->h5_type());
        // TODO: read attrs
    }
}


strv_t RowHdf5ReaderBase::stored_field_names(
        const Row &row, const hdf5::NodeReader &parent, str_t name) {
    hdf5::GroupReader group(parent, name);
    strv_t out;
    out.reserve(row.m_fields.size());
    for (FieldBase *field: row.m_fields) {
        REQUIRE_FALSE(field->m_name.empty(), "All fields selected for HDF5 write must be named");
        if (group.child_exists(field->m_name)) out.emplace_back(field->m_name);
    }
    return out;
}


uint_t RowHdf5ReaderBase::get_nitem(const Row &row, const hdf5::NodeReader& parent, str_t name,
                                    strv_t field_names) {
    if (row.m_fields.empty()) return 0ul;
    hdf5::GroupReader group(parent, name);
    // take the number of items in the first selected column
    auto ifield = group.first_existing_child(field_names);
    REQUIRE_LE(ifield, field_names.size(), "no selected field not found in HDF5 group");
    return hdf5::NdDistListReader::local_nitem(group.m_handle, field_names[ifield]);
}

void RowHdf5ReaderBase::read_all_fields() {
    for (uint_t ifield = 0ul; ifield < m_column_readers.size(); ++ifield) {
        auto &reader = m_column_readers[ifield];
        auto &field = m_row.m_fields[ifield];
        field->load(reader, m_iitem);
    }
}

void RowHdf5ReaderBase::null_read_all_fields() {
    for (auto & reader : m_column_readers) reader.read_h5item_bytes(m_iitem, nullptr);
}

bool RowHdf5ReaderBase::read(const uint_t &n) {
    logging::debug_("beginning HDF5 read loop over multidimensional list dataset");
    if (!m_iitem) {
        // first call, so reset table and row
        m_row.m_table->clear();
        m_row.m_table->push_back(m_nitem);
        m_row.restart();
    }
    auto limit = std::min(m_iitem + n, m_nitem_max);
    while (m_iitem < limit) {
        if (m_iitem >= m_nitem) null_read_all_fields();
        else read_all_fields();
        if (m_iitem < m_nitem) ++m_row;
        ++m_iitem;
    }
    mpi::barrier();
    logging::debug_("ending HDF5 read loop over multidimensional list dataset");
    return mpi::all_lor(m_iitem >= m_nitem);
}

void RowHdf5ReaderBase::read() {
    read(m_nitem_max);
}