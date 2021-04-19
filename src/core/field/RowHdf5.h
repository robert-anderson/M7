//
// Created by rja on 17/03/2021.
//

#ifndef M7_ROWHDF5_H
#define M7_ROWHDF5_H

#include "Row.h"
#include "FieldBase.h"

struct RowHdf5Base {
    const size_t m_nitem, m_nitem_max, m_nitem_total;
    std::vector<std::string> m_field_names;
    defs::inds m_selected_field_inds;

    RowHdf5Base(const Row &row, size_t nitem, std::vector<std::string> field_names) :
            m_nitem(nitem), m_nitem_max(mpi::all_max(m_nitem)),
            m_nitem_total(mpi::all_sum(m_nitem)), m_field_names(field_names) {
        m_selected_field_inds.reserve(m_field_names.size());
        for (const auto &field_name: m_field_names) {
            MPI_REQUIRE(!field_name.empty(), "Selected field name must be non-zero in length");
            bool match_found = false;
            size_t i =0ul;
            for (FieldBase *field_ptr : row.m_fields) {
                if (field_ptr->m_name == field_name) match_found = true, m_selected_field_inds.emplace_back(i);
                ++i;
            }
            MPI_REQUIRE(match_found, "Invalid field name \"" + field_name + "\"");
        }
    }

protected:
    std::vector<std::string> get_all_field_names(const Row &row) const {
        std::vector<std::string> out;
        out.reserve(row.m_fields.size());
        for (FieldBase *field_ptr: row.m_fields) {
            MPI_REQUIRE(!field_ptr->m_name.empty(), "All fields selected for HDF5 write must be named");
            out.emplace_back(field_ptr->m_name);
        }
        return out;
    }
};

template<typename row_t>
struct RowHdf5Writer : RowHdf5Base, row_t {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    hdf5::GroupWriter m_group;
    std::vector<hdf5::NdListWriter> m_writers;

    RowHdf5Writer(const row_t &row, hdf5::GroupWriter &parent, std::string name, size_t nitem,
                  std::vector<std::string> field_names) :
            RowHdf5Base(row, nitem, field_names), row_t(row), m_group(name, parent) {
        m_writers.reserve(m_selected_field_inds.size());
        for (auto &ifield: m_selected_field_inds) {
            auto field = Row::m_fields[ifield];
            std::string field_name = field->m_name;
            m_writers.emplace_back(m_group, field_name, field->h5_shape(), nitem, field->h5_type(),
                                   field->h5_dim_names());
            /*
             * the writer creates the dataset, to which the field's attributes,
             * if any, must be added
             */
            field->h5_write_attrs(m_writers.back().m_dataset_handle);
        }
    }


public:

    RowHdf5Writer(const row_t &row, hdf5::GroupWriter &parent, std::string name, size_t nitem) :
            RowHdf5Writer(row, parent, name, nitem, get_all_field_names(row)) {}

    virtual void write(const size_t &iitem) {
        for (size_t ifield = 0ul; ifield < m_selected_field_inds.size(); ++ifield) {
            auto &writer = m_writers[ifield];
            auto field = Row::m_fields[ifield];
            if (iitem < m_nitem) field->h5_write(writer, iitem);
            else {
                // runoff write operation for collective I/O
                writer.write_h5item_bytes(iitem, nullptr);
            }
        }
    }
};


template<typename row_t>
struct RowHdf5Reader : RowHdf5Base, row_t {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    hdf5::GroupReader m_group;
    std::vector<hdf5::NdListReader> m_readers;

private:

    size_t get_nitem(const Row &row, hdf5::GroupReader &parent, std::string name) const {
        if (row.m_fields.empty()) return 0ul;
        hdf5::GroupReader group(name, parent);
        FieldBase *field = row.m_fields[0];
        hdf5::NdListReader reader(group, field->m_name, field->h5_type());
        return reader.m_nitem_local;
    }

public:

    RowHdf5Reader(const row_t &row, hdf5::GroupReader &parent, std::string name, std::vector<std::string> field_names) :
            RowHdf5Base(row, get_nitem(row, parent, name), field_names), row_t(row), m_group(name, parent) {
        m_readers.reserve(m_selected_field_inds.size());
        for (auto &ifield: m_selected_field_inds) {
            auto &field = Row::m_fields[ifield];
            std::string field_name = field->m_name;
            m_readers.emplace_back(m_group, field_name, field->h5_type());
            // TODO: read attrs
        }
    }

    RowHdf5Reader(const row_t &row, hdf5::GroupReader &parent, std::string name) :
            RowHdf5Reader(row, parent, name, get_all_field_names(row)) {}

    virtual void read(const size_t &iitem) {
        for (size_t ifield = 0ul; ifield < m_readers.size(); ++ifield) {
            auto &reader = m_readers[ifield];
            auto &field = Row::m_fields[ifield];
            if (iitem < m_nitem) field->h5_read(reader, iitem);
            else {
                // runoff read operation for collective I/O
                reader.read_h5item_bytes(iitem, nullptr);
            }
        }
    }
};


#if 0
template<typename row_t>
struct RowHdf5Reader : row_t {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    hdf5::GroupReader m_group;
    std::vector<hdf5::NdListReader> m_readers;

    size_t nitem() const {
        return m_readers[0].m_nitem_local;
    }

    RowHdf5Reader(const row_t &row, hdf5::GroupReader& parent, std::string name):
            row_t(row), m_group(name, parent){
        m_readers.reserve(Row::m_fields.size());
        //TODO: give fields names
        size_t ifield = 0ul;
        for (auto& field: Row::m_fields)
            m_readers.emplace_back(m_group, "field_"+std::to_string(ifield++), field->h5_type());
    }
    void read(const size_t& iitem){
        for (size_t ifield=0ul; ifield<Row::m_fields.size(); ++ifield){
            auto& reader = m_readers[ifield];
            auto& field = Row::m_fields[ifield];
            field->h5_read(reader, iitem);
        }
    }
};
#endif


#endif //M7_ROWHDF5_H
