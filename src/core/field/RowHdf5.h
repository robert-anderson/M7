//
// Created by rja on 17/03/2021.
//

#ifndef M7_ROWHDF5_H
#define M7_ROWHDF5_H

#include "Row.h"
#include "FieldBase.h"

template<typename row_t>
struct RowHdf5Writer : row_t {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    const size_t m_nitem, m_nitem_max;
    hdf5::GroupWriter m_group;
    std::vector<hdf5::NdListWriterBase> m_writers;
    RowHdf5Writer(const row_t &row, hdf5::GroupWriter& parent, std::string name, size_t nitem):
            row_t(row), m_nitem(nitem), m_nitem_max(mpi::all_max(m_nitem)), m_group(name, parent){
        m_writers.reserve(Row::m_fields.size());
        size_t iuntitled = 0ul;
        for (auto& field: Row::m_fields) {
            std::string field_name = field->m_name;
            if (field_name.empty()) field_name = "untitled_" + std::to_string(iuntitled++);
            m_writers.emplace_back(m_group, field_name, field->h5_shape(), nitem, field->h5_type(), field->h5_dim_names());
            /*
             * the writer creates the dataset, to which the fields attributes,
             * if any, must be added
             */
            field->h5_write_attrs(m_writers.back().m_dataset_handle);
        }
    }
    void write(const size_t& iitem){
        for (size_t ifield=0ul; ifield<Row::m_fields.size(); ++ifield){
            auto& writer = m_writers[ifield];
            auto& field = Row::m_fields[ifield];
            if (iitem<m_nitem) field->h5_write(writer, iitem);
            else {
                // runoff write operation for collective I/O
                writer.write_h5item_bytes(iitem, nullptr);
            }
        }
    }
};

template<typename row_t>
struct RowHdf5Reader : row_t {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    hdf5::GroupReader m_group;
    std::vector<hdf5::NdListReaderBase> m_readers;

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


#endif //M7_ROWHDF5_H
