//
// Created by Robert J. Anderson on 17/03/2021.
//

#ifndef M7_ROWHDF5_H
#define M7_ROWHDF5_H

#include <utility>
#include <M7_lib/hdf5/Group.h>

#include "Row.h"
#include "FieldBase.h"

/**
 * For the sake of efficient communication, the rows of a table contain diverse data types, stored contiguously on a
 * byte-typed buffer. HDF5 stores typed data, and for this we must transpose the data stored in the table. This is to
 * make the HDF5 archives produced by M7 more useful as an output format. Thus, it is necessary to transpose the table
 * prior to output.
 */
struct RowHdf5Base {
    /**
     * number of rows stored on this rank
     */
    const uint_t m_nitem;
    /**
     * largest number of items stored on a single rank
     */
    const uint_t m_nitem_max;
    /**
     * total number of items across all ranks
     */
    const uint_t m_nitem_total;
    /**
     * names of the HDF5 datasets / row fields to be read or written
     */
    const std::vector<std::string> m_field_names;
    /**
     * field names are used to extract the positions in the Row::m_fields member, and then subsequent accesses to
     * selected fields are made via this resulting index vector.
     */
    const uintv_t m_selected_field_inds;
    /**
     * index of next item to be read or written
     */
    uint_t m_iitem = 0ul;

    RowHdf5Base(const Row& row, uint_t nitem, std::vector<std::string> field_names);

protected:

    uintv_t make_selected_field_inds(const std::vector<FieldBase *>& fields) const;

};

struct RowHdf5ReaderBase : RowHdf5Base {
    /**
     * reference to superclass instance cast to a Row
     */
    Row& m_row;
    hdf5::GroupReader m_group;
    std::vector<hdf5::NdDistListReader> m_column_readers;

    RowHdf5ReaderBase(Row& row, const hdf5::NodeReader& parent, std::string name, std::vector<std::string> field_names);

protected:

    /**
     * when no field_names vector is specified, we assume all datasets in the HDF5 group which are named after a field
     * are to be selected
     * @param row
     *  row from which to extract all field names
     * @param parent
     *  HDF5 group in which to look for field names among the dataset names
     * @return
     *  all field names
     */
    static std::vector<std::string> stored_field_names(const Row& row, const hdf5::NodeReader& parent, std::string name);

private:

    static uint_t get_nitem(const Row& row, const hdf5::NodeReader& parent,
                            std::string name, std::vector<std::string> field_names);

    /**
     * reads actual data. called in the case that m_iitem < m_nitem
     */
    void read_all_fields();

    /**
     * runoff read operation required for collective I/O. called in the case that m_iitem >= m_nitem
     */
    void null_read_all_fields();

public:

    /**
     * read a block of rows from HDF5 group
     * @param n
     *  number of rows to read
     * @return
     *  true if there remain unread rows on any process
     */
    bool read(const uint_t& n);

    /**
     * read all rows from HDF5 group
     */
    void read();
};

struct RowHdf5WriterBase : RowHdf5Base {
    /**
     * reference to superclass instance cast to a Row
     */
    Row& m_row;
    hdf5::GroupWriter m_group;
    std::vector<hdf5::NdDistListWriter> m_column_writers;

    RowHdf5WriterBase(Row& row, const hdf5::NodeWriter& parent, std::string name, uint_t nitem,
                      std::vector<std::string> field_names);

protected:
    /**
     * when no field_names vector is specified, we assume all fields in the row are to be selected
     * @param row
     *  row from which to extract all field names
     * @return
     *  all field names
     */
    static std::vector<std::string> get_all_field_names(const Row& row);

private:
    /**
     * writes actual data. called in the case that m_iitem < m_nitem
     */
    void write_all_fields();

    /**
     * runoff write operation required for collective I/O. called in the case that m_iitem >= m_nitem
     */
    void null_write_all_fields();

public:

    /**
     * write a block of rows from HDF5 group
     * @param n
     *  number of rows to write
     * @return
     *  true if there remain unwritten rows on any process
     */
    bool write(const uint_t& n);

    void write();

};



template<typename row_t>
struct RowHdf5Reader : row_t, RowHdf5ReaderBase {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");

    RowHdf5Reader(const row_t& row, const hdf5::NodeReader& parent, std::string name, std::vector<std::string> field_names) :
            row_t(row), RowHdf5ReaderBase(*this, parent, name, field_names) {}

    RowHdf5Reader(const row_t& row, const hdf5::NodeReader& parent, std::string name) :
            RowHdf5Reader(row, parent, name, stored_field_names(row, parent, name)) {}
};

/**
 * creates a copy of a row in the table to be output, and a vector of column writers. Then each call to write method
 * writes another item of each column to the archive.
 * @tparam row_t
 *  specific row-derived type
 */
template<typename row_t>
struct RowHdf5Writer : row_t, RowHdf5WriterBase {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");

    RowHdf5Writer(const row_t& row, const hdf5::NodeWriter& parent, std::string name, uint_t nitem,
                  std::vector<std::string> field_names) :
            row_t(row), RowHdf5WriterBase(*this, parent, name, nitem, field_names) {}

    RowHdf5Writer(const row_t& row, const hdf5::NodeWriter& parent, std::string name, uint_t nitem) :
            RowHdf5Writer(row, parent, name, nitem, get_all_field_names(row)) {}
};

#endif //M7_ROWHDF5_H
