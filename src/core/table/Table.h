//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_TABLE_H
#define SANDBOX2_TABLE_H


#include <vector>
#include <cstring>
#include <functional>
#include "src/core/util/defs.h"
#include "src/core/util/utils.h"
#include "src/core/thread/AlignedAllocator.h"
#include "Field.h"

class Table {
    const std::string m_name;
protected:
    // data buffer
    std::vector<defs::data_t, AlignedAllocator<defs::data_t, defs::cache_line_size>> m_data;
    // all associated fields
    std::vector<Field *> m_fields;
    // size of one row in bytes
    size_t m_row_size = 0;
    // size of each row in bytes padded to cache line width
    size_t m_padded_row_size = 0;
    // size of each row in data_t words padded to cache line width
    size_t m_padded_row_dsize = 0;
    // data buffers may be divided into segments of an equal number of rows
    const size_t m_nsegment;
    // the number of rows in a table segment, may be increased by a call to "expand"
    size_t m_nrow_per_segment = 0;
    // size of each segment in bytes (recomputed each expand)
    size_t m_segment_size = 0;
    // size of each segment in data_t words (recomputed each expand)
    size_t m_segment_dsize = 0;
    // offset of the beginning of each segment from the beginning of the table
    defs::inds m_segment_doffsets;
    // for stats output, debugging and other reporting tasks:
    std::vector<std::function<std::string(const Element *)>> m_printers;
public:

    Table(std::string name, size_t nsegment = 1);

    char *field_begin(const Field *field, const size_t &irow, const size_t isegment = 0);

    char *row_begin(const size_t &irow, const size_t isegment = 0);

    virtual void expand(size_t delta_nrow);

    virtual void resize(size_t nrow);

    size_t irow(const size_t &irow, const size_t &isegment = 0) const;

    virtual void zero();

    void zero_row(const size_t &irow, const size_t &isegment);

    size_t add_field(Field *field);

    void update_last_field();

    void print();

    const size_t &nrow_per_segment() const;

    bool compatible_with(const Table &other) const;

    bool is_allocated() const;

    std::string row_to_string(size_t irow, size_t isegment);

    virtual std::string to_string();

    virtual std::string to_string(const defs::inds &nrows);

    void print_row(size_t irow, size_t isegment);

    size_t dsize() const;

private:
    void update_row_size(size_t size);

    void increment_row_size(size_t delta);

    void roundup_row_size();

    void update_nrow_per_segment(size_t nrow);

    void increment_nrow_per_segment(size_t delta);

};


#endif //SANDBOX2_TABLE_H
