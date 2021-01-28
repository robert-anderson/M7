//
// Created by rja on 21/10/2020.
//

#include "Table.h"
#include "src/core/io/Logging.h"

Table::Table() {}

Table::Table(const Table &other) :
        m_row_size(other.m_row_size), m_row_dsize(other.m_row_dsize),
        m_current_byte_offset(other.m_current_byte_offset){
    other.m_last_copied = this;
}

void Table::set_buffer(Buffer *buffer) {
    ASSERT(buffer);
    ASSERT(!m_bw.allocated())
    buffer->append_window(&m_bw);
}

bool Table::is_full() const {
    return m_hwm==m_nrow;
}

size_t Table::push_back(size_t nrow) {
    if (m_hwm+nrow>m_nrow) expand(nrow);
    auto tmp = m_hwm;
    m_hwm+=nrow;
    return tmp;
}

size_t Table::get_free_row() {
    if (m_free_rows.empty()) return push_back();
    auto irow = m_free_rows.top();
    m_free_rows.pop();
    return irow;
}

defs::data_t *Table::dbegin() {
    return m_bw.dbegin();
}

const defs::data_t *Table::dbegin() const {
    return m_bw.dbegin();
}

defs::data_t *Table::dbegin(const size_t &irow) {
    ASSERT(irow<m_hwm)
    ASSERT(m_bw.dbegin())
    return m_bw.dbegin() + irow * m_row_dsize;
}

const defs::data_t *Table::dbegin(const size_t &irow) const {
    ASSERT(irow<m_hwm)
    ASSERT(m_bw.dbegin())
    return m_bw.dbegin() + irow * m_row_dsize;
}

size_t Table::add_column(const ColumnBase *column) {
    // returns the offset in bytes for the column being added
    auto offset = 0ul;
    if(!m_columns.empty()){
        offset = m_columns.back()->m_offset+m_columns.back()->m_size;
        if (!m_columns.back()->is_same_type_as(*column)){
            // go to next whole dataword
            offset = integer_utils::divceil(offset, defs::nbyte_data)*defs::nbyte_data;
        }
    }

    m_current_byte_offset = offset + column->m_size;
    m_row_dsize = integer_utils::divceil(m_current_byte_offset, defs::nbyte_data);
    m_row_size = m_row_dsize*defs::nbyte_data;

    m_columns.push_back(column);
    return offset;
}

void Table::clear() {
    if (!m_bw.allocated()) return;
    std::memset(dbegin(), 0, m_row_size * m_hwm);
    m_hwm = 0ul;
    while (!m_free_rows.empty()) m_free_rows.pop();
}

void Table::clear(const size_t &irow) {
    std::memset(dbegin(irow), 0, m_row_size);
    m_free_rows.push(irow);
}

bool Table::is_cleared() const {
    return std::all_of(dbegin(), dbegin()+m_row_dsize*m_nrow, [](const defs::data_t& i){return i==0;});
}

bool Table::is_cleared(const size_t &irow) const {
    return std::all_of(dbegin(irow), dbegin(irow)+m_row_dsize, [](const defs::data_t& i){return i==0;});
}

std::string Table::column_details(size_t width) const {
    std::string res;
    for (size_t i = 0ul; i < m_columns.size(); ++i) {
        std::string desc;
        if (!m_columns[i]->m_description.empty()) desc = " (\""+m_columns[i]->m_description+"\")";
        res += "\nField " + std::to_string(i) + desc + ":\n";
        for (auto pair:m_columns[i]->m_data.m_details)
            res += "\t" + utils::padded_string(pair.first, width) + ": " + utils::padded_string(pair.second, width)+"\n";
    }
    return res;
}

void Table::print_column_details(size_t width) const {
    std::cout << column_details(width) << std::endl;
}

size_t Table::bw_dsize() const {
    return m_bw.dsize();
}

void Table::print_contents(const defs::inds *ordering) const {
    const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
    for (size_t iirow=0ul; iirow<n; ++iirow){
        auto irow = ordering ? (*ordering)[iirow] : iirow;
        std::cout << irow << ". ";
        for (auto column: m_columns){
            std::cout << column->to_string(irow)+" ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void Table::print_contents(const ExtremalValues &xv) const {
    defs::inds tmp;
    tmp.reserve(xv.nfound());
    for (size_t i=0ul; i<xv.nfound(); ++i) tmp.push_back(xv[i]);
    print_contents(&tmp);
}

void Table::resize(size_t nrow) {
    assert(nrow>m_nrow);
    m_bw.resize(nrow);
    m_nrow = nrow;
}

void Table::expand(size_t nrow, double expansion_factor) {
    m_bw.expand(nrow, expansion_factor);
    m_nrow = m_bw.nrow();
}

void Table::expand(size_t nrow) {
    m_bw.expand(nrow);
    m_nrow = m_bw.nrow();
}

void Table::erase_rows(const defs::inds &irows, const cb_list_t& callbacks) {
    for (auto irow : irows) {
        for (auto f: callbacks) f(irow);
        clear(irow);
    }
}

void Table::insert_rows(const Buffer &recv, const cb_list_t& callbacks) {
    const auto nrow = recv.dsize()/m_row_dsize;
    auto irow_first = push_back(nrow);
    std::memcpy(dbegin(irow_first), recv.dbegin(), nrow*m_row_size);
    for (size_t irow = irow_first; irow<irow_first+nrow; ++irow){
        for (auto f: callbacks) f(irow);
    }
}

void Table::send_rows(const defs::inds &irows, size_t irank_dst, const cb_list_t& callbacks) {
    size_t nrow;
    nrow = irows.size();
    log::info_("Transferring {} rows outward to rank {}", nrow, irank_dst);
    // make rows to be sent contiguous in memory
    Buffer send("Outward transfer buffer", 1, m_row_dsize);
    send.resize(nrow);
    for (auto iirow = 0ul; iirow < nrow; ++iirow) {
        const auto &irow = irows[iirow];
        std::memcpy(send.dbegin() + iirow * m_row_dsize, dbegin(irow), m_row_size);
    }

    mpi::send(&nrow, 1, irank_dst, 0);
    mpi::send(send.dbegin(), m_row_dsize * nrow, irank_dst, 1);
    /*
     * sent rows can now be erased, after all callbacks are called
     */
    erase_rows(irows, callbacks);
}

void Table::recv_rows(size_t irank_src, const cb_list_t& callbacks){
    size_t nrow;
    mpi::recv(&nrow, 1, irank_src, 0);
    log::info_("Transferring {} rows inward from rank {}", nrow, irank_src);
    Buffer recv("Inward transfer buffer", 1, m_row_dsize*nrow);
    mpi::recv(recv.dbegin(), m_row_dsize * nrow, irank_src, 1);
    /*
     * now emplace received rows in table buffer window, and call all callbacks for each
     */
    insert_rows(recv, callbacks);
}

bool Table::has_compatible_format(const Table &other) {
    if (other.m_row_size != m_row_size) return false;
    if (other.m_columns.size() != m_columns.size()) return false;
    if (other.m_current_byte_offset != m_current_byte_offset) return false;
    for (size_t ifield=0ul; ifield<m_columns.size(); ++ifield){
        auto this_field = m_columns[ifield];
        auto other_field = m_columns[ifield];
        if (!other_field->is_same_type_as(*this_field)) return false;
        if (other_field->m_nelement != this_field->m_nelement) return false;
        if (other_field->m_size != this_field->m_size) return false;
        if (other_field->m_offset != this_field->m_offset) return false;
    }
    return true;
}

void Table::copy_row_in(const Table &src, size_t irow_src, size_t irow_dst) {
    ASSERT(irow_dst<m_hwm);
    ASSERT(has_compatible_format(src));
    std::memcpy(dbegin(irow_dst), src.dbegin(irow_src), m_row_size);
}