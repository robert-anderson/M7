//
// Created by rja on 19/05/2021.
//

#ifndef M7_ROWPROTECTOR_H
#define M7_ROWPROTECTOR_H

#include <src/core/parallel/RankAllocator.h>

/**
 * keeps track of rows which must not be erased
 */
struct RowProtector {
    /**
     * Table containing rows which can be protected by this object
     */
    const TableBase& m_table;
    /**
     * stores protected status of every row in m_table. should be implemented as a bitmap
     */
    std::vector<bool> m_flags;
    /**
     * number of true values in m_flags
     */
    size_t m_nprotected = 0ul;
    /**
     * iterator in m_table::m_row_protectors associated with this object
     */
    typename std::list<RowProtector *>::iterator m_it;

    RowProtector(const TableBase& table);

    ~RowProtector();
    /**
     * mark row as protected
     * @param irow
     *  row index to protect
     */
    void protect(const size_t &irow);
    /**
     * cease treating row as protected
     * @param irow
     *  row index to release
     */
    void release(const size_t &irow);
    /**
     * called when the associated m_table is resized
     * @param nrow
     */
    void on_resize(size_t nrow);
    /**
     * @param irow
     *  row index
     * @return
     *  true if indexed row is protected
     */
    bool is_protected(const size_t &irow) const;
    /**
     * @return
     *  true if this object protects any rows
     */
    bool is_protected() const;
};


#endif //M7_ROWPROTECTOR_H
