//
// Created by rja on 18/05/2021.
//

#ifndef M7_SUBSPACE_H
#define M7_SUBSPACE_H

#include <src/core/wavefunction/Wavefunction.h>

/**
 * Set of MappedTable keys distributed over all MPI ranks.
 * @tparam row_t
 *  data layout of the mapped table whose keys are to be tracked
 */
template<typename row_t>
struct ProtectedKeySet : RowProtector {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    typedef typename KeyField<row_t>::type key_field_t;
    struct KeyRow : Row {
        /**
         * key field
         */
        key_field_t m_key;
        /**
         * row index of key within the mapped table ()
         */
        field::Number<size_t> m_row_id;
    };

    BufferedTable<KeyRow, true> m_local;
    BufferedTable<KeyRow, true> m_global;

    void add_(row_t& row){
        protect(row.m_i);
    }

    ProtectedKeySet(MappedTable<row_t>& table): RowProtector(table){}
};




/**
 * Network of a specific excitation level
 */
struct SubspaceExlvl {

};

/**
 * All connections are stored between elements in the DynamicRowSet.
 *
 * The required information is
 */
struct Subspace : Wavefunction::DynamicRowSet {

    //std::vector<std::forward_list<>> m_

    Subspace(const Wavefunction& wf, std::string name): Wavefunction::DynamicRowSet(wf, name){}

    void update() override;

};



#endif //M7_SUBSPACE_H
