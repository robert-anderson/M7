//
// Created by rja on 19/02/23.
//

#include "HfExcitHists.h"

hf_excit_hist::IndVals::IndVals(const hdf5::NodeReader &parent, str_t name) :
        m_gr(parent, name),
        m_inds(m_gr, "indices", mpi::on_node_i_am_root()),
        m_vals(m_gr, "values", mpi::on_node_i_am_root()) {
    if (mpi::on_node_i_am_root()){
        auto order = m_vals.sort_inds(false, true);
        m_inds.reorder_rows(order);
        m_vals.reorder(order);
    }
}
