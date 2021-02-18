//
// Created by rja on 29/11/2020.
//

#ifndef M7_INTERMEDIATETABLE_H
#define M7_INTERMEDIATETABLE_H

#if 0
#include "src/core/table/Table.h"
#include "src/core/field/Fields.h"

struct IntermediateTable : public Table {
    fields::Number<size_t> m_irow_parent;
    fields::Number<unsigned short> m_iorb;
    fields::Number<unsigned short> m_jorb;
    fields::Number<unsigned short> m_aorb;
    fields::Number<unsigned short> m_borb;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;

    IntermediateTable(size_t nroot, size_t nreplica):
    m_irow_parent(this, "index of the parent row in the WalkerTable"),
    m_iorb(this, "first occupied spin orbital index"),
    m_jorb(this, "second occupied spin orbital index"),
    m_aorb(this, "first vacant spin orbital index"),
    m_borb(this, "second vacant spin orbital index"),
    m_weight(this, "value of the vector element", nroot, nreplica)
    {}


};


#endif //M7_INTERMEDIATETABLE_H
#endif //M7_INTERMEDIATETABLE_H
