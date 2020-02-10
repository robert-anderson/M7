//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_DATASYSTEM_H
#define M7_DATASYSTEM_H


#include "src/data/Table.h"
#include "../parallel/MPIWrapper.h"

class DataSystem {
    Table &m_store;
    Table &m_send_buffer;
    Table &m_recv_buffer;
public:
    DataSystem(Table&, Table&, Table&);
    void communicate();
};


#endif //M7_DATASYSTEM_H
