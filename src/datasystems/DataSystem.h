//
// Created by Robert John Anderson on 2020-02-07.
//

#ifndef M7_DATASYSTEM_H
#define M7_DATASYSTEM_H


#include "../data/DataTable.h"
#include "../parallel/MPIWrapper.h"

class DataSystem {
    DataTable &m_store;
    DataTable &m_send_buffer;
    DataTable &m_recv_buffer;
public:
    DataSystem(DataTable&, DataTable&, DataTable&);
    void communicate();
};


#endif //M7_DATASYSTEM_H
