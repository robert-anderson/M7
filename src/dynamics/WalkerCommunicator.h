//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERCOMMUNICATOR_H
#define M7_WALKERCOMMUNICATOR_H

#include <memory>
#include "SpawnList.h"
#include "src/data/List.h"
#include "src/data/TableCommunicator.h"


class WalkerCommunicator : public TableCommunicator<SpawnList> {
public:
    WalkerCommunicator(size_t nsite, size_t nrow_send, size_t nrow_recv);
};


#endif //M7_WALKERCOMMUNICATOR_H
