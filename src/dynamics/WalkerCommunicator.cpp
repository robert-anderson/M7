//
// Created by Robert John Anderson on 2020-02-24.
//

#include "WalkerCommunicator.h"

WalkerCommunicator::WalkerCommunicator(size_t nsite, size_t nrow_send, size_t nrow_recv) :
        TableCommunicator<SpawnList>(SpawnList::spec_T(nsite), nrow_send, nrow_recv){}