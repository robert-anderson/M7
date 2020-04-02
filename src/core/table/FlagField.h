//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_FLAGFIELD_H
#define SANDBOX2_FLAGFIELD_H


#include "src/defs.h"
#include "Flag.h"
#include "BitsetField.h"

class FlagField : public BitsetField{
    std::vector<Flag*> m_flags;
public:
    FlagField(Table* table, size_t nelement): BitsetField(table, nelement, 0){}

    defs::pair add_flag(Flag* flag);

    friend class Flag;
};


#endif //SANDBOX2_FLAGFIELD_H
