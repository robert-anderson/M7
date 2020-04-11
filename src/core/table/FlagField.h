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
    FlagField(Table* table, size_t nelement, const std::string& description):
    BitsetField(table, nelement, 0, description){}

    defs::pair add_flag(Flag* flag);

    const std::string description() const override {
        std::string result = m_description;
        size_t i=0ul;
        for (auto flag : m_flags) {
            result+="\n\t\t"+std::to_string(i)+": "+flag->m_description;
            ++i;
        }
        return result;
    }

    friend class Flag;
};


#endif //SANDBOX2_FLAGFIELD_H
