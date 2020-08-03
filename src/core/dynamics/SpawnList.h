//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H

#include <src/core/table/FlagField.h>
#include <src/core/table/DeterminantField.h>
#include <src/core/table/NumericField.h>
#include <src/core/list/List.h>


struct SpawnList : public List {
    DeterminantField m_determinant;
    NumericField <defs::wf_t> m_weight;
private:
    struct Flags : public FlagField {
        Flag m_parent_initiator;
        Flag m_parent_deterministic;
        Flags(Table *table, size_t nelement = 1, const std::string &description="") :
            FlagField(table, nelement, description),
            m_parent_initiator(this, nelement),
            m_parent_deterministic(this, nelement)
            {}
    };

public:
    Flags m_flags;

    SpawnList(std::string name, size_t nsite, size_t nsegment) : List(name, nsegment),
        m_determinant(this, 1, nsite, "Determinant"),
        m_weight(this, 1, "Weight"),
        m_flags(this, 1, "Flags") {}

    size_t add(const size_t isegment, const DeterminantElement &determinant, const defs::wf_t &weight,
               bool flag_parent_initiator, bool flag_parent_deterministic){
        auto irow = List::push(isegment);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << "reserved row of segment: " << irow << std::endl;
#endif
        m_determinant(irow, isegment) = determinant;
        m_weight(irow, isegment) = weight;
        m_flags.m_parent_initiator(irow, isegment) = flag_parent_initiator;
        m_flags.m_parent_deterministic(irow, isegment) = flag_parent_deterministic;
        return irow;
    }
};


#endif //M7_SPAWNLIST_H
