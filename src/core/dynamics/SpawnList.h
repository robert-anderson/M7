//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNLIST_H
#define M7_SPAWNLIST_H

#include <src/core/table/FlagField.h>
#include <src/core/table/DeterminantField.h>
#include <src/core/table/NumericField.h>
#include <src/core/list/PerforableMappedList.h>


struct SpawnList : public PerforableMappedList<DeterminantElement> {
    DeterminantField determinant;
    NumericField <defs::wf_t> weight;
private:
    struct Flags : public FlagField {
        Flag parent_initiator;
        Flags(Table *table, size_t nelement = 1) :
            FlagField(table, nelement),
            parent_initiator(this, nelement){}
    };

public:
    Flags flags;

    SpawnList(size_t nsite, size_t nbucket) :
        PerforableMappedList<DeterminantElement>(determinant, nbucket),
        determinant(this, 1, nsite), weight(this), flags(this) {}
};


#endif //M7_SPAWNLIST_H
