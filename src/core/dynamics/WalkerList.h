//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H


#include <src/core/list/PerforableMappedList.h>
#include <src/core/table/DeterminantField.h>
#include <src/core/table/FlagField.h>

struct WalkerList : public PerforableMappedList<DeterminantElement> {
    DeterminantField determinant;
    NumericField<defs::wf_t> weight;
    NumericField<defs::ham_comp_t> hdiag;

private:
    struct Flags : public FlagField {
        Flag reference_connection, initiator, deterministic;

        Flags(Table *table, size_t nelement = 1) :
            FlagField(table, nelement),
            reference_connection(this, nelement),
            initiator(this, nelement),
            deterministic(this, nelement) {}
    };

public:
    Flags flags;

    WalkerList(size_t nsite, size_t nbucket) :
        PerforableMappedList<DeterminantElement>(determinant, nbucket),
        determinant(this, 1, nsite), weight(this), hdiag(this), flags(this) {}
};


#endif //M7_WALKERLIST_H
