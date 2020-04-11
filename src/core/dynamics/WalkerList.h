//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H


#include <src/core/list/PerforableMappedList.h>
#include <src/core/table/DeterminantField.h>
#include <src/core/table/FlagField.h>

struct WalkerList : public PerforableMappedList<DeterminantElement> {
    DeterminantField m_determinant;
    NumericField<defs::wf_t> m_weight;
    NumericField<defs::ham_comp_t> m_hdiag;

private:
    struct Flags : public FlagField {
        Flag m_initiator, m_reference_connection, m_deterministic;

        Flags(Table *table, size_t nelement = 1, const std::string &description="") :
            FlagField(table, nelement, description),
            m_initiator(this, nelement, "Initiator status"),
            m_reference_connection(this, nelement, "Connected to reference determinant"),
            m_deterministic(this, nelement, "In deterministic subspace") {}
    };

public:
    Flags m_flags;

    WalkerList(size_t nsite, size_t nbucket) :
        PerforableMappedList<DeterminantElement>(m_determinant, nbucket),
        m_determinant(this, 1, nsite, "Determinant"),
        m_weight(this, 1, "Weight"),
        m_hdiag(this, 1, "Diagonal Hamiltonian matrix element"),
        m_flags(this, 1, "Flags") {}

    using PerforableMappedList<DeterminantElement>::push;

    size_t add(Mutex &mutex, const DeterminantElement &key, const defs::wf_t &weight, const defs::ham_comp_t &hdiag,
               bool initiator = false, bool reference_connection = false, bool deterministic = false) {
        auto irow = PerforableMappedList::push(mutex, key);
        m_weight(irow) = weight;
        m_hdiag(irow) = hdiag;
        m_flags.m_initiator(irow) = initiator;
        assert(m_flags.m_initiator(irow)==true);
        m_flags.m_reference_connection(irow) = reference_connection;
        m_flags.m_deterministic(irow) = deterministic;
        return irow;
    }

    size_t add(const DeterminantElement &key, const defs::wf_t &weight, const defs::ham_comp_t &hdiag,
               bool initiator = false, bool reference_connection = false, bool deterministic = false) {
        auto mutex = key_mutex(key);
        return add(mutex, key, weight, hdiag, initiator, reference_connection, deterministic);
    }
};


#endif //M7_WALKERLIST_H
