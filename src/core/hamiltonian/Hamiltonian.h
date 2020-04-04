//
// Created by rja on 27/02/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <cstddef>
#include <src/core/fermion/Determinant.h>
#include <src/core/list/MappedList.h>
#include <src/core/fermion/Connection.h>
#include "src/consts.h"
#include "src/defs.h"
#include "src/core/table/DeterminantField.h"

class Hamiltonian {
protected:
    const size_t m_nsite;
    size_t m_nci;

public:
    Hamiltonian(const size_t &nsite);

    virtual consts::component_t<defs::ham_t>::type get_energy(const DeterminantElement &det) const = 0;

    virtual defs::ham_t get_element_0(const DeterminantElement &det) const = 0;

    virtual defs::ham_t get_element_1(const DeterminantElement &bra, const size_t &removed, const size_t &inserted) const = 0;

    virtual defs::ham_t get_element_2(const size_t &removed1, const size_t &removed2,
                              const size_t &inserted1, const size_t &inserted2) const = 0;

    virtual defs::ham_t get_element(const DeterminantElement &ket, const AntisymConnection &connection) const = 0;

    defs::ham_t get_element(const DeterminantElement &bra, const DeterminantElement &ket) const {
        return get_element(ket, AntisymConnection(ket, bra));
    }

    size_t nsite() const {
        return m_nsite;
    }

    size_t nci() const {
        return m_nci;
    }

    virtual bool spin_conserving() const = 0;

    virtual size_t nelec() const = 0;

    virtual bool spin_resolved() const = 0;

    Determinant guess_reference(const int &spin_level) const;

    Determinant refine_guess_reference(const DeterminantElement& ref) const;

    Determinant choose_reference(const int &spin_level) const;

    class ConnectionList : public MappedList<DeterminantElement> {
    public:
        DeterminantField determinant;
        NumericField<defs::ham_t> helement;

        ConnectionList(size_t nsite, size_t nbucket) :
            MappedList(determinant, nbucket),
            determinant(this, 1, nsite), helement(this) {}
    };
    ConnectionList all_connections_of_det(const Determinant &ref, const defs::ham_comp_t eps=0.0) const;

};

#endif //M7_HAMILTONIAN_H
