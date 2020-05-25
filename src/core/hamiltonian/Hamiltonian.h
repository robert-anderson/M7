//
// Created by rja on 27/02/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <cstddef>
#include <src/core/fermion/Determinant.h>
#include <src/core/list/MappedList.h>
#include <src/core/fermion/Connection.h>
#include <src/core/fermion/DecodedDeterminant.h>
#include "src/core/util/consts.h"
#include "src/core/util/defs.h"
#include "src/core/table/DeterminantField.h"

class Hamiltonian {
protected:
    const size_t m_nsite;
    size_t m_nci;

public:
    Hamiltonian(const size_t &nsite);

    consts::component_t<defs::ham_t>::type get_energy(const DeterminantElement &det) const {
        return consts::real(get_element_0(det));
    }

    virtual ~Hamiltonian()= default;

    virtual defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const = 0;

    defs::ham_t get_element_0(const OccupiedOrbitals &occs) const {
        return get_element_0(occs.m_inds, occs.m_nind);
    }

    defs::ham_t get_element_0(const DeterminantElement &det) const {
        OccupiedOrbitals occs(det);
        return get_element_0(occs.m_inds, occs.m_nind);
    }

    defs::ham_t get_element_0(const AntisymConnection &connection) const {
        return get_element_0(connection.com(), connection.ncom());
    }

    virtual defs::ham_t get_element_1(const AntisymConnection &connection) const = 0;

    virtual defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const = 0;

    defs::ham_t get_element_2(const Connection &connection) const {
        return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    }

    defs::ham_t get_element_2(const AntisymConnection &connection) const {
        const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
        return connection.phase() ? -element : element;
    }

    defs::ham_t get_element(const AntisymConnection &connection) const {
        switch (connection.nexcit()) {
            case 0:
                return get_element_0(connection);
            case 1:
                return get_element_1(connection);
            case 2:
                return get_element_2(connection);
            default:
                return 0;
        }
    }

    defs::ham_t get_element(const DeterminantElement &bra, const DeterminantElement &ket) const {
        return get_element(AntisymConnection(ket, bra));
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

    Determinant refine_guess_reference(const DeterminantElement &ref) const;

    Determinant choose_reference(const int &spin_level) const;

    class ConnectionList : public MappedList<DeterminantElement> {
    public:
        DeterminantField determinant;
        NumericField<defs::ham_t> helement;

        ConnectionList(size_t nsite, size_t nbucket) :
            MappedList(determinant, nbucket),
            determinant(this, 1, nsite), helement(this) {}
    };

    ConnectionList all_connections_of_det(const Determinant &ref, const defs::ham_comp_t eps = 1e-1) const;

};

#endif //M7_HAMILTONIAN_H
