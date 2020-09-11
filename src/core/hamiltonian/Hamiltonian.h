//
// Created by rja on 27/02/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <cstddef>
#include "src/core/basis/Determinant.h"
#include "src/core/list/MappedList.h"
#include "src/core/dynamics/WalkerList.h"
#include "src/core/basis/Connection.h"
#include "src/core/basis/DecodedDeterminant.h"
#include "src/core/util/consts.h"
#include "src/core/util/defs.h"
#include "src/core/basis/DeterminantField.h"
#include <src/core/parallel/RankAllocator.h>
#include <src/core/io/FcidumpFileReader.h>

class Hamiltonian {
protected:
    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;

public:
    Hamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e, bool complex_valued);

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
                ASSERT(connection.ncom()+connection.nexcit()==nelec());
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

    size_t nci() const {
        return integer_utils::combinatorial(2*nsite(), nelec());
    }

    const size_t& nsite() const {
        return m_nsite;
    }

    bool spin_conserving_1e() const {
        return m_spin_conserving_1e;
    }

    bool spin_conserving_2e() const {
        return m_spin_conserving_2e;
    }

    bool spin_conserving() const {
        return m_spin_conserving_1e && m_spin_conserving_2e;
    }

    const size_t& nelec() const {
        return m_nelec;
    }

    Determinant guess_reference(const int &spin_level) const;

    Determinant refine_guess_reference(const DeterminantElement &ref) const;

    Determinant choose_reference(const int &spin_level) const;

    class DeterminantList : public MappedList<DeterminantElement> {
    public:
        DeterminantField determinant;

        DeterminantList(std::string name, size_t nsite, size_t nbucket) :
                MappedList(name, determinant, nbucket),
                determinant(this, 1, nsite){}
    };

    void generate_ci_space(WalkerList* list, RankAllocator<DeterminantElement>& ra, const int &spin_level) const;

    const bool& complex_valued() const {
        return m_complex_valued;
    }
};

#endif //M7_HAMILTONIAN_H
