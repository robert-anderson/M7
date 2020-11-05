//
// Created by rja on 27/02/2020.
//

#ifndef M7_FERMIONHAMILTONIAN_H
#define M7_FERMIONHAMILTONIAN_H

#include <cstddef>
#include <src/core/basis/DeterminantConnection.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/field/Elements.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"


class FermionHamiltonian {
protected:
    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

public:
    FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                       bool complex_valued, bool spin_resolved);

    FermionHamiltonian(const FcidumpFileReader<defs::ham_t> &file_reader);


    FermionHamiltonian(const std::string& fname, bool spin_major);

    consts::component_t<defs::ham_t>::type get_energy(const views::FermionOnv &det) const;

    defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const;

    defs::ham_t get_element_0(const OccupiedOrbitals &occs) const;

    defs::ham_t get_element_0(const views::FermionOnv &det) const;

    defs::ham_t get_element_0(const AntisymConnection &connection) const;

    defs::ham_t get_element_1(const AntisymConnection &connection) const;

    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;

    defs::ham_t get_element_2(const DeterminantConnection &connection) const;

    defs::ham_t get_element_2(const AntisymConnection &connection) const;

    defs::ham_t get_element(const AntisymConnection &connection) const;

    defs::ham_t get_element(const views::FermionOnv &bra, const views::FermionOnv &ket) const;

    size_t nci() const {
        return ci_utils::fermion_dim(nsite(), nelec());
    }

    const size_t &nsite() const {
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

    const size_t &nelec() const {
        return m_nelec;
    }

    const bool &complex_valued() const {
        return m_complex_valued;
    }

//    elements::FermionOnv guess_reference(const int &spin_level) const;
//
//    elements::FermionOnv refine_guess_reference(const views::FermionOnv &ref) const;
//
//    elements::FermionOnv choose_reference(const int &spin_level) const;
//
//    class DeterminantList : public MappedList<DeterminantElement> {
//    public:
//        DeterminantField determinant;
//
//        DeterminantList(std::string name, size_t nsite, size_t nbucket) :
//                MappedList(name, determinant, nbucket),
//                determinant(this, 1, nsite){}
//    };
//
//    void generate_ci_space(WalkerList* list, RankAllocator<DeterminantElement>& ra, const int &spin_level) const;

};

#endif //M7_FERMIONHAMILTONIAN_H
