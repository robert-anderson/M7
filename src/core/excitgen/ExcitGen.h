//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITGEN_H
#define M7_EXCITGEN_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/connection/Connections.h>
#include "src/core/sample/PRNG.h"

class ExcitGen {
protected:
    const Hamiltonian &m_h;
    PRNG &m_prng;
    const size_t m_nspinorb;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
public:
    ExcitGen(const Hamiltonian &h, PRNG &prng);

    virtual bool draw(const fields::FrmOnv &src_onv,
                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    virtual bool draw(const fields::FrmBosOnv &src_onv,
                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn);

    virtual size_t approx_nconn() const {
        return 0ul;
    }

    virtual std::string description() const {
        return "";
    }

    virtual ~ExcitGen(){}
};

/**
 * Base class for stochastic Fermion sector excitations
 */
class FrmExcitGen : public ExcitGen {
protected:
    const size_t m_nexcit;
    const bool m_spin_conserving;

public:
    FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t nexcit);

    bool draw(const FrmOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmBosOnv &conn) override;

    std::string description() const override {
        return log::format("Fermion annihilate {}, create {}", m_nexcit, m_nexcit);
    }
};

/**
 * Base class for stochastic excitations which may couple Fermion and Boson sectors excitations
 */
class FrmBosExcitGen : public ExcitGen {
public:
    FrmBosExcitGen(const Hamiltonian &h, PRNG &prng) : ExcitGen(h, prng){}

    bool draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return ExcitGen::draw(src_onv, occs, vacs, prob, helem, conn);
    }

    std::string description() const override {
        return "Fermion-Boson coupling";
    }
};

#endif //M7_EXCITGEN_H
