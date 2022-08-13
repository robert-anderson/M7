//
// Created by rja on 13/08/22.
//

#ifndef M7_GUIDE_H
#define M7_GUIDE_H

#include "M7_lib/field/Mbf.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"

namespace guide {
    struct Wavefunction {
    protected:
        virtual ham_t frm_overlap(const field::FrmOnv& /*mbf*/) const {
            return 0.0;
        }

        virtual ham_t bos_overlap(const field::BosOnv& /*mbf*/) const {
            return 0.0;
        }

        virtual ham_t frmbos_overlap(const field::FrmBosOnv& /*mbf*/) const {
            return 0.0;
        }

    public:
        ham_t overlap(const field::FrmOnv& mbf) const;
        ham_t overlap(const field::BosOnv& mbf) const;
        ham_t overlap(const field::FrmBosOnv& mbf) const;
    };

    struct EnergyDependent : Wavefunction {
        const Hamiltonian& m_h;
        EnergyDependent(const Hamiltonian& h);
    };

    struct GutzwillerLike : EnergyDependent {
        const double m_fac;
        GutzwillerLike(const Hamiltonian& h, double fac);

    private:
        template<typename mbf_t>
        ham_t get(const mbf_t& mbf) const {
            return std::exp(-m_fac*m_h.get_element(mbf));
        }
    protected:
        ham_t frm_overlap(const FrmOnv& mbf) const override {
            return get(mbf);
        }

        ham_t bos_overlap(const BosOnv& mbf) const override {
            return get(mbf);
        }

        ham_t frmbos_overlap(const FrmBosOnv& mbf) const override {
            return get(mbf);
        }
    };

    /**
     * exponentially suppress multiple occupancies
     */
    struct SuppressMultiOcc : Wavefunction {
        const double m_fac;
        SuppressMultiOcc(double fac): m_fac(fac) {}

    protected:
        ham_t frm_overlap(const FrmOnv& mbf) const override {
            return std::exp(-m_fac*mbf.nclosed_shell());
        }

        ham_t bos_overlap(const BosOnv& mbf) const override {
            return std::exp(-m_fac*mbf.occ_npair());
        }

        ham_t frmbos_overlap(const FrmBosOnv& mbf) const override {
            return SuppressMultiOcc::frm_overlap(mbf.m_frm)*SuppressMultiOcc::bos_overlap(mbf.m_bos);
        }
    };
}


#endif //M7_GUIDE_H
