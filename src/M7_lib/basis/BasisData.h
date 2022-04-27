//
// Created by Robert J. Anderson on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"


namespace sys {
    /**
     * 2*Ms value to which the fermionic many-body basis is constrained if conserved, else the value only serves
     * as a hint for the purpose of reference / initial MBF selection
     */
    template<typename T>
    class OptionallyConserved {
        const T m_v;
        const std::string m_name;
        const bool m_conserve;
    public:
        OptionallyConserved(T v, bool conserve, std::string name):
            m_v(v), m_name(std::move(name)), m_conserve(conserve){}
        /*
         * set an undefined value
         */
        explicit OptionallyConserved(std::string name): OptionallyConserved(std::numeric_limits<T>::max(), false, name){}
        explicit OptionallyConserved(const OptionallyConserved& oc1, const OptionallyConserved& oc2) :
            OptionallyConserved(oc1.defined() ? oc1 : oc2){
            if (oc1.defined() && oc2.defined())
                REQUIRE_EQ(oc1, oc2, "incompatible "+m_name+" values");
        }
        operator T() const {
            // accesses should be guarded by checking if the quantity is defined first
            DEBUG_ASSERT_TRUE(defined(), "accessing undefined "+m_name+" value");
            return m_v;
        }
        bool conserve() const {
            // accesses should be guarded by checking if the quantity is defined first
            DEBUG_ASSERT_TRUE(defined(), "accessing undefined "+m_name+" value");
            return m_conserve;
        }
        bool defined() const {
            return m_v!=std::numeric_limits<T>::max();
        }
        bool operator==(const OptionallyConserved& other) const {
            return m_v==other.m_v && m_conserve==other.m_conserve;
        }
    };

    namespace frm {
        /**
         * extent of the fermionic single particle basis
         */
        struct Size {
            const size_t m_nsite;
            const size_t m_nspinorb;
            /**
             * number of pairs of distinct spin orbitals
             */
            const size_t m_nspinorb_pair;

            Size(size_t nsite):
                m_nsite(nsite), m_nspinorb(2 * nsite), m_nspinorb_pair(integer_utils::nspair(m_nspinorb)){}

            operator size_t () const {
                return m_nsite;
            }

            size_t isite(size_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb < m_nsite ? ispinorb : ispinorb - m_nsite;
            }

            size_t ispin(size_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb >= m_nsite;
            }

            size_t ispinorb(size_t ispin, size_t isite) const {
                DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
                DEBUG_ASSERT_LT(isite, m_nsite, "site index OOB");
                return ispin ? isite + m_nsite : isite;
            }

            size_t ispinorb(std::pair<size_t, size_t> pair) const {
                return ispinorb(pair.first, pair.second);
            }

            int ms2(size_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb < m_nsite ? 1 : -1;
            }

            /**
             * @param restricted_orbs
             *  true if the basis is spin-resolved e.g. UHF
             * @return
             *  number of indices needed in the access of coefficients
             */
            size_t ncoeff_ind(bool restricted_orbs) const {
                return restricted_orbs ? m_nspinorb : m_nsite;
            }
        };

        /**
         * properties of the fermionic single-particle basis
         */
        struct Basis : Size {
            /**
             * mapping from fermion site indices to Abelian group labels
             */
            const AbelianGroupMap m_abgrp_map;
            /**
             * true if the two spin orbitals corresponding to the same site have identical functional form e.g. in UHF basis
             */
            const bool m_restricted_orbs;

        private:
            operator size_t() const {
                return 0ul;
            }
        public:
            operator bool() const {
                return m_nsite;
            }

            Basis(size_t nsite, AbelianGroupMap abgrp_map, bool restricted_orbs):
                Size(nsite), m_abgrp_map(std::move(abgrp_map)), m_restricted_orbs(restricted_orbs){}
            /*
             * non-resolved spin, C1 point group (no spatial symmetry)
             */
            Basis(size_t nsite): Basis(nsite, {nsite}, false){}

            /*
             * combine the properties of two single-particle bases
             */
            explicit Basis(const Basis& basis1, const Basis& basis2):
                Basis(basis1.m_nsite ? basis1.m_nsite : basis2.m_nsite,
                      basis1.m_abgrp_map ? basis1.m_abgrp_map : basis2.m_abgrp_map,
                      basis1.m_restricted_orbs && basis2.m_restricted_orbs) {
                if (basis1 && basis2) {
                    REQUIRE_EQ(basis1.m_nsite, basis2.m_nsite, "incompatible numbers of sites");
                    if (basis1.m_abgrp_map && basis2.m_abgrp_map) {
                        // both Hilbert spaces define point groups, so check compatibility
                        REQUIRE_EQ(basis1.m_abgrp_map, basis2.m_abgrp_map, "incompatible Abelian group maps");
                    }
                }
            }

            bool operator==(const Basis& other) const {
                return m_nsite==other.m_nsite &&
                        m_abgrp_map==other.m_abgrp_map &&
                        m_restricted_orbs==other.m_restricted_orbs;
            }

            size_t ncoeff_ind() const {
                return Size::ncoeff_ind(m_restricted_orbs);
            }
        };

        struct Electrons {
        private:
            /**
             * number of electrons in the system (0ul if not conserved)
             */
            const size_t m_n;
        public:
            /**
             * number of pairs of distinct electrons in the system (0ul if not conserved)
             */
            const size_t m_npair;
            /**
             * 2*Ms value to which the fermionic many-body basis is constrained if conserved, else the value only serves
             * as a hint for the purpose of reference / initial MBF selection
             */
            struct Ms2 : public OptionallyConserved<int>{
                Ms2(int v, bool conserve): OptionallyConserved<int>(v, conserve, "2*Ms"){}
                /*
                 * set either 0 or 1 based on the evenness of the electron number
                 */
                Ms2(size_t nelec): Ms2(!(nelec&1ul), false){}
                Ms2(): OptionallyConserved<int>("2*Ms"){}
            };
            const Ms2 m_ms2;
            /**
             * numbers of occupied alpha and beta electrons (0 if spin unconserved)
             */
            const size_t m_nalpha, m_nbeta;

            Electrons(size_t n, Ms2 ms2): m_n(n), m_npair(integer_utils::nspair(m_n)), m_ms2(ms2),
                                          m_nalpha(m_ms2.conserve() ? (m_n+m_ms2)/2 : 0ul), m_nbeta(m_n-m_nalpha) {
                if (m_ms2.conserve() && m_n)
                    REQUIRE_EQ(size_t(std::abs(m_ms2) % 2), m_n % 2,
                               "2*Ms quantum number given incompatible with number of electrons");
            }

            Electrons(size_t n): Electrons(n, Ms2(n)){}

            Electrons(const Electrons& e1, const Electrons& e2):
                    Electrons(e1.m_n ? e1.m_n : e2.m_n, Ms2(e1.m_ms2, e2.m_ms2)){
                if (e1 && e2) REQUIRE_EQ(e1, e2, "incompatible numbers of electrons");
            }

            operator size_t () const {
                return m_n;
            }
        };

        /**
         * properties of a fermionic many-body Hilbert space sector
         */
        struct Sector {
            const Basis m_basis;
            const Electrons m_elecs;
            /**
             * number of vacant spin orbitals
             */
            const size_t m_nvac;
            /**
             * numbers of occupied alpha and beta vacant orbitals (0 if spin unconserved)
             */
            const size_t m_nvac_alpha, m_nvac_beta;

            explicit Sector(Basis basis, Electrons elecs):
                m_basis(basis), m_elecs(elecs),
                m_nvac(m_basis.m_nspinorb - m_elecs),
                m_nvac_alpha(m_elecs.m_ms2.conserve() ? m_basis.m_nsite - m_elecs.m_nalpha : 0ul),
                m_nvac_beta(m_elecs.m_ms2.conserve() ? m_basis.m_nsite - m_elecs.m_nbeta : 0ul){}

            /*
             * combine the properties of two Hilbert space sectors
             */
            explicit Sector(const Sector& sector1, const Sector& sector2):
                    Sector(Basis(sector1.m_basis, sector2.m_basis),
                           Electrons(sector1.m_elecs, sector2.m_elecs)){}

            bool operator==(const Sector& other) const {
                return m_basis==other.m_basis && m_elecs==other.m_elecs;
            }

            operator bool() const {
                return m_basis;
            }

            size_t size() const {
                if (m_elecs.m_ms2.conserve()) {
                    const auto na = integer_utils::combinatorial(m_basis.m_nsite, m_elecs.m_nalpha);
                    const auto nb = integer_utils::combinatorial(m_basis.m_nsite, m_elecs.m_nbeta);
                    return na * nb;
                }
                return integer_utils::combinatorial(m_basis.m_nspinorb, m_elecs);
            }
        };
    }

    namespace bos {
        /**
         * extent of the bosonic single particle basis
         */
        struct Size {
            const size_t m_nmode;

            Size(size_t nmode) : m_nmode(nmode){}

            operator size_t() const {
                return m_nmode;
            }
        };

        struct Basis : Size {
            const size_t m_occ_cutoff;
        private:
            operator size_t() const {
                return 0ul;
            }
        public:
            operator bool() const {
                return m_nmode;
            }
            Basis(size_t nmode, size_t occ_cutoff=defs::max_bos_occ): Size(nmode), m_occ_cutoff(occ_cutoff){}

            /*
             * combine the properties of two single-particle bases
             */
            explicit Basis(const Basis& basis1, const Basis& basis2):
                    Basis(basis1.m_nmode ? basis1.m_nmode : basis2.m_nmode, basis1.m_occ_cutoff){
                if (basis1 && basis2) {
                    REQUIRE_EQ(basis1.m_nmode, basis2.m_nmode, "incompatible numbers of modes");
                    REQUIRE_EQ(basis1.m_occ_cutoff, basis2.m_occ_cutoff, "incompatible occupation cutoffs");
                }
            }

            bool operator==(const Basis& other) const {
                return m_nmode==other.m_nmode && m_occ_cutoff==other.m_occ_cutoff;
            }
        };

        struct Bosons : public OptionallyConserved<size_t> {
            Bosons(size_t v, bool conserve): OptionallyConserved<size_t>(v, conserve, "nboson"){}
            Bosons(): OptionallyConserved<size_t>("nboson"){}
        };

        struct Sector {
            const Basis m_basis;
            const Bosons m_bosons;
            explicit Sector(Basis basis, Bosons bosons): m_basis(basis), m_bosons(bosons){}
            /*
             * combine the properties of two Hilbert space sectors
             */
            explicit Sector(const Sector& sector1, const Sector& sector2):
                Sector(Basis(sector1.m_basis, sector2.m_basis),
                       Bosons(sector1.m_bosons, sector2.m_bosons)){}

            bool operator==(const Sector& other) const {
                return m_basis==other.m_basis && m_bosons==other.m_bosons;
            }

            operator bool() const {
                return m_basis;
            }
        };

    }

    struct Size {
        const frm::Size m_frm;
        const bos::Size m_bos;
        Size(size_t nsite, size_t nmode): m_frm(nsite), m_bos(nmode){}
        void require_pure_frm() const {
            REQUIRE_FALSE(m_bos, "Single particle basis specification is not purely fermionic");
        }
        void require_pure_bos() const {
            REQUIRE_FALSE(m_frm, "Single particle basis specification is not purely bosonic");
        }
    };

    struct Basis {
        const frm::Basis m_frm;
        const bos::Basis m_bos;

        Basis(frm::Basis frm, bos::Basis bos): m_frm(std::move(frm)), m_bos(bos){}

        void require_pure_frm() const {
            size().require_pure_frm();
        }
        void require_pure_bos() const {
            size().require_pure_bos();
        }

        Size size() const {
            return {m_frm.m_nsite, m_bos.m_nmode};
        }

        operator bool() const {
            return m_frm || m_bos;
        }
    };

    struct Particles {
        const frm::Electrons m_frm;
        const bos::Bosons m_bos;
    };

    struct Sector {
        const frm::Sector m_frm;
        const bos::Sector m_bos;

        Sector(frm::Sector frm, bos::Sector bos): m_frm(std::move(frm)), m_bos(std::move(bos)){}
        explicit Sector(frm::Sector frm): Sector(frm, bos::Sector({0ul}, {})){}
        explicit Sector(bos::Sector bos): Sector(frm::Sector({0ul}, {0ul}), bos){}
        Sector(Basis basis, Particles particles):
            Sector(frm::Sector(basis.m_frm, particles.m_frm), bos::Sector(basis.m_bos, particles.m_bos)){}

        Basis basis() const {
            return {m_frm.m_basis, m_bos.m_basis};
        }

        Size size() const {
            return basis().size();
        }

        Particles particles() const {
            return {m_frm.m_elecs, m_bos.m_bosons};
        }
    };
}


#endif //M7_BASISDATA_H