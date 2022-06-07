//
// Created by Robert J. Anderson on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"
#include "Lattice.h"

namespace conservation {
    enum State {Undefined, Yes, No, Hint};
    /**
     * wrapper class to handle physical values that are allowed to be in the following states:
     *  - undefined
     *  - conserved with a defined value
     *  - unconserved with a defined value (value can be used as a hint when "slight" symmetry breaking occurs)
     *  - unconserved without a defined value
     */
    template<typename T>
    class Optional {
        const T m_value;
        const std::string m_name;
        const bool m_conserve;
    public:
        Optional(T v, bool conserve, std::string name): m_value(v), m_name(std::move(name)), m_conserve(conserve){}
        /*
         * set an undefined value (but state is defined if conserve is false)
         */
        explicit Optional(bool conserve, std::string name): Optional(std::numeric_limits<T>::max(), conserve, name){}
        /*
         * set an undefined state
         */
        explicit Optional(std::string name): Optional(true, name){}
        explicit Optional(const Optional& o1, const Optional& o2) : Optional(o1.defined() ? o1 : o2){
            if (o1.defined() && o2.defined()) {
                REQUIRE_EQ(o1, o2, "incompatible "+m_name+" values");
                REQUIRE_EQ(o1.m_conserve, o2.m_conserve, "incompatible "+m_name+" conservation statuses");
            }
        }

        State state() const {
            if (m_value==std::numeric_limits<T>::max()) return m_conserve ? Undefined : No;
            else return m_conserve ? Yes : Hint;
        }
        bool has_value() const {
            auto tmp = state();
            return tmp==Yes || tmp==Hint;
        }
        bool defined() const {
            return state()!=Undefined;
        }

        operator T() const {
            // value accesses should be guarded by checking if the quantity is defined first
            DEBUG_ASSERT_TRUE(has_value(), "accessing undefined "+m_name+" value");
            return m_value;
        }
        bool conserve() const {
            DEBUG_ASSERT_NE(state(), Undefined, "accessing "+m_name+" in undefined state");
            return m_conserve;
        }
        bool operator==(const Optional& other) const {
            return m_value==other.m_value && m_conserve==other.m_conserve;
        }
    };
};

namespace sys {
    /**
     * optionally conservable physical values are allowed to be in the following states:
     *  - undefined
     *  - conserved with a defined value
     *  - unconserved without a defined value
     *  - unconserved with a defined value (value can be used as a hint when "slight" symmetry breaking occurs)
     */


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
            static size_t isite(size_t ispinorb, size_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb < nsite ? ispinorb : ispinorb - nsite;
            }

            size_t ispin(size_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb >= m_nsite;
            }
            static size_t ispin(size_t ispinorb, size_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb >= nsite;
            }

            size_t ispinorb(size_t ispin, size_t isite) const {
                DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
                DEBUG_ASSERT_LT(isite, m_nsite, "site index OOB");
                return ispin ? isite + m_nsite : isite;
            }
            static size_t ispinorb(size_t ispin, size_t isite, size_t nsite) {
                DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
                DEBUG_ASSERT_LT(isite, nsite, "site index OOB");
                return ispin ? isite + nsite : isite;
            }

            size_t ispinorb(std::pair<size_t, size_t> pair) const {
                return ispinorb(pair.first, pair.second);
            }
            static size_t ispinorb(std::pair<size_t, size_t> pair, size_t nsite) {
                return ispinorb(pair.first, pair.second, nsite);
            }

            int ms2(size_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb < m_nsite ? 1 : -1;
            }
            static int ms2(size_t ispinorb, size_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb < nsite ? 1 : -1;
            }

            /**
             * @param spin_resolved
             *  true if the basis is spin-resolved e.g. UHF
             * @return
             *  number of indices needed in the access of coefficients
             */
            size_t ncoeff_ind(bool spin_resolved) const {
                return spin_resolved ? m_nspinorb : m_nsite;
            }
            static size_t ncoeff_ind(bool spin_resolved, size_t nsite) {
                return spin_resolved ? 2*nsite : nsite;
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
            const bool m_spin_resolved;
            /**
             * set to non-null if the basis has lattice structure
             */
            const std::shared_ptr<lattice::Base> m_lattice;

        private:
            using Size::operator unsigned long;

        public:

            explicit operator bool() const {
                return m_nsite;
            }

            Basis(size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, std::shared_ptr<lattice::Base> lattice):
                    Size(nsite), m_abgrp_map(std::move(abgrp_map)), m_spin_resolved(spin_resolved), m_lattice(lattice){}

            Basis(std::shared_ptr<lattice::Base> lattice):
                Basis(lattice->m_nsite, {lattice->m_nsite}, false, lattice){}

            Basis(size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved):
                Basis(nsite, abgrp_map, spin_resolved, nullptr){}
            /*
             * non-resolved spin, C1 point group (no spatial symmetry)
             */
            Basis(size_t nsite): Basis(nsite, {nsite}, false) {}

            bool operator==(const Basis& other) const {
                return (m_nsite == other.m_nsite) && (m_abgrp_map == other.m_abgrp_map) && (m_spin_resolved == other.m_spin_resolved);
            }

            size_t ncoeff_ind() const {
                return Size::ncoeff_ind(m_spin_resolved);
            }
        };

        /**
         * 2*Ms value to which the fermionic many-body basis is constrained if conserved, else the value only serves
         * as a hint for the purpose of reference / initial MBF selection
         */
        struct Ms2 : public conservation::Optional<int>{
            /**
             * either 0 or 1 based on the evenness of the electron number
             */
            static int lowest_value(size_t nelec){
                return nelec&1ul;
            }
            Ms2(int v, bool conserve=true): conservation::Optional<int>(v, conserve, "2*Ms"){}
            Ms2(): conservation::Optional<int>("2*Ms"){}
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
            const Ms2 m_ms2;
            /**
             * numbers of occupied alpha and beta electrons (0 if spin unconserved)
             */
            const size_t m_nalpha, m_nbeta;

            Electrons(size_t n, Ms2 ms2): m_n(n), m_npair(integer_utils::nspair(m_n)), m_ms2(ms2),
                                          m_nalpha(m_ms2.conserve() ? (m_n+m_ms2)/2 : 0ul),
                                          m_nbeta(m_ms2.conserve() ? m_n-m_nalpha : 0ul) {
                if (m_ms2.conserve() && m_n)
                    REQUIRE_EQ(size_t(std::abs(m_ms2) % 2), m_n % 2,
                               "2*Ms quantum number given incompatible with number of electrons");
            }

            Electrons(size_t n): Electrons(n, Ms2(Ms2::lowest_value(n))){}

            Electrons(const Electrons& e1, const Electrons& e2):
                    Electrons(e1.m_n ? e1.m_n : e2.m_n, Ms2(e1.m_ms2, e2.m_ms2)){
                if (e1 && e2) REQUIRE_EQ(e1, e2, "incompatible numbers of electrons");
            }

            operator size_t () const {
                return m_n;
            }

            bool operator==(const Electrons& other){
                return m_n==other.m_n && m_ms2==other.m_ms2;
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
                m_nvac_beta(m_elecs.m_ms2.conserve() ? m_basis.m_nsite - m_elecs.m_nbeta : 0ul){
                REQUIRE_LE_ALL(m_elecs, m_basis.m_nspinorb, "unphysical number of electrons");
            }

            /*
             * combine the properties of two Hilbert space sectors
             */
            bool operator==(const Sector& other) const {
                return m_basis==other.m_basis && m_elecs==other.m_elecs;
            }

            operator bool() const {
                return m_basis.m_nsite;
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
            using Size::operator unsigned long;
        public:
            explicit operator bool() const {
                return m_nmode;
            }
            Basis(size_t nmode, size_t occ_cutoff=defs::max_bos_occ): Size(nmode), m_occ_cutoff(occ_cutoff){}

            bool operator==(const Basis& other) const {
                return (m_occ_cutoff == other.m_occ_cutoff) && (m_nmode == other.m_nmode);
            }
        };

        struct Bosons : public conservation::Optional<size_t> {
            Bosons(size_t v, bool conserve): conservation::Optional<size_t>(v, conserve, "nboson"){}
            Bosons(): Bosons(~0ul, true){}
            Bosons(const Bosons& b1, const Bosons& b2): Bosons(b1.defined() ? b1 : b2){}
        };

        struct Sector {
            const Basis m_basis;
            const Bosons m_bosons;
            explicit Sector(Basis basis, Bosons bosons): m_basis(basis), m_bosons(bosons){}

            bool operator==(const Sector& other) const {
                return m_basis==other.m_basis && m_bosons==other.m_bosons;
            }

            operator bool() const {
                return bool(m_basis);
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
            return {static_cast<const frm::Size&>(m_frm), static_cast<const bos::Size&>(m_bos)};
        }

        operator bool() const {
            return bool(m_frm) || bool(m_bos);
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

        // TODO: rename - "size" makes it sound like the hilbert space dimension is returned
        Size size() const {
            return basis().size();
        }

        Particles particles() const {
            return {m_frm.m_elecs, m_bos.m_bosons};
        }
    };
}

static std::ostream &operator<<(std::ostream &os, const sys::frm::Basis &v) {
    os << log::format("[nsite: {}, spin resolved: {}, site irreps: {}]",
                      v.m_nsite, v.m_spin_resolved, utils::to_string(v.m_abgrp_map.m_site_irreps));
    return os;
}

static std::ostream &operator<<(std::ostream &os, const sys::bos::Basis &v) {
    os << log::format("[nmode: {}, max mode occ: {}]", v.m_nmode, v.m_occ_cutoff);
    return os;
}

#endif //M7_BASISDATA_H