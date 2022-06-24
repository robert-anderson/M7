//
// Created by Robert J. Anderson on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"
#include "Lattice.h"

using namespace defs;

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

        defs::info_map_t info() const {
            return {{"value", m_value}, {"conserved", m_conserve}};
        }

        std::string to_string() const {
            return convert::to_string(info());
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
        static constexpr int c_undefined_ms2 = std::numeric_limits<int>::max();
        /**
         * extent of the fermionic single particle basis
         */
        struct Size {
            const uint_t m_nsite;
            const uint_t m_nspinorb;
            /**
             * number of pairs of distinct spin orbitals
             */
            const uint_t m_nspinorb_pair;

            Size(uint_t nsite);

            operator uint_t () const {
                return m_nsite;
            }

            uint_t isite(uint_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb < m_nsite ? ispinorb : ispinorb - m_nsite;
            }
            static uint_t isite(uint_t ispinorb, uint_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb < nsite ? ispinorb : ispinorb - nsite;
            }

            uint_t ispin(uint_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb >= m_nsite;
            }
            static uint_t ispin(uint_t ispinorb, uint_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb >= nsite;
            }

            uint_t ispinorb(uint_t ispin, uint_t isite) const {
                DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
                DEBUG_ASSERT_LT(isite, m_nsite, "site index OOB");
                return ispin ? isite + m_nsite : isite;
            }
            static uint_t ispinorb(uint_t ispin, uint_t isite, uint_t nsite) {
                DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
                DEBUG_ASSERT_LT(isite, nsite, "site index OOB");
                return ispin ? isite + nsite : isite;
            }

            uint_t ispinorb(std::pair<uint_t, uint_t> pair) const {
                return ispinorb(pair.first, pair.second);
            }
            static uint_t ispinorb(std::pair<uint_t, uint_t> pair, uint_t nsite) {
                return ispinorb(pair.first, pair.second, nsite);
            }

            int ms2(uint_t ispinorb) const {
                DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
                return ispinorb < m_nsite ? 1 : -1;
            }
            static int ms2(uint_t ispinorb, uint_t nsite) {
                DEBUG_ASSERT_LT(ispinorb, 2*nsite, "spin orbital index OOB");
                return ispinorb < nsite ? 1 : -1;
            }

            /**
             * @param spin_resolved
             *  true if the basis is spin-resolved e.g. UHF
             * @return
             *  number of indices needed in the access of coefficients
             */
            uint_t ncoeff_ind(bool spin_resolved) const;
            static uint_t ncoeff_ind(bool spin_resolved, uint_t nsite);
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

            Basis(uint_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, std::shared_ptr<lattice::Base> lattice);

            Basis(std::shared_ptr<lattice::Base> lattice);

            Basis(uint_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved);
            /*
             * non-resolved spin, C1 point group (no spatial symmetry)
             */
            Basis(uint_t nsite): Basis(nsite, {nsite}, false) {}

            bool operator==(const Basis& other) const;

            uint_t ncoeff_ind() const;

            defs::info_map_t info() const;

            std::string to_string() const;
        };

        /**
         * 2*Ms value to which the fermionic many-body basis is constrained if conserved, else the value only serves
         * as a hint for the purpose of reference / initial MBF selection
         */
        struct Ms2 : public conservation::Optional<int>{
            /**
             * either 0 or 1 based on the evenness of the electron number
             */
            static int lowest_value(uint_t nelec);
            Ms2(int v, bool conserve=true);
            Ms2();
        };

        struct Electrons {
        private:
            /**
             * number of electrons in the system (0ul if not conserved)
             */
            const uint_t m_n;

        public:
            /**
             * number of pairs of distinct electrons in the system (0ul if not conserved)
             */
            const uint_t m_npair;
            const Ms2 m_ms2;
            /**
             * numbers of occupied alpha and beta electrons (0 if spin unconserved)
             */
            const uint_t m_nalpha, m_nbeta;

            Electrons(uint_t n, Ms2 ms2);

            Electrons(uint_t n): Electrons(n, Ms2(Ms2::lowest_value(n))){}

            Electrons(const Electrons& e1, const Electrons& e2);

            operator uint_t () const {
                return m_n;
            }

            bool operator==(const Electrons& other) const;

            defs::info_map_t info() const;
            std::string to_string() const;
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
            const uint_t m_nvac;
            /**
             * numbers of occupied alpha and beta vacant orbitals (0 if spin unconserved)
             */
            const uint_t m_nvac_alpha, m_nvac_beta;

            explicit Sector(Basis basis, Electrons elecs);

            /*
             * combine the properties of two Hilbert space sectors
             */
            bool operator==(const Sector& other) const;

            operator bool() const {
                return m_basis.m_nsite;
            }

            uint_t size() const;
        };
    }

    namespace bos {
        /**
         * maximum occupation supported by the mode occupation type
         */
        static constexpr uint_t c_max_occ = std::numeric_limits<bos_occ_t>::max();
        /**
         * extent of the bosonic single particle basis
         */
        struct Size {
            const uint_t m_nmode;

            Size(uint_t nmode);

            operator uint_t() const {
                return m_nmode;
            }
        };

        struct Basis : Size {
            const uint_t m_occ_cutoff;
        private:
            using Size::operator unsigned long;
        public:
            explicit operator bool() const {
                return m_nmode;
            }
            Basis(uint_t nmode, uint_t occ_cutoff=c_max_occ);

            bool operator==(const Basis& other) const;
            
            defs::info_map_t info() const;

            std::string to_string() const;
        };

        struct Bosons : public conservation::Optional<uint_t> {
            Bosons(uint_t v, bool conserve);
            Bosons();
            Bosons(const Bosons& b1, const Bosons& b2);
        };

        struct Sector {
            const Basis m_basis;
            const Bosons m_bosons;
            explicit Sector(Basis basis, Bosons bosons);

            bool operator==(const Sector& other) const;

            operator bool() const {
                return bool(m_basis);
            }
        };

    }

    struct Size {
        const frm::Size m_frm;
        const bos::Size m_bos;
        Size(uint_t nsite, uint_t nmode);
        void require_pure_frm() const;
        void require_pure_bos() const;
    };

    struct Basis {
        const frm::Basis m_frm;
        const bos::Basis m_bos;

        Basis(frm::Basis frm, bos::Basis bos);

        void require_pure_frm() const;
        void require_pure_bos() const;

        Size size() const;

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

        Sector(frm::Sector frm, bos::Sector bos);
        explicit Sector(frm::Sector frm);
        explicit Sector(bos::Sector bos);
        Sector(Basis basis, Particles particles);

        Basis basis() const;

        // TODO: rename - "size" makes it sound like the hilbert space dimension is returned
        Size size() const;

        Particles particles() const;
    };
}

/*
static std::ostream &operator<<(std::ostream &os, const sys::frm::Basis &v) {
    os << log::format("[nsite: {}, spin resolved: {}, site irreps: {}]",
                      v.m_nsite, v.m_spin_resolved, convert::to_string(v.m_abgrp_map.m_site_irreps));
    return os;
}

static std::ostream &operator<<(std::ostream &os, const sys::bos::Basis &v) {
    os << log::format("[nmode: {}, max mode occ: {}]", v.m_nmode, v.m_occ_cutoff);
    return os;
}
 */

#endif //M7_BASISDATA_H