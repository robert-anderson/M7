//
// Created by rja on 07/04/2021.
//

#ifndef M7_FRMONVFIELD_H
#define M7_FRMONVFIELD_H

#include <src/core/basis/BasisData.h>
#include <src/core/caches/DecodedMbf.h>
#include "BitsetField.h"

struct FrmOnvField : BitsetField<size_t, 2> {
    typedef BitsetField<size_t, 2> base_t;
    using base_t::get;
    using base_t::set;
    using base_t::clr;
    using base_t::put;
    using base_t::operator=;
    using base_t::inds_t;

    const size_t m_nsite;
    const size_t& m_nspinorb;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmOnv m_decoded;

    FrmOnvField(Row* row, BasisData bd, std::string name="");
    FrmOnvField(Row* row, size_t nsite, std::string name="");

    FrmOnvField(const FrmOnvField& other);

    FrmOnvField& operator=(const FrmOnvField& other) {
        base_t::operator=(other);
        return *this;
    }

    FrmOnvField &operator=(std::pair<const defs::inds &, const defs::inds &> setbits);

    bool operator==(const FrmOnvField& other) const {
        return base_t::operator==(other);
    }

    bool operator==(const defs::inds& inds) const;

    void set(const size_t& bit_offset, const defs::inds& setbits);

    void set(const size_t& site_offset, const defs::inds& setbits_alpha, const defs::inds& setbits_beta);

    void set(const defs::inds& setbits_alpha, const defs::inds& setbits_beta);

    void set_spins(const defs::inds& alpha_sites);

    void put_spin_channel(const size_t& ispin, bool set);

    void clr(const size_t& bit_offset, const defs::inds& clrbits);

    void clr(const size_t& site_offset, const defs::inds& clrbits_alpha, const defs::inds& clrbits_beta);

    void excite(const size_t &i, const size_t &j) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, i);
        set(dptr, j);
    }
    void excite(inds_t ann, inds_t cre){
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, ann);
        set(dptr, cre);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, i);
        clr(dptr, j);
        set(dptr, k);
        set(dptr, l);
    }
    void excite(inds_t ann1, inds_t ann2, inds_t cre1, inds_t cre2) {
        auto* dptr = reinterpret_cast<size_t *>(begin());
        clr(dptr, ann1);
        clr(dptr, ann1);
        set(dptr, cre1);
        set(dptr, cre2);
    }

    void foreach(const std::function<void(const size_t&)>& body_fn) const;

private:
    void foreach_pair_inner(const size_t& ibit, const std::function<void(const size_t&, const size_t&)>& body_fn) const;

public:
    void foreach_pair(const std::function<void(const size_t&, const size_t&)>& body_fn) const;

    void foreach_pair(const std::function<void(const size_t&)>& body_fn_outer,
                      const std::function<void(const size_t&, const size_t&)>& body_fn_inner) const;

    /**
     * 2 x the total z-axis projection of total spin
     * @return
     *  nalpha-nbeta where alpha is conventionally the 0 spin channel
     */
    int ms2() const;

    int nalpha() const;

    size_t site_nocc(const size_t& isite) const;

    bool all_sites_single_occ() const;

    std::string to_string() const override;

    const size_t& nsite() const;


    static size_t isite(const size_t& ibit, const size_t& nsite);

    static size_t ispin(const size_t& ibit, const size_t& nsite);

    static size_t ibit(const size_t& ispin, const size_t& isite, const size_t& nsite);

    /**
     * convenient conversion between spin-site and site indices
     * @param ibit
     *  bit index
     * @return
     *  associated site index
     */
    size_t isite(const size_t& ibit) const;
    /**
     * convenient conversion between spin-site and spin indices
     * @param ibit
     *  bit index
     * @return
     *  associated spin index (0 or 1)
     */
    size_t ispin(const size_t& ibit) const;
    /**
     * convenient conversion between spin and site indices and the flat "bit" index
     * @param ispin
     *  spin channel index
     * @param isite
     *  site index
     * @return
     *  flat index od spin-site
     */
    size_t ibit(const size_t& ispin, const size_t& isite) const;
};


#endif //M7_FRMONVFIELD_H
