//
// Created by rja on 07/04/2021.
//

#ifndef M7_FRMONVFIELD_H
#define M7_FRMONVFIELD_H

#include <src/core/basis/BasisDims.h>
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

    FrmOnvField(Row* row, size_t nsite, std::string name="");

    FrmOnvField(Row* row, BasisDims bd, std::string name="");

    FrmOnvField(const FrmOnvField& other);

    FrmOnvField& operator=(const FrmOnvField& other) {
        FieldBase::operator=(other);
        return *this;
    }

    FrmOnvField &operator=(std::pair<const defs::inds &, const defs::inds &> setbits);

    void set(const size_t& bit_offset, const defs::inds& setbits) {
        for(auto& i: setbits) set(bit_offset+i);
    }

    void set(const size_t& site_offset, const defs::inds& setbits_alpha, const defs::inds& setbits_beta) {
        for(auto& i: setbits_alpha) set({0, site_offset+i});
        for(auto& i: setbits_beta) set({1, site_offset+i});
    }

    void set(const defs::inds& setbits_alpha, const defs::inds& setbits_beta) {
        zero();
        for(auto& i: setbits_alpha) set({0, i});
        for(auto& i: setbits_beta) set({1, i});
    }

    void clr(const size_t& bit_offset, const defs::inds& clrbits) {
        for(auto& i: clrbits) clr(bit_offset+i);
    }

    void clr(const size_t& site_offset, const defs::inds& clrbits_alpha, const defs::inds& clrbits_beta) {
        for(auto& i: clrbits_alpha) clr({0, site_offset+i});
        for(auto& i: clrbits_beta) clr({1, site_offset+i});
    }

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

    std::string to_string() const override;

    const size_t& nsite() const;

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

    static size_t isite(const size_t& ibit, const size_t& nsite);

    static size_t ispin(const size_t& ibit, const size_t& nsite);
};


#endif //M7_FRMONVFIELD_H
