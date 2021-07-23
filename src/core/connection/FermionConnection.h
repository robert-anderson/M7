//
// Created by rja on 23/07/2021.
//

#ifndef M7_FERMIONCONNECTION_H
#define M7_FERMIONCONNECTION_H

#include "src/core/field/Fields.h"

class FrmOnvConnection {
protected:
    defs::inds m_ann, m_cre;

public:
    explicit FrmOnvConnection(size_t nsite) {
        m_ann.reserve(2 * nsite);
        m_cre.reserve(2 * nsite);
    }

    virtual operator bool() const {
        return nexcit();
    }

    const defs::inds &ann() const { return m_ann; }

    const size_t &ann(const size_t &i) const {
        ASSERT(i < nann());
        return m_ann[i];
    }

    size_t nann() const { return m_ann.size(); }

    const defs::inds &cre() const { return m_cre; }

    const size_t &cre(const size_t &i) const {
        ASSERT(i < ncre());
        return m_cre[i];
    }

    size_t ncre() const { return m_cre.size(); }

    void connect(const fields::Onv<0> &in, const fields::Onv<0> &out) {
        DEBUG_ASSERT_FALSE(in.is_zero(), "should not be computing connection from zero ONV");
        DEBUG_ASSERT_FALSE(out.is_zero(), "should not be computing connection to zero ONV");
        zero();

        size_t in_work, out_work, work;
        for (size_t idataword = 0ul; idataword < in.m_dsize; ++idataword) {
            in_work = in.get_dataword(idataword);
            out_work = out.get_dataword(idataword);
            work = in_work & ~out_work;
            while (work) add_ann(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
            work = out_work & ~in_work;
            while (work) add_cre(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
        }
    }

    void apply(const fields::Onv<0> &in, fields::Onv<0> &out) {
        DEBUG_ASSERT_FALSE(in.is_zero(), "should not be computing connection from zero ONV");
#ifndef NDEBUG
        for (size_t i = 0ul; i < nann(); ++i) ASSERT(in.get(ann(i)));
        for (size_t i = 0ul; i < ncre(); ++i) ASSERT(!in.get(cre(i)));
#endif
        out = in;
        for (size_t i = 0ul; i < nann(); ++i) out.clr(ann(i));
        for (size_t i = 0ul; i < ncre(); ++i) out.set(cre(i));
        DEBUG_ASSERT_EQ(in.nsetbit(), out.nsetbit(),
                        "currently, all excitations are particle number conserving, so something has gone wrong here");
    }

    void zero() {
        m_cre.clear();
        m_ann.clear();
    }

    void add_cre(const size_t &i) {
        ASSERT(m_cre.size() < m_cre.capacity());
        m_cre.push_back(i);
    }

    void add_ann(const size_t &i) {
        ASSERT(m_ann.size() < m_ann.capacity());
        m_ann.push_back(i);
    }

    void add(const size_t &ann, const size_t &cre) {
        add_ann(ann);
        add_cre(cre);
    }

    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2) {
        add_ann(ann1);
        add_ann(ann2);
        add_cre(cre1);
        add_cre(cre2);
    }

    void sort() {
        std::sort(m_cre.begin(), m_cre.end());
        std::sort(m_ann.begin(), m_ann.end());
    }

    size_t nexcit() const;

    bool connected() const {
        return nexcit() <= 2;
    }
};


#endif //M7_FERMIONCONNECTION_H
