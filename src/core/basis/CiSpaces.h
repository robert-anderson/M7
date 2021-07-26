//
// Created by rja on 08/05/2021.
//

#ifndef M7_CISPACES_H
#define M7_CISPACES_H


#include <src/core/enumerator/Enumerator.h>
#include <src/core/table/BufferedFields.h>
#include <src/core/parallel/RankAllocator.h>
#include <src/core/wavefunction/WalkerTable.h>
#include <src/core/wavefunction/Wavefunction.h>
#include <src/core/util/Foreach.h>
#include "Suites.h"

namespace ci_gen {

    typedef std::function<bool(const fields::mbf_t &)> include_fn_t;

    static include_fn_t default_include_fn() {
        return [](const fields::mbf_t &) { return true; };
    }

    static include_fn_t default_include_fn(const RankAllocator<WalkerTableRow> &ra) {
        return [&ra](const fields::mbf_t &mbf) { return mpi::i_am(ra.get_rank(mbf)); };
    }

    static include_fn_t default_include_fn(const Wavefunction &wf) {
        return default_include_fn(wf.m_ra);
    }

    struct Base {
        const size_t m_nsite, m_nelec;
        mutable suite::Mbfs m_mbf_work;
        const include_fn_t m_include_fn;


        Base(size_t nsite, size_t nelec, include_fn_t include_fn = default_include_fn());

        void add_if_included(Row &row, fields::mbf_t &mbf) {
            if (m_include_fn(m_mbf_work[mbf])) {
                row.push_back_jump();
                mbf = m_mbf_work[mbf];
                row.m_table->post_insert(row.index());
            }
        }

    protected:
        static void set_from_inds(fields::FrmOnv& onv, const defs::inds& inds){
            onv = inds;
        }
        static void set_from_inds(fields::FrmBosOnv& onv, const defs::inds& inds){
            onv.m_frm = inds;
        }
        static void set_from_inds(fields::FrmOnv& onv, const defs::inds& alpha_inds, const defs::inds& beta_inds){
            onv.set(alpha_inds, beta_inds);
        }
        static void set_from_inds(fields::FrmBosOnv& onv, const defs::inds& alpha_inds, const defs::inds& beta_inds){
            onv.m_frm.set(alpha_inds, beta_inds);
        }
    };

    struct NoSym : Base {
        foreach::rtnd::Ordered<> m_foreach;

        NoSym(size_t nsite, size_t nelec, const include_fn_t &include_fn = default_include_fn());

    public:

        void operator()(Row &row, fields::mbf_t &mbf){
            ASSERT(mbf.belongs_to_row(&row));
            row.m_table->clear();
            auto body = [&](){
                set_from_inds(m_mbf_work[mbf], m_foreach.inds());
                add_if_included(row, mbf);
            };
            m_foreach(body);
        }
    };

    struct SpinSym : Base {
        foreach::rtnd::Ordered<> m_foreach_alpha;
        foreach::rtnd::Ordered<> m_foreach_beta;

        SpinSym(size_t nsite, size_t nelec, int spin, const include_fn_t &include_fn = default_include_fn());

        void operator()(Row &row, fields::mbf_t &mbf);
    };

};


#endif //M7_CISPACES_H
