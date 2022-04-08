//
// Created by rja on 08/05/2021.
//

#ifndef M7_CISPACES_H
#define M7_CISPACES_H


#include <M7_lib/enumerator/Enumerator.h>
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/wavefunction/WalkerTable.h>
#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/foreach/Foreach.h>

#include "Suites.h"

namespace ci_gen {

    typedef std::function<bool(const field::Mbf &)> include_fn_t;

    static include_fn_t default_include_fn() {
        return [](const field::Mbf &) { return true; };
    }

    static include_fn_t default_include_fn(const RankAllocator<WalkerTableRow> &ra) {
        return [&ra](const field::Mbf &mbf) { return mpi::i_am(ra.get_rank(mbf)); };
    }

    static include_fn_t default_include_fn(const Wavefunction &wf) {
        return default_include_fn(wf.m_ra);
    }

    struct Base {
        const BasisData m_bd;
        const size_t m_nelec;
        mutable suite::Mbfs m_mbf_work;
        const include_fn_t m_include_fn;


        Base(const BasisData& bd, size_t nelec, include_fn_t include_fn = default_include_fn());

        void add_if_included(Row &row, field::Mbf &mbf) {
            if (m_include_fn(m_mbf_work[mbf])) {
                row.push_back_jump();
                mbf = m_mbf_work[mbf];
                row.m_table->post_insert(row.index());
            }
        }

    protected:
        static void set_from_inds(field::FrmOnv& onv, const defs::inds& inds){
            onv = inds;
        }
        static void set_from_inds(field::FrmBosOnv& onv, const defs::inds& inds){
            onv.m_frm = inds;
        }
        static void set_from_inds(field::BosOnv& onv, const defs::inds& inds){
            onv.zero();
        }
        static void set_from_inds(field::FrmOnv& onv, const defs::inds& alpha_inds, const defs::inds& beta_inds){
            onv.set(alpha_inds, beta_inds);
        }
        static void set_from_inds(field::FrmBosOnv& onv, const defs::inds& alpha_inds, const defs::inds& beta_inds){
            onv.m_frm.set(alpha_inds, beta_inds);
        }
        static void set_from_inds(field::BosOnv& onv, const defs::inds& alpha_inds, const defs::inds& beta_inds){
            onv.zero();
        }
    };

    struct NoSym : Base {
        foreach::rtnd::Ordered<> m_foreach;

        NoSym(const BasisData& bd, size_t nelec, const include_fn_t &include_fn = default_include_fn());

    public:

        void operator()(Row &row, field::Mbf &mbf){
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

        SpinSym(const BasisData& bd, size_t nelec, int spin, const include_fn_t &include_fn = default_include_fn());

        void operator()(Row &row, field::Mbf &mbf);
    };

};


#endif //M7_CISPACES_H
