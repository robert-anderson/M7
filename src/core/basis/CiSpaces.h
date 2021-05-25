//
// Created by rja on 08/05/2021.
//

#ifndef M7_CISPACES_H
#define M7_CISPACES_H


#include <src/core/enumerator/Enumerator.h>
#include <src/core/table/BufferedFields.h>
#include <src/core/parallel/RankAllocator.h>
#include <src/core/dynamics/WalkerTable.h>
#include <src/core/dynamics/Wavefunction.h>

namespace ci_gen {

    typedef std::function<bool(const fields::Onv<>&)> include_fn_t;
    static include_fn_t default_include_fn() {
        return [](const fields::Onv<>&){return true;};
    }

    static include_fn_t default_include_fn(const RankAllocator<WalkerTableRow>& ra) {
        return [&ra](const fields::Onv<>& onv){return mpi::i_am(ra.get_rank(onv));};
    }

    static include_fn_t default_include_fn(const Wavefunction& wf){
        return default_include_fn(wf.m_ra);
    }

    struct Base {
        const size_t m_nsite, m_nelec;
        mutable buffered::Onv<> m_onv_work;
        const include_fn_t m_include_fn;


        Base(size_t nsite, size_t nelec, include_fn_t include_fn=default_include_fn());

        void add_if_included(Row& row, fields::Onv<>& onv){
            if (m_include_fn(m_onv_work)) {
                row.push_back_jump();
                onv = m_onv_work;
                row.m_table->post_insert(row.m_i);
            }
        }
    };

    struct NoSym : Base {
        foreach::rtnd::Ordered<> m_foreach;
        NoSym(size_t nsite, size_t nelec, const include_fn_t& include_fn=default_include_fn());

        void operator()(Row& row, fields::Onv<>& onv);
    };

    struct SpinSym : Base {
        foreach::rtnd::Ordered<> m_foreach_alpha;
        foreach::rtnd::Ordered<> m_foreach_beta;
        SpinSym(size_t nsite, size_t nelec, int spin, const include_fn_t& include_fn=default_include_fn());

        void operator()(Row& row, fields::Onv<>& onv);
    };
};


#endif //M7_CISPACES_H
