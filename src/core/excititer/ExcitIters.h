//
// Created by rja on 25/08/2021.
//

#ifndef M7_EXCITITERS_H
#define M7_EXCITITERS_H

#include <src/core/util/Foreach.h>
#include "ExcitIter.h"

/**
 * classes of different iterator implementations, defined in a namespace to allow for names to be shared with
 * equivalently-connected excitation generation implementations without conflict
 */
namespace excititers {

    struct Frm : public ExcitIter {
        using ExcitIter::foreach;
    protected:
        foreach::rtnd::Ordered<> m_cre_loop;
        foreach::rtnd::Ordered<> m_ann_loop;

        fn_c_t<field::FrmOnv> convert(conn::FrmBosOnv &work_conn, const fn_c_t<field::FrmBosOnv> &fn);

    public:
        Frm(const Hamiltonian &ham, size_t exsig);

        void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<field::FrmOnv> &body) override;

        void foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) override;

        void foreach(const BosOnv &src, conn::BosOnv &conn, const fn_c_t<BosOnv> &body) override {}
    };


    struct Ladder : public ExcitIter {
        const bool m_cre;

        Ladder(const Hamiltonian &ham, size_t exsig);

        void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<field::FrmOnv> &body) override {}

        void foreach(const field::FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<field::FrmBosOnv> &body) override {}

        void foreach(const BosOnv &src, conn::BosOnv &conn, const fn_c_t<BosOnv> &body) override {}
    };

    struct Bos : public ExcitIter {
    protected:
        fn_c_t<field::BosOnv> convert(conn::FrmBosOnv &work_conn, const fn_c_t<field::FrmBosOnv> &fn);
    public:

        Bos(const Hamiltonian &ham, size_t exsig) : ExcitIter(ham, exsig) {
            REQUIRE_TRUE(exsig_utils::is_pure_bos(exsig), "excitation signature should not have fermionic operators")
            REQUIRE_TRUE(exsig_utils::conserves_nbos(exsig), "excitation signature should conserve boson number");
            REQUIRE_EQ(exsig_utils::decode_nbos_ann(exsig), 2ul, "only currently implemented for doubles");
        }

        void foreach(const FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<FrmOnv> &body) override {}

        void foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) override;

        void foreach(const BosOnv &src, conn::BosOnv &conn, const fn_c_t<BosOnv> &body) override;
    };
}


#endif //M7_EXCITITERS_H
