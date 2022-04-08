//
// Created by rja on 26/07/2021.
//

#ifndef M7_SUITES_H
#define M7_SUITES_H

#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/connection/Connections.h>

/**
 * classes for "workspace" objects that are MBF type dependent. Objects for every implemented MBF type are allocated as
 * fields in a row, and then the appropriate object for the MBF type in use can be either accessed directly from these
 * fields, or by overloading
 */
//TODO (RJA): rename to "working"
namespace suite {
    struct MbfsRow : Row {
        field::FrmOnv m_frm;
        field::FrmBosOnv m_frmbos;
        field::BosOnv m_bos;
        MbfsRow(const BasisData& bd):
            m_frm(this, bd.m_frm, "fermion ONV"),
            m_frmbos(this, bd, "fermion-boson ONV"),
            m_bos(this, bd.m_bos, "boson ONV"){}
    };

    struct Mbfs : BufferedTable<MbfsRow>{

        Mbfs(const BasisData& bd): BufferedTable<MbfsRow>("Work space for MBFs", {{bd}}){
            m_row.push_back_jump();
        }

        field::FrmOnv& operator[](const field::FrmOnv& mbf){
            return m_row.m_frm;
        }
        field::FrmBosOnv& operator[](const field::FrmBosOnv& mbf){
            return m_row.m_frmbos;
        }
        field::BosOnv& operator[](const field::BosOnv& mbf){
            return m_row.m_bos;
        }
    };

    struct Conns {
        conn::FrmOnv m_frmonv;
        conn::BosOnv m_bosonv;
        conn::FrmBosOnv m_frmbosonv;
        Conns(const BasisData& bd): m_frmonv(bd.m_frm), m_bosonv(bd.m_bos), m_frmbosonv(bd){}

        conn::FrmOnv& operator[](const field::FrmOnv& mbf){
            return m_frmonv;
        }
        conn::BosOnv& operator[](const field::BosOnv& mbf){
            return m_bosonv;
        }
        conn::FrmBosOnv& operator[](const field::FrmBosOnv& mbf){
            return m_frmbosonv;
        }
    };

    struct ComOps {
        com_ops::Frm m_frm;
        com_ops::FrmBos m_frmbos;
        com_ops::Bos m_bos;
        ComOps(const BasisData& bd): m_frm(bd.m_frm.m_nsite), m_frmbos(bd), m_bos(bd.m_bos.m_nmode){}

        com_ops::Frm& operator[](const field::FrmOnv& mbf){
            return m_frm;
        }
        com_ops::FrmBos& operator[](const field::FrmBosOnv& mbf){
            return m_frmbos;
        }
        com_ops::Bos& operator[](const field::BosOnv& mbf){
            return m_bos;
        }
    };

}

#endif //M7_SUITES_H
