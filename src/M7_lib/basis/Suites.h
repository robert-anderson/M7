//
// Created by Robert J. Anderson on 26/07/2021.
//

#ifndef M7_SUITES_H
#define M7_SUITES_H

#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/connection/Connections.h>
#include "BasisData.h"

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
        MbfsRow(const sys::Sector& sector);
    };

    struct Mbfs : BufferedTable<MbfsRow>{

        Mbfs(const sys::Sector& sector);

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
        explicit Conns(const sys::Size& size);

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
        ComOps(const sys::Size& size);

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
