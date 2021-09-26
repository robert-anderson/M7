//
// Created by rja on 26/07/2021.
//

#ifndef M7_SUITES_H
#define M7_SUITES_H

#include "src/core/table/BufferedFields.h"
#include "src/core/connection/Connections.h"

namespace suite {

    struct Mbfs {
        buffered::FrmOnv m_frmonv;
        buffered::BosOnv m_bosonv;
        buffered::FrmBosOnv m_frmbosonv;
        //buffered::FrmCsf m_frmcsf;
        Mbfs(BasisDims bd): m_frmonv(bd), m_bosonv(bd), m_frmbosonv(bd){}//, m_frmcsf(nsite){}

        field::FrmOnv& operator[](const field::FrmOnv& mbf){
            return m_frmonv;
        }
        field::BosOnv& operator[](const field::BosOnv& mbf){
            return m_bosonv;
        }
        field::FrmBosOnv& operator[](const field::FrmBosOnv& mbf){
            return m_frmbosonv;
        }
    };

    struct Conns {
        conn::FrmOnv m_frmonv;
        conn::BosOnv m_bosonv;
        conn::FrmBosOnv m_frmbosonv;
        Conns(BasisDims bd): m_frmonv(bd), m_bosonv(bd), m_frmbosonv(bd){}

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

}

#endif //M7_SUITES_H
