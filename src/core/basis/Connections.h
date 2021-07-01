//
// Created by rja on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"
#include "FermiBosConnection.h"

namespace conn {

    template <bool enable_bosons=defs::enable_bosons>
    using Basic = typename std::conditional<enable_bosons, AntisymFermiBosConnection, FermionOnvConnection>::type;

    template <bool enable_bosons=defs::enable_bosons>
    using Antisym = typename std::conditional<enable_bosons, AntisymFermiBosConnection, AntisymFermionOnvConnection>::type;

    using Boson = BosonOnvConnection;

    struct ExlvlSector {
        defs::exlvl_t m_nann;
        defs::exlvl_t m_ncre;
        ExlvlSector(const size_t& nann, const size_t& ncre){
            m_nann = utils::safe_narrow<defs::exlvl_t>(nann);
            m_ncre = utils::safe_narrow<defs::exlvl_t>(ncre);
        }
        ExlvlSector(const defs::exlvl_t& nann, const defs::exlvl_t& ncre): m_nann(nann), m_ncre(ncre){}
    };

    struct FermionExlvl : ExlvlSector {
        FermionExlvl(const size_t& nann, const size_t& ncre): ExlvlSector(nann, ncre){}
        FermionExlvl(const defs::exlvl_t& nann, const defs::exlvl_t& ncre): ExlvlSector(nann, ncre){}
        FermionExlvl(const fields::Number<defs::exlvl_t>& field):
            ExlvlSector(field[0], field[1]){}
    };

    struct FermiBosExlvl {
        ExlvlSector m_frm, m_bos;
        FermiBosExlvl(ExlvlSector&& frm, ExlvlSector&& bos): m_frm(frm), m_bos(bos){}
        FermiBosExlvl(const fields::Number<defs::exlvl_t>& field):
                FermiBosExlvl({field[0], field[1]}, {field[2], field[3]}){}
    };

    template <bool enable_bosons=defs::enable_bosons>
    using Exlvl = typename std::conditional<enable_bosons, FermiBosExlvl, FermionExlvl>::type;

    template <bool enable_bosons=defs::enable_bosons>
    static constexpr size_t nelement_exlvl() {return enable_bosons ? 4ul:2ul;}

    static Exlvl<0> exlvl(const Basic<0>& conn){
        return {conn.nann(), conn.nann()};
    }
//    Exlvl<1> exlvl(const Basic<1>& conn){
//        return {conn.nann(), conn.nann()}, conn.m_bonvconn.nexcit()};
//    }

//    using FermionOnv = FermionOnvConnection;
//    using AsFermionOnv = AntisymFermionOnvConnection;
//
//    //using FermiBosOnv = FermiBosConnection;
//    using AsFermiBosOnv = AntisymFermiBosConnection;
//    //using Onv = std::conditional<defs::enable_bosons, FermiBosOnv, FermionOnv>::type;
//    using AsOnv = std::conditional<defs::enable_bosons, AsFermiBosOnv, AsFermionOnv>::type;
}

#endif //M7_CONNECTIONS_H
