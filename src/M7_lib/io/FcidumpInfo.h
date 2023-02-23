//
// Created by rja on 25/06/22.
//

#ifndef M7_FCIDUMPINFO_H
#define M7_FCIDUMPINFO_H

#include <utility>
#include <M7_lib/integrals/IntegralArray2e.h>

#include "M7_lib/hdf5/File.h"
#include "FortranNamelistReader.h"

/**
 * this only involves reading the header of the FCIDUMP file, so it can be done on all ranks - no need for comms
 */
struct FcidumpInfo {
    /**
     * Supported I/O formats for Fcidumps
     */
    enum Implementation {
        CSV, MolcasHDF5
    };
    /**
     * Supported protocols for specifying H coefficients in a spin-unrestricted basis
     */
    enum UnrestrictStyle {
        SpinMinor,    // consecutive pairs of spinorbs correspond to the up and down spin functions of the same orbital
                      // (this is the style used internally by NECI), default.
        SpinMajor,    // the first [1, nsite] spin orbitals are the up spinorbs, the remainder [nsite+1, 2*nsite]
                      // (1-based counting) are the down spinorbs
        SpinBlocks,   // all indices are spatial, but the integrals are split into blocks. for the 2-body ints:
                      // (uu|uu), (uu|dd), (dd|dd), and for the 1-body ints: h_uu, h_dd
                      // blocks are delimited by a line of zeros
    };

    /**
     * convenient conversion for (non-quantum chemical) coefficient reading contexts where blocks are not applicable
     */
    static UnrestrictStyle ur_style(bool spin_major) {
        return spin_major ? SpinMajor : SpinMinor;
    }
    static UnrestrictStyle ur_style(const str_t& str){
        if (str=="minor") return SpinMinor;
        if (str=="major") return SpinMajor;
        if (str=="blocks") return SpinBlocks;
        ABORT_ALL("invalid unrestrict style");
        return {};
    }
    static str_t ur_desc(UnrestrictStyle style){
        switch (style) {
            case SpinMinor: return "spin minor";
            case SpinMajor: return "spin major";
            case SpinBlocks: return "spin blocks (Molpro)";
        }
        return {};
    }

    const str_t m_fname;
    const Implementation m_impl;
    const UnrestrictStyle m_ur_style;
    const integrals_2e::syms::Sym m_init_2e_perm_sym;
    /*
     * metadata fields
     */
    bool m_uhf, m_relativistic, m_spin_resolved;
    uint_t m_nelec, m_nsite, m_nspinorb, m_norb_distinct;
    int m_ms2;
    uintv_t m_orbsym;

private:
    void set_metadata(bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, const uintv_t& orbsym);

public:

    FcidumpInfo(str_t fname, UnrestrictStyle ur_style=SpinMinor,
                integrals_2e::syms::Sym init_2e_perm_sym=integrals_2e::syms::DHR);

};

struct EbdumpInfo : FcidumpInfo {
    uint_t m_nmode;
    EbdumpInfo(str_t fname, UnrestrictStyle ur_style=SpinMinor);

};

#endif //M7_FCIDUMPINFO_H
