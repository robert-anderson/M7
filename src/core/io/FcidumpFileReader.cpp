//
// Created by rja on 05/11/2020.
//

#include "FcidumpFileReader.h"

FcidumpFileReader::FcidumpFileReader(const std::string &fname, bool spin_major) :
        HamiltonianFileReader(fname, 4, false),
        m_spin_major(spin_major),
        m_nelec(read_header_int(fname, "NELEC")),
        m_orbsym(read_header_array(fname, "ORBSYM", -1))
{
    REQUIRE_EQ(m_orbsym.size(), m_norb, "invalid ORBSYM specified in FCIDUMP file");
    set_symm_and_rank(fname);

    if (m_spin_resolved) {
        defs::inds inds(4);
        defs::ham_t v;
        while (next(inds, v)) {
            if (!consts::float_is_zero(std::abs(v))) {
                if (((inds[0] < m_nspatorb) != (inds[1] < m_nspatorb)) ||
                    ((inds[2] < m_nspatorb) != (inds[3] < m_nspatorb))) {
                    // spin non-conserving example found
                    if (nset_ind(inds)==2) m_spin_conserving_1e = false;
                    else m_spin_conserving_2e = false;
                }
            }
        }
        SparseArrayFileReader<defs::ham_t>::reset(); // go back to beginning of entries
    }
    if (m_spin_conserving_1e) log::info("FCIDUMP file conserves spin in 1 particle integrals");
    else log::info("FCIDUMP file does NOT conserve spin in 1 particle integrals");
    if (m_spin_conserving_2e) log::info("FCIDUMP file conserves spin in 2 particle integrals");
    else log::info("FCIDUMP file does NOT conserve spin in 2 particle integrals");
    log::info("FCIDUMP file contains 2 particle integrals of maximum excitation rank " + std::to_string(m_int_2e_rank));
}

bool FcidumpFileReader::spin_conserving() const {
    return m_spin_conserving_1e && m_spin_conserving_2e;
}

void FcidumpFileReader::inds_to_orbs(defs::inds &inds) {
    m_inds_to_orbs(inds);
}


void FcidumpFileReader::set_symm_and_rank(const std::string &filename) {
    m_int_2e_rank = 0;
    auto get_rank = [](const defs::inds& inds){
        // [ij|kl]
        const size_t i=std::min(inds[0], inds[1]);
        const size_t j=std::max(inds[0], inds[1]);
        const size_t k=std::min(inds[2], inds[3]);
        const size_t l=std::max(inds[2], inds[3]);
        if (i==k && j==l) return 0ul;
        else if (i!=k && j!=l) return 2ul;
        else return 1ul;
    };

    log::info("Determining permutational symmetry of integral file entries");
    SparseArrayFileReader<defs::ham_t> reader(filename, 4, false);
    defs::inds inds(4);
    defs::ham_t value;
    // this will eventually hold all orderings of the first example of an
    // integral with 4 distinct indices
    std::array<defs::inds, 8> inds_distinct;
    m_isymm = 0ul;
    while (reader.next(inds, value)) {
        if (std::all_of(inds.begin(), inds.end(), [](size_t i){return i>0;})) {
            // we have a two body integral
            if (m_int_2e_rank<2) {
                size_t rank = get_rank(inds);
                if (rank > m_int_2e_rank) m_int_2e_rank = rank;
            }

            if (!m_isymm) {
                inds_distinct[0].assign(inds.begin(), inds.end());
                std::sort(inds.begin(), inds.end());
                // still looking for an example of four distinct indices
                if (std::adjacent_find(inds.begin(), inds.end()) == inds.end()) {
                    for (size_t i = 1ul; i < orderings.size(); ++i) {
                        for (size_t j = 0ul; j < 4; ++j)
                            inds_distinct[i].push_back(inds_distinct[0][orderings[i][j]]);
                    }
                    m_isymm++;
                }
            } else {
                for (auto tmp : inds_distinct) {
                    if (std::equal(inds.begin(), inds.end(), tmp.begin())){
                        m_isymm++;
                        break;
                    }
                }
            }
        }
    }
    if (!m_isymm){
        log::info("Permutational symmetry of integral file could not be determined");
        m_isymm = 1;
    } else {
        log::info("Permutational symmetry of integral file found to be {}", 8/m_isymm);
        m_isymm/=8;    // 8->1; 4->2; 2->4; 1->8
    }
}

size_t FcidumpFileReader::ranksig(const defs::inds &inds) const {
    auto nset_inds = HamiltonianFileReader::nset_ind(inds);
    return exsig_utils::encode(nset_inds/2, nset_inds/2, 0, 0);
}

size_t FcidumpFileReader::exsig(const defs::inds &inds, const size_t& ranksig) const {
    switch (ranksig) {
        case 0ul: return 0ul;
        case exsig_utils::ex_single:
            return inds[0]==inds[1] ? 0ul : exsig_utils::ex_single;
            case exsig_utils::ex_double:
            return inds[0]==inds[2] ?
                (inds[1]==inds[3] ? 0ul : exsig_utils::ex_single):
                (inds[1]==inds[3] ? exsig_utils::ex_single : exsig_utils::ex_double);
        default:
            return ~0ul;
    }
}
