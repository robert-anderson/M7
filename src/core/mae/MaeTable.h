//
// Created by rja on 11/08/2021.
//

#ifndef M7_MAETABLE_H
#define M7_MAETABLE_H

#include "src/core/field/Fields.h"

using namespace fields;

namespace mae_inds {
    struct Frm : MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> {
        typedef MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> base_t;
        const std::string m_name;
        Numbers<defs::mev_ind_t, 1> &m_ann;
        Numbers<defs::mev_ind_t, 1> &m_cre;

        Frm(Row *row, size_t nann, size_t ncre, std::string name = "");

        Frm(const Frm &other);

        Frm &operator=(const Frm &other);

        Frm &operator=(const std::pair<defs::inds, defs::inds> &inds);

        /**
         * The RDM entry keys represent the ascending-ordered sector of the normal-ordered pure expectation values,
         * the full spin-resolved RDM can be constructed from this via permutations of the indices.
         * @return
         *  true if the SQ creation and annihilation operator index vectors are properly ordered
         */
        bool is_ordered() const;

        /**
         * all elements of the RDM have the same rank, but not the same excitation level. If the number of indices in
         * common between the creation and annihilation operators is zero, the excitation level is the same as the rank
         * @return
         *  number of indices in common between ascending ordered creation and annihilation spin orbital operator strings
         */
        size_t ncommon_sq_op_ind() const;

        void common_sq_op_inds(defs::inds &common) const;
    };
}

#endif //M7_MAETABLE_H
