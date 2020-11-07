//
// Created by rja on 07/11/2020.
//

#include "StatsColumn.h"
#include "StatsFile.h"

StatsColumnBase::StatsColumnBase(StatsSpecifier* spec, size_t nelement, defs::inds shape, size_t nsubcolumn, std::string description) :
        m_nelement(nelement), m_shape(shape), m_nsubcolumn(nsubcolumn), m_description(description) {
    spec->add_column(this);
}
