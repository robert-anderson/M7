//
// Created by rja on 27/09/2021.
//

#include "BosOnvField.h"

BosOnvField::BosOnvField(Row *row, size_t nmode, std::string name) :
        NdNumberField<defs::bos_occ_t, 1>(row, {{nmode}, {"boson mode occupations"}}, name) {}

BosOnvField::BosOnvField(Row *row, BasisDims bd, std::string name) : BosOnvField(row, bd.m_nmode, name) {
    bd.require_pure_bos();
}