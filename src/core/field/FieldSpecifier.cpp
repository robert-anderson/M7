//
// Created by rja on 21/10/2020.
//

#include "FieldSpecifier.h"
#include "src/core/table/Table.h"

FieldSpecifier::FieldSpecifier(size_t element_size, const std::type_info &type_info) :
        m_data{element_size, type_info, {{"element size (bytes)", std::to_string(element_size)}}},
        m_null_buffer(element_size, 0){}