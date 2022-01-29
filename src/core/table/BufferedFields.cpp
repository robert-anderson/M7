//
// Created by anderson on 1/27/22.
//

#include "BufferedFields.h"

BufferedFieldRow::BufferedFieldRow() :m_internal_row(), m_internal_buffer("", 1){}

BufferedFieldRow::BufferedFieldRow(const Row &row) : m_internal_row(row), m_internal_buffer("", 1){}

