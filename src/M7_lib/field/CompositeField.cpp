//
// Created by anderson on 1/25/22.
//

#include "CompositeField.h"

std::string CompositeFieldBase::prefix(std::string base, std::string prefix) {
    return string_utils::prefix(base, prefix, '.');
}
