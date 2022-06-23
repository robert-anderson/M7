//
// Created by Robert J. Anderson on 1/25/22.
//

#include "CompositeField.h"
#include "M7_lib/util/String.h"

std::string CompositeFieldBase::prefix(std::string base, std::string prefix) {
    return string::prefix(base, prefix, '.');
}