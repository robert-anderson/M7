//
// Created by Robert J. Anderson on 1/25/22.
//

#include "CompositeField.h"
#include "M7_lib/util/String.h"

str_t CompositeFieldBase::prefix(str_t base, str_t prefix) {
    return string::prefix(base, prefix, '.');
}