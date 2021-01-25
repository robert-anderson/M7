//
// Created by rja on 04/11/2020.
//

#include "BosonOnvSpecifier.h"

const size_t &BosonOnvSpecifier::nmode() const {
    return m_format.extent(0);
}

BosonOnvSpecifier::View BosonOnvSpecifier::operator()(char *ptr) const {
    return View(*this, ptr);
}

BosonOnvSpecifier::View::View(const BosonOnvSpecifier &field, char *ptr) : NumericArraySpecifier<uint8_t, 1>::View(field, ptr) {}

size_t BosonOnvSpecifier::View::nboson() const {
    return std::accumulate(
            reinterpret_cast<uint8_t *>(m_ptr),
            reinterpret_cast<uint8_t *>(m_ptr) + nelement(), 0ul);
}

size_t BosonOnvSpecifier::View::nmode() const {
    return spec().nmode();
}