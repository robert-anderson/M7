//
// Created by rja on 02/11/2020.
//

#include "BitsetSpecifier.h"

/*
 * View::BitView
 */
BitsetSpecifier::View::BitView::BitView(const BitsetSpecifier::View &view, const size_t &ibit):
        m_view(std::unique_ptr<View>(new View(view))), m_ibit(ibit){}

BitsetSpecifier::View::BitView::BitView(const BitsetSpecifier::View::BitView &bv):
        m_view(std::unique_ptr<View>(new View(*bv.m_view))), m_ibit(bv.m_ibit) {
}

BitsetSpecifier::View::BitView &BitsetSpecifier::View::BitView::operator=(bool t) {
    if (t) m_view->set(m_ibit);
    else m_view->clr(m_ibit);
    return *this;
}

BitsetSpecifier::View::BitView::operator bool() const {
    return m_view->get(m_ibit);
}


/*
 * View
 */
BitsetSpecifier::View::View(const BitsetSpecifier &field, char *ptr) : ColumnSpecifier::View(field, ptr){}

BitsetSpecifier::View::View(const BitsetSpecifier::View &other) : ColumnSpecifier::View(other){
    ASSERT(other.nbit()==nbit());
}

BitsetSpecifier::View::bitview_t BitsetSpecifier::View::operator[](const size_t &ibit) {
    return BitView(*this, ibit);
}

BitsetSpecifier::View::const_bitview_t BitsetSpecifier::View::operator[](const size_t &ibit) const {
    return BitView(*this, ibit);
}

const size_t &BitsetSpecifier::View::nbit() const {
    return static_cast<const BitsetSpecifier &>(m_spec).m_nbit;
}

const size_t &BitsetSpecifier::View::ndataword() const {
    return static_cast<const BitsetSpecifier &>(m_spec).m_ndataword;
}

bool BitsetSpecifier::View::get(const size_t &ibit) const {
    ASSERT(ibit < nbit());
    return bit_utils::get(*cdptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
}

defs::data_t BitsetSpecifier::View::get_dataword(const size_t &idataword) const {
    auto tmp = *cdptr(idataword);
    if (idataword + 1 == ndataword()) {
        auto n = nbit() - defs::nbit_data * idataword;
        tmp = bit_utils::truncate(tmp, n);
    }
    return tmp;
}

defs::data_t BitsetSpecifier::View::get_antidataword(const size_t &idataword) const {
    auto tmp = ~*cdptr(idataword);
    if (idataword + 1 == ndataword()) {
        auto n = nbit() - defs::nbit_data * idataword;
        tmp = bit_utils::truncate(tmp, n);
    }
    return tmp;
}

size_t BitsetSpecifier::View::nsetbit() const {
    size_t result = 0;
    for (size_t idataword = 0ul; idataword<ndataword(); ++idataword){
        result+=bit_utils::nsetbit(get_dataword(idataword));
    }
    return result;
}

std::string BitsetSpecifier::View::to_string() const {
    std::string res;
    res.reserve(nbit());
    for (size_t i = 0ul; i < nbit(); ++i) res += get(i) ? "1" : "0";
    return res;
}

void BitsetSpecifier::View::set(const size_t &ibit) {
    ASSERT(ibit < nbit());
    bit_utils::set(*dptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
}

void BitsetSpecifier::View::set(const defs::inds &setinds) {
    for (auto &i: setinds) set(i);
}

void BitsetSpecifier::View::clr(const size_t &ibit) {
    ASSERT(ibit < nbit());
    bit_utils::clr(*dptr(ibit / defs::nbit_data), ibit % defs::nbit_data);
}

BitsetSpecifier::View &BitsetSpecifier::View::operator=(const defs::inds &ibits) {
    zero();
    for (auto& ibit : ibits) set(ibit);
    return *this;
}

BitsetSpecifier::View &BitsetSpecifier::View::operator=(const BitsetSpecifier::View &other) {
    ASSERT(other.nbit()==nbit())
    ColumnSpecifier::View::operator=(other);
    return *this;
}


/*
 * BitsetSpecifier
 */
BitsetSpecifier::BitsetSpecifier(size_t nbit) :
        ColumnSpecifier(defs::nbyte_data * integer_utils::divceil(nbit, defs::nbit_data),
                        typeid(BitsetSpecifier)), m_nbit(nbit),
        m_ndataword(integer_utils::divceil(nbit, defs::nbit_data)){
    m_data.m_details["type"] = "Bitset";
    m_data.m_details["number of bits"] = std::to_string(m_nbit);
}

std::string BitsetSpecifier::element_string(char *ptr) const {
    return View(*this, ptr).to_string();
}

BitsetSpecifier::cview_t BitsetSpecifier::operator()(char *ptr) const {
    return View(*this, ptr);
}

BitsetSpecifier::view_t BitsetSpecifier::operator()(char *ptr) {
    return View(*this, ptr);
}
