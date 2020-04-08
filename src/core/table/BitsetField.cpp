//
// Created by Robert John Anderson on 2020-03-26.
//

#include "BitsetField.h"

BitsetElement::BitsetElement(BitsetField *field, char *begin): Element(field, begin){}

void BitsetElement::set(const defs::pair &pair) {
    assert(pair.first < size());
    const auto offset = BitsetField::rectify_offset(pair);
    bit_utils::set(*(m_begin + offset.first), offset.second);
}

void BitsetElement::set(const size_t &ibit) {
    assert(ibit<nbit());
    set(defs::pair{0, ibit});
}

void BitsetElement::set(const defs::inds &inds) {
    for (auto const &ind: inds) set(ind);
}

void BitsetElement::clr(const defs::pair &pair) {
    assert(pair.first < size());
    const auto offset = BitsetField::rectify_offset(pair);
    bit_utils::clr(*(m_begin + offset.first), offset.second);
}

void BitsetElement::clr(const size_t &ibit) {
    clr(defs::pair{0, ibit});
}

bool BitsetElement::get(const defs::pair &pair) const {
    assert(pair.first < size());
    const auto offset = BitsetField::rectify_offset(pair);
    return bit_utils::get(*(m_begin + offset.first), offset.second);
}

bool BitsetElement::get(const size_t &ibit) const {
    return get({0, ibit});
}

std::string BitsetElement::to_string(size_t ibegin, size_t iend) {
    std::string result;
    for (size_t ibit = ibegin; ibit < iend; ++ibit) {
        result += get(ibit) ? "1" : "0";
    }
    return result;
}

size_t BitsetElement::nsetbit() const {
    size_t result = 0;
    for (size_t idataword = 0ul; idataword<m_field->element_dsize(); ++idataword){
        result+=bit_utils::nsetbit(get_dataword(idataword));
    }
    return result;
}

bool BitsetElement::is_zero() const {
    return nsetbit()==0;
}

std::string BitsetElement::to_string() {
    return to_string(0, nbit());
}

defs::pair BitsetField::rectify_offset(const defs::pair &pair) {
    /*
     * if the bit-offset in a (byte-offset, bit-offset) pair exceeds CHAR_BIT
     * this function will advance the byte-offset accordingly
     */
    return {pair.first + pair.second / CHAR_BIT, pair.second % CHAR_BIT};
}
