////
//// Created by RJA on 26/10/2020.
////
//
//#ifndef M7_FLAG_H
//#define M7_FLAG_H
//
//#include "src/core/nd/NdFormat.h"
//
//struct FlagFieldX;
//
//struct Flag {
//    const size_t m_offset;
//    const size_t m_nbit;
//    const std::string m_description;
//    Flag(FlagFieldX* field, size_t offset, size_t nbit, std::string description):
//    m_offset(offset), m_nbit(nbit), m_description(std::move(description)){}
//};
//
//
//template <size_t nind>
//struct NdFlag : Flag {
//    NdFormat<nind> m_format;
//    Flag(FlagFieldX* field, size_t offset, std::array<size_t, nind> shape, std::string description):
//};
//
//
//#endif //M7_FLAG_H
