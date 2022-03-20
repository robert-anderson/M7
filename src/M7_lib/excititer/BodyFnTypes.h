//
// Created by rja on 25/08/2021.
//

#ifndef M7_BODYFNTYPES_H
#define M7_BODYFNTYPES_H

#include <M7_lib/connection/Connections.h>

namespace body_fn_types {
    /*
     * body function definitions
     * the names of these contain one-letter abbreviations of the arguments they expect:
     *  "c": the generated connection
     *  "d": the destination MBF
     *  "h": the associated hamiltonian matrix element
     */
    /**
     * this is the "canonical" prototype of the body function. all other prototypes are constructed as needed by the
     * adaptors
     */
    template<typename mbf_t>
    using fn_c_t = std::function<void(const conn::from_field_t<mbf_t> &)>;
    /**
     * body function which accepts connection and dst MBF as const refs
     */
    template<typename mbf_t>
    using fn_cd_t = std::function<void(const conn::from_field_t<mbf_t> &, const mbf_t &)>;
    /**
     * body function which accepts connection as const ref and matrix element by value
     */
    template<typename mbf_t>
    using fn_ch_t = std::function<void(const conn::from_field_t<mbf_t> &, defs::ham_t)>;
    /**
     * most complete prototype of body function which accepts connection and dst MBF as const refs and the matrix
     * element by value
     */
    template<typename mbf_t>
    using fn_cdh_t = std::function<void(const conn::from_field_t<mbf_t> &, const mbf_t &, defs::ham_t)>;
    /**
     * body function which accepts dst MBF as const ref only
     */
    template<typename mbf_t>
    using fn_d_t = std::function<void(const mbf_t &)>;
    /**
     * untemplated body function which accepts matrix element by value only
     */
    using fn_h_t = std::function<void(defs::ham_t)>;
    /**
     * body function which accepts dst MBF as const ref and matrix element by value
     */
    template<typename mbf_t>
    using fn_dh_t = std::function<void(const mbf_t &, defs::ham_t)>;
}

#endif //M7_BODYFNTYPES_H
