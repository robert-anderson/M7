//
// Created by Robert J. Anderson on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include <M7_lib/table/BufferedFields.h>

#include "FrmBosOnvConnection.h"

namespace conn {

    typedef FrmOnvConnection FrmOnv;
    typedef BosOnvConnection BosOnv;
    typedef FrmBosOnvConnection FrmBosOnv;

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<uint_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_type_ind> Mbf;

    template<typename T=void>
    struct selector {
        typedef void type;
    };
    template<>
    struct selector<FrmOnvField> {
        typedef FrmOnv type;
    };
    template<>
    struct selector<FrmBosOnvField> {
        typedef FrmBosOnv type;
    };
    template<>
    struct selector<BosOnvField> {
        typedef BosOnv type;
    };
    template<>
    struct selector<buffered::FrmOnv> {
        typedef FrmOnv type;
    };
    template<>
    struct selector<buffered::FrmBosOnv> {
        typedef FrmBosOnv type;
    };
    template<>
    struct selector<buffered::BosOnv> {
        typedef BosOnv type;
    };

    template<typename T>
    using from_field_t = typename selector<T>::type;
}

#endif //M7_CONNECTIONS_H
