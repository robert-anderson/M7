//
// Created by rja on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include "FrmBosOnvConnection.h"

namespace conn {

    typedef FrmOnvConnection FrmOnv;
    typedef BosOnvConnection BosOnv;
    typedef FrmBosOnvConnection FrmBosOnv;

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_ind> Mbf;

    template<typename T=void> struct selector {typedef void type;};
    template<> struct selector<FrmOnvField> {typedef FrmOnv type;};
    template<> struct selector<FrmBosOnvField> {typedef FrmBosOnv type;};
    template<> struct selector<BosOnvField> {typedef BosOnv type;};

    template<typename T>
    using from_field_t = typename selector<T>::type;
}

#endif //M7_CONNECTIONS_H
