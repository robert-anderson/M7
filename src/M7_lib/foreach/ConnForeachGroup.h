//
// Created by Robert J. Anderson on 12/04/2022.
//

#ifndef M7_CONNFOREACHGROUP_H
#define M7_CONNFOREACHGROUP_H

#include "ConnForeach.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"

class ConnForeachGroup {
    /**
     * forward list containing all iterators
     */
    conn_foreach::base_list_t m_list;

public:
    explicit ConnForeachGroup(const Hamiltonian &ham);

    using function_t = conn_foreach::Base::function_t;

    template<typename mbf_t>
    void loop(conn::from_field_t<mbf_t> &conn, const mbf_t &src, const function_t &fn) const {
        for(const auto& foreach : m_list) {
            conn.clear();
            foreach->loop(conn, src, fn);
        }
    }

    void log() const;

};

#endif //M7_CONNFOREACHGROUP_H