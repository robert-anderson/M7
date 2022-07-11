//
// Created by rja on 10/07/22.
//

#include "ConfComponents.h"

conf_components::Path::Path(std::list<str_t> list) :
        m_list(std::move(list)),
        m_string(string::join(v_t<str_t>(m_list.cbegin(), m_list.cend()), ".")){}

conf_components::Path conf_components::Path::operator+(const str_t& name) const {
    auto tmp = m_list;
    tmp.push_back(name);
    return {tmp};
}
