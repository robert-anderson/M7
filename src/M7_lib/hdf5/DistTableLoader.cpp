//
// Created by rja on 02/05/23.
//

#include "DistTableLoader.h"

uint_t DistTableLoader::max_item_size() const {
    uint_t max = 0ul;
    for (auto& field_loader: m_fields_loaders) {
        const auto size = field_loader.m_loader->m_format.m_local.m_item.m_size;
        if (size > max) max = size;
    }
    return max;
}

uint_t DistTableLoader::nitem_next(uint_t max_nitem_per_op) const {
    return m_fields_loaders.front().m_loader->nitem_next(max_nitem_per_op);
}

v_t<DistTableLoader::NameFieldPair> DistTableLoader::make_pairs(const v_t<FieldBase*>& fields) {
    v_t<NameFieldPair> tmp;
    tmp.reserve(fields.size());
    for (auto field: fields) tmp.emplace_back(field->m_name, field);
    return tmp;
}

uint_t DistTableLoader::nitem() const {
    return m_fields_loaders.front().m_loader->m_format.m_nitem;
}

uint_t DistTableLoader::nitem_local() const {
    return m_fields_loaders.front().m_loader->m_format.m_local.m_nitem;
}

DistTableLoader::DistTableLoader(const hdf5::NodeReader& nr, v_t<DistTableLoader::NameFieldPair> pairs) {
    m_loaders.reserve(pairs.size());
    m_fields_loaders.reserve(pairs.size());
    for (auto& pair: pairs) {
        m_loaders.emplace_back(nr, pair.m_name, true, true);
        REQUIRE_EQ_ALL(m_loaders.back().m_format.m_nitem, m_loaders.front().m_format.m_nitem,
                       "number of elements should be constant for all fields");
        m_fields_loaders.push_back({pair.m_field, &m_loaders.back()});
    }
}

DistTableLoader::DistTableLoader(const hdf5::NodeReader& nr, const v_t<FieldBase*>& fields) :
        DistTableLoader(nr, make_pairs(fields)){}
