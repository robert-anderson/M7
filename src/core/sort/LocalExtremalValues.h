//
// Created by rja on 09/05/2021.
//

#ifndef M7_LOCALEXTREMALVALUES_H
#define M7_LOCALEXTREMALVALUES_H


#include <src/core/util/utils.h>

template<typename row_t, typename T, size_t nind>
struct LocalExtremalValues {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    row_t m_work_row;
    fields::Numbers <T, nind> &m_work_row_field;
    row_t m_work_row_cmp;
    fields::Numbers <T, nind> &m_work_row_cmp_field;
    const bool m_largest, m_absval;
    const std::function<bool(const T &, const T &)> m_value_cmp_fn;
    ExtremalIndices m_xinds;

private:

    std::function<bool(const T &, const T &)> make_value_cmp_fn(dispatch_utils::BoolTag<false> is_unsigned) const {
        if (m_absval) {
            if (m_largest)
                return [](const T &v, const T &v_cmp) { return std::abs(v) < std::abs(v_cmp); };
            else
                return [](const T &v, const T &v_cmp) { return std::abs(v) > std::abs(v_cmp); };
        } else {
            if (m_largest)
                return [](const T &v, const T &v_cmp) { return v < v_cmp; };
            else
                return [](const T &v, const T &v_cmp) { return v > v_cmp; };
        }

    }

    std::function<bool(const T &, const T &)> make_value_cmp_fn(dispatch_utils::BoolTag<true> is_unsigned) const {
        if (m_largest)
            return [](const T &v, const T &v_cmp) { return v < v_cmp; };
        else
            return [](const T &v, const T &v_cmp) { return v > v_cmp; };
    }

    ExtremalIndices::cmp_fn_t make_cmp_fn() const {
        return [&](const size_t &irow, const size_t &irow_cmp) {
            static_cast<const Row &>(m_work_row).jump(irow);
            static_cast<const Row &>(m_work_row_cmp).jump(irow_cmp);
            for (size_t i = 0ul; i < m_work_row_field.nelement(); ++i) {
                if (m_work_row_field[i] != m_work_row_cmp_field[i])
                    return m_value_cmp_fn(m_work_row_field[i], m_work_row_cmp_field[i]);
            }
            return false;
        };
    }

public:

    LocalExtremalValues(row_t &row, fields::Numbers<T, nind> &field, bool largest, bool absval) :
            m_work_row(row),
            m_work_row_field(fields::identify(m_work_row, row, field)),
            m_work_row_cmp(row),
            m_work_row_cmp_field(fields::identify(m_work_row_cmp, row, field)),
            m_largest(largest), m_absval(absval),
            m_value_cmp_fn(make_value_cmp_fn(dispatch_utils::BoolTag<std::is_unsigned<T>::value>())),
            m_xinds(make_cmp_fn()) {
        ASSERT(static_cast<const FieldBase&>(m_work_row_field).belongs_to_row(m_work_row));
        ASSERT(static_cast<const FieldBase&>(m_work_row_cmp_field).belongs_to_row(m_work_row_cmp));
        m_xinds.reset(static_cast<const Row &>(row).m_table->m_hwm);
    }

    void find(size_t nfind) {
        m_xinds.find(nfind);
    }

};


#endif //M7_LOCALEXTREMALVALUES_H
