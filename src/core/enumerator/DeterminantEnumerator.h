//
// Created by rja on 04/11/2020.
//

#ifndef M7_DETERMINANTENUMERATOR_H
#define M7_DETERMINANTENUMERATOR_H

#include "SpinConDetEnumerator.h"
#include "SpinNonConDetEnumerator.h"

class DeterminantEnumerator : public Enumerator<views::Determinant> {
    std::unique_ptr<Enumerator<views::Determinant>> m_internal_enum;

public:
    DeterminantEnumerator(size_t nsite, size_t nelec) :
            m_internal_enum(new SpinNonConDetEnumerator(nsite, nelec)) {}

    DeterminantEnumerator(size_t nsite, size_t nelec, int spin) :
            m_internal_enum(new SpinConDetEnumerator(nsite, nelec, spin)) {}

    bool next_element(views::Determinant &result) override {
        return m_internal_enum->next(result);
    }
};


#endif //M7_DETERMINANTENUMERATOR_H
