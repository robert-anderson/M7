//
// Created by Robert J. Anderson on 06/02/2021.
//

#ifndef M7_NDSEQUENCE_H
#define M7_NDSEQUENCE_H

#include "NdFormat.h"
#include <M7_lib/enumerator/ProductEnumerator.h"
#include <M7_lib/io/Logging.h"
#include <M7_lib/parallel/MPIAssert.h"

template<size_t nind>
struct NdSequence {
    NdFormat<nind> m_format;
    typedef uinta_t<nind> inds_t;
    const std::vector<inds_t> m_inds_vector;
    NdSequence(uinta_t<nind> shape): m_format(shape),
    m_inds_vector(make_inds_vector(m_format)){}

    struct Cursor {
        const NdSequence<nind>& m_sequence;
        typename std::vector<inds_t>::const_iterator m_inds;

        Cursor(const NdSequence<nind>& sequence): m_sequence(sequence){
            to_front();
        }

        size_t ielement() const {
            return std::distance(m_sequence.m_inds_vector.cbegin(), m_inds);
        }

        const inds_t& inds() const {
            return *m_inds;
        }

        operator bool() const {
            return m_inds!=m_sequence.m_inds_vector.end();
        }

        size_t nelement() const {
            return m_sequence.m_inds_vector.size();
        }

        template<typename ...Args>
        Cursor& to(Args... inds) {
            m_inds = m_sequence.m_inds_vector.cbegin() + m_sequence.m_format.flatten(inds...);
            return *this;
        }

        Cursor& operator++(){
            m_inds++;
            return *this;
        }

        Cursor& operator++(int){
            ++m_inds;
            return *this;
        }

        void to_front() {
            m_inds = m_sequence.m_inds_vector.cbegin();
        }
    };

    Cursor cursor() const{
        return Cursor(*this);
    }

    static std::vector<uinta_t<nind>> make_inds_vector(const NdFormat<nind>& format) {
        const size_t limit = 2<<10;
        if (format.m_nelement > limit){
            log::warn("Attempting to cache {} (> {}) Nd index arrays, this class is intended for small products of Nd extents",
                      format.m_nelement, limit);
        }
        REQUIRE_LE(format.m_nelement, limit, "shape product is too large, highly likely NdSequence is misapplied");
        std::vector<uinta_t<nind>> tmp;
        defs::uintv_t shape(nind);
        std::copy_n(format.m_shape.cbegin(), nind, shape.begin());
        ProductEnumerator enumerator(std::move(shape));
        tmp.resize(format.m_nelement);
        defs::uintv_t inds(nind);
        size_t iflat = ~0ul;
        while (enumerator.next(inds, iflat)) std::copy_n(inds.begin(), nind, tmp[iflat].begin());
        ASSERT(iflat == tmp.size());
        return tmp;
    }

};


#endif //M7_NDSEQUENCE_H
