//
// Created by rja on 06/02/2021.
//

#ifndef M7_NDSEQUENCE_H
#define M7_NDSEQUENCE_H

#include "NdFormat.h"
#include "src/core/enumerator/ProductEnumerator.h"
#include "src/core/io/Logging.h"
#include "src/core/parallel/MPIAssert.h"

template<size_t nind>
struct NdSequence {
    const NdFormat<nind>& m_format;
    typedef std::array<size_t, nind> inds_t;
    const std::vector<inds_t> m_inds_vector;
    NdSequence(const NdFormat<nind>& format): m_format(format),
    m_inds_vector(make_inds_vector(format)){}

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

    static std::vector<std::array<size_t, nind>> make_inds_vector(const NdFormat<nind>& format) {
        const size_t limit = 2<<10;
        if (format.nelement() > limit){
            log::warn("Attempting to cache {} (> {}) Nd index arrays, this class is intended for small products of Nd extents",
                      format.nelement(), limit);
        }
        MPI_REQUIRE(format.nelement() <= limit, "shape product is too large, highly likely NdSequence is misapplied");
        std::vector<std::array<size_t, nind>> tmp;
        defs::inds shape(nind);
        std::copy_n(format.shape().begin(), nind, shape.begin());
        ProductEnumerator enumerator(std::move(shape));
        tmp.resize(format.nelement());
        defs::inds inds(nind);
        size_t iflat = ~0ul;
        while (enumerator.next(inds, iflat)) std::copy_n(inds.begin(), nind, tmp[iflat].begin());
        ASSERT(iflat == tmp.size());
        return tmp;
    }

};


#endif //M7_NDSEQUENCE_H
