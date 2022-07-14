//
// Created by rja on 14/07/22.
//

#ifndef M7_FCI_H
#define M7_FCI_H

#include "M7_lib/hamiltonian/Hamiltonian.h"
#include "M7_lib/foreach/MbfForeach.h"


namespace fci {

    struct BasisIters {
        typedef mbf_foreach::Base single_poly_t;
        typedef std::unique_ptr<single_poly_t> single_poly_ptr_t;
        typedef mbf_foreach::PairBase pair_poly_t;
        typedef std::unique_ptr<pair_poly_t> pair_poly_ptr_t;

        single_poly_ptr_t m_single = nullptr;
        pair_poly_ptr_t m_pair = nullptr;

    private:

        template<typename single_t>
        void set_pair(const single_t& single, const mbf_foreach::frm::Base&) {
            typedef mbf_foreach::frm::Pair<single_t> pair_t;
            m_pair = smart_ptr::make_poly_unique<pair_poly_t, pair_t>(single);
        }
        template<typename single_t>
        void set_pair(const single_t& single, const mbf_foreach::bos::Base&) {
            typedef mbf_foreach::bos::Pair<single_t> pair_t;
            m_pair = smart_ptr::make_poly_unique<pair_poly_t, pair_t>(single);
        }
        template<typename single_t>
        void set_pair(const single_t& single, const mbf_foreach::frm_bos::Base&) {
            typedef mbf_foreach::frm_bos::Pair<single_t> pair_t;
            m_pair = smart_ptr::make_poly_unique<pair_poly_t, pair_t>(single);
        }

    public:
        template<typename single_t>
        BasisIters(const single_t& single) {
            m_single = smart_ptr::make_poly_unique<single_poly_t, single_t>(single);
            /*
             * use overloading of the second parameter to statically dispatch to the correct ctor within either
             *  mbf_foreach::frm, mbf_foreach::bos, or mbf_foreach::frm_bos
             */
            set_pair<single_t>(single, single);
        }

        BasisIters(){}

        uint_t niter_single() const;

        /**
         * examine the requirements of the given Hamiltonian, and build an iterator over the row and column spaces
         * of its matrix representation in a physically appropriate, complete space of MBFs
         * @param h
         *  hamiltonian object
         * @param force_general
         *  if true, use the full combinatorial basis with no symmetries assumed or imposed (for testing)
         * @return
         *  a type-agnostic iterator which fills the dense matrix representation on a call to its loop method
         */
        static BasisIters make(const Hamiltonian& h, sys::Particles particles, bool force_general);

        static BasisIters make(const Hamiltonian& h);

        /**
         * fall-through case if no specific iterator is found or if general iterator is requested
         */
        static BasisIters make_general(const Hamiltonian &h, sys::Particles particles);
    };
}

#endif //M7_FCI_H
