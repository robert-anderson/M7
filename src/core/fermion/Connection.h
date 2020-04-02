//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_CONNECTION_H
#define M7_CONNECTION_H

#include <src/core/fermion/Determinant.h>

/*
 * <bra| ...cre... ...des... |ket>
 * the m_des array contains the spin orbital indices in the ket but not the bra
 * the m_cre array contains the spin orbital indices in the bra but not the ket
 */

struct Connection {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::inds m_des, m_cre;
    size_t m_ndes, m_ncre;

    Connection(const DeterminantElement &ket, const DeterminantElement &bra) :
        m_nbit(ket.nbit()), m_element_dsize(ket.dsize()),
        m_des(m_nbit), m_cre(m_nbit) {
        assert(ket.compatible_with(bra));
        update_connection(ket, bra);
    }

    virtual void update(const DeterminantElement &ket, const DeterminantElement &bra) {
        update_connection(ket, bra);
    }

protected:
    void update_connection(const DeterminantElement &ket, const DeterminantElement &bra) {
        assert(ket.nbit() == m_nbit);
        assert(ket.dsize() == m_element_dsize);
        assert(ket.compatible_with(bra));
        m_ndes = 0ul;
        m_ncre = 0ul;

        DeterminantElement::DatawordEnumerator ket_enumerator(ket);
        DeterminantElement::DatawordEnumerator bra_enumerator(bra);
        defs::data_t ket_work, bra_work, work;
        size_t idataword = ~0ul;
        while(ket_enumerator.next(ket_work, idataword) && bra_enumerator.next(bra_work)){
            work = ket_work &~ bra_work;
            while (work) m_des[m_ndes++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
            work = bra_work &~ ket_work;
            while (work) m_cre[m_ncre++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
        }
    }
};

/*
 * a connection in which the common indices and antisymmetric phase is computed
 */
struct AntisymConnection : Connection {
    defs::inds m_com;
    size_t m_ncom;
    bool m_phase;

    AntisymConnection(const DeterminantElement &ket, const DeterminantElement &bra) :
        Connection(ket, bra), m_com(m_nbit) {
        update_antisym(ket, bra);
    }

    void update(const DeterminantElement &ket, const DeterminantElement &bra) override {
        update_connection(ket, bra);
        update_antisym(ket, bra);
    }

protected:
    void update_antisym(const DeterminantElement &ket, const DeterminantElement &bra) {
        m_ncom = 0ul;
        size_t nperm = 0ul;

        DeterminantElement::DatawordEnumerator ket_enumerator(ket);
        DeterminantElement::DatawordEnumerator bra_enumerator(bra);
        defs::data_t ket_work, bra_work, work;
        size_t idataword = ~0ul;

        auto des_iter = m_des.begin();
        const auto des_end = m_des.begin()+m_ndes;
        auto cre_iter = m_cre.begin();
        const auto cre_end = m_cre.begin()+m_ncre;

        while(ket_enumerator.next(ket_work, idataword) && bra_enumerator.next(bra_work)){
            work = ket_work & bra_work;
            while (work) {
                auto &com = m_com[m_ncom];
                com = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
                if (des_iter!=des_end && *des_iter > com){des_iter++; nperm+=m_ncom-1;}
                if (cre_iter!=cre_end && *cre_iter > com){cre_iter++; nperm+=m_ncom-1;}
                m_ncom++;
            }
        }
        while (des_iter!=des_end) {des_iter++; nperm+=m_ncom;}
        while (cre_iter!=cre_end) {cre_iter++; nperm+=m_ncom;}
        m_phase = nperm&1ul;
    }
};

#endif //M7_CONNECTION_H
