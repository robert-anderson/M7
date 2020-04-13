//
// Created by rja on 29/02/2020.
//

#ifndef M7_DETERMINANTSAMPLER_H
#define M7_DETERMINANTSAMPLER_H


#include "HeatBathSampler.h"

class alignas(defs::cache_line_size) DeterminantSampler {
public:
    enum Outcome{no_excitations, single_excitation, double_excitation, both_excitations};
private:
    const HeatBathSampler &m_precomputed;
    PRNG &m_prng;
    Determinant m_src_det;
    AntisymConnection m_anticonn;
    OccupiedOrbitals m_occ;
    VacantOrbitals m_vac;
    const size_t &m_nspinorb;
    std::vector<defs::prob_t> m_P1;
    Aliaser m_P1_aliaser;
    Aliaser m_P2_aliaser;
    Aliaser m_P3_aliaser;
    Aliaser m_P4_aliaser;
    std::vector<defs::prob_t> m_P2_qp; // prob of picking q given p picked first
    std::vector<defs::prob_t> m_P2_pq; // prob of picking p given q picked first

    AntisymConnection m_single_excitation, m_double_excitation;
    Determinant m_single_dst_det, m_double_dst_det;
    defs::prob_t m_single_prob, m_double_prob;

    Outcome m_outcome = no_excitations;

public:

    DeterminantSampler(const HeatBathSampler &precomputed);

    void update(const DeterminantElement& det);

    void set_P1(std::vector<defs::prob_t> &P1);

    void set_P2(std::vector<defs::prob_t> &P2, const size_t &p);

    void draw_pq(size_t &p, size_t &q);

    void draw_pq(size_t &p, size_t &q, defs::prob_t &prob);

    void draw_r(const size_t &p, const size_t &q, size_t &r);

    void draw_r(const size_t &p, const size_t &q, size_t &r, defs::prob_t &prob);

    void draw_pqr(size_t &p, size_t &q, size_t &r);

    void draw_pqr(size_t &p, size_t &q, size_t &r, defs::prob_t &prob);

    void draw_s(const size_t &p, const size_t &q, const size_t &r, size_t &s);

    void draw_s(const size_t &p, const size_t &q, const size_t &r, size_t &s, defs::prob_t &prob);

    void draw(size_t &p, size_t &q, size_t &r, size_t &s,
              defs::prob_t &prob_single, defs::prob_t &prob_double,
              defs::ham_t &helement_single, defs::ham_t &helement_double);

    void draw();

    bool single_generated() const {return m_outcome==single_excitation || m_outcome==both_excitations;}
    bool double_generated() const {return m_outcome==double_excitation || m_outcome==both_excitations;}

    AntisymConnection& get_single() {
        assert(single_generated());
        return m_single_excitation;
    }

    AntisymConnection& get_double() {
        assert(double_generated());
        return m_double_excitation;
    }

    const Determinant &get_single_dst_det(){
        assert(single_generated());
        return m_single_dst_det;
    }

    const Determinant &get_double_dst_det(){
        assert(double_generated());
        return m_double_dst_det;
    }

    const defs::prob_t& get_single_prob() const {return m_single_prob;}
    const defs::prob_t& get_double_prob() const {return m_double_prob;}

private:
    defs::prob_t proposal(const size_t &ip, const size_t &r, const defs::ham_t &helement_single);

    defs::prob_t proposal(const size_t &ip, const size_t &iq,
                          const size_t &r, const size_t &s,
                          const defs::ham_t &h);
};


#endif //M7_DETERMINANTSAMPLER_H
