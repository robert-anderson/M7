//
// Created by rja on 26/01/23.
//

#ifndef M7_HFC2ACCUMULATION_H
#define M7_HFC2ACCUMULATION_H

#include "M7_lib/table/GlobalAccumulation.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/communication/SharedRows.h"
#include "M7_lib/parallel/Reduction.h"

namespace hf_excit_coeffs{
    struct Row : ::Row {
        field::MaeInds m_excit_inds;
        field::Numbers<wf_t, c_ndim_wf> m_weight;
        Row(NdFormat<c_ndim_wf> format):
            m_excit_inds(this, opsig::c_doub, "excit_indices"), m_weight(this, format){}
        field::MaeInds &key_field() {
            return m_excit_inds;
        };
        field::Numbers<wf_t, c_ndim_wf> &value_field() {
            return m_weight;
        };
    };

    class HfExcitCoeffs : public GlobalAccumulation<Row> {
        static const NdFormat<c_ndim_wf>& fmt(const shared_rows::Walker* hf) {
            return hf->gathered().m_row.m_weight.m_format;
        }
    public:
        Reduction<wf_t> m_coherent_c4_l1;
        Reduction<wf_t> m_total_c4_l1;
        const shared_rows::Walker* m_hf;
        mutable conn::Mbf m_work_conn;
        mutable buffered::MaeInds m_work_key;
        mutable buffered::Numbers<wf_t, c_ndim_wf> m_work_weight;
        HfExcitCoeffs(const shared_rows::Walker* hf):
            GlobalAccumulation<Row>("HF excit coeffs", {fmt(hf)}), m_hf(hf),
            m_work_conn(hf->mbf().m_basis), m_work_key(opsig::c_doub), m_work_weight(fmt(hf)){}

        void add(const Walker& walker) {
            m_work_conn.connect(m_hf->mbf(), walker.m_mbf);
            if (m_work_conn.exsig()!=opsig::c_doub) return;
            auto phase = m_work_conn.phase(m_hf->mbf());
            m_work_key = m_work_conn;
            m_work_weight = walker.m_weight;
            if (phase) m_work_weight*=-1.0;
            GlobalAccumulation<Row>::add(m_work_key, m_work_weight);
        }

    private:
        template<typename fn_t>
        static void part_predict_loop(const fn_t& fn) {
            functor::assert_prototype<void(uinta_t<4>, bool)>(fn);
            fn({0, 1, 2, 3}, false);
            fn({0, 2, 1, 3}, true);
            fn({0, 3, 1, 2}, false);
            fn({1, 2, 0, 3}, false);
            fn({1, 3, 0, 2}, true);
            fn({2, 3, 0, 1}, false);
        }

        template<typename fn_t>
        static void predict_loop(const fn_t& fn) {
            functor::assert_prototype<void(uinta_t<8>, bool)>(fn);
            auto fn_outer = [&](uinta_t<4> cre_inds, bool cre_par){
                uinta_t<8> inds;
                std::copy(cre_inds.cbegin(), cre_inds.cend(), inds.begin());
                auto fn_inner = [&](uinta_t<4> ann_inds, bool ann_par) {
                    auto par = (cre_par!=ann_par);
                    std::copy(ann_inds.cbegin(), ann_inds.cend(), inds.begin()+4);
                    fn(inds, par);
                };
                part_predict_loop(fn_inner);
            };
            part_predict_loop(fn_outer);
        }

    public:
        bool predict(const field::FrmOnv& mbf, wf_t& v) const {
            v = 0.0;
            m_work_conn.connect(m_hf->mbf(), mbf);
            const auto exsig = m_work_conn.exsig();
            if (exsig != opsig::c_quad) return false;
            auto& conn = m_work_conn;
            auto& key = m_work_key.m_frm;
            auto fn = [&](uinta_t<8> inds, bool par) {
                wf_t prod;
                // first pair of inds are creation
                key.m_cre[0] = conn.m_cre[inds[0]];
                key.m_cre[1] = conn.m_cre[inds[1]];
                // third pair of inds are annihilation
                key.m_ann[0] = conn.m_ann[inds[4]];
                key.m_ann[1] = conn.m_ann[inds[5]];
                {
                    auto lookup_row = current().lookup(m_work_key);
                    prod = lookup_row ? lookup_row.m_weight[0] : 0.0;
                }
                if (!prod) return;
                // second pair of inds are creation
                key.m_cre[0] = conn.m_cre[inds[2]];
                key.m_cre[1] = conn.m_cre[inds[3]];
                // last pair of inds are annihilation
                key.m_ann[0] = conn.m_ann[inds[6]];
                key.m_ann[1] = conn.m_ann[inds[7]];
                {
                    auto lookup_row = current().lookup(m_work_key);
                    prod *= lookup_row ? lookup_row.m_weight[0] : 0.0;
                }
                v += par ? -prod : prod;
            };
            predict_loop(fn);
            /*
             * in intermediate normalization:
             * c4 ~= 0.5 * c2^2
             * let nk be the number of walkers arrays for an excitation level k
             * v contains (n2 x naccum)^2 (naccum is number of samples)
             */
            v /= naccum() * naccum();
            /*
             * v now contains n2^2
             * ck = nk / n0
             * nk = n0 * ck
             *
             * n4 ~= 0.5 * n0 * (n2 / n0)^2
             * n4 ~= n2^2 / (2*n0)
             */
            v /= m_hf->weight(0);
            if (m_work_conn.phase(m_hf->mbf())) v *= -1;
            return true;
        }

        wf_t predict(const field::BosOnv&) const {
            return 0.0;
        }
        wf_t predict(const field::FrmBosOnv&) const {
            return 0.0;
        }

        bool is_initiator(const Walker& walker, double fac) {
            wf_t p;
            auto is_c4 = predict(walker.m_mbf, p);
            if (!is_c4) return false;
            const auto weight = walker.m_weight[0];
            m_total_c4_l1.m_local += std::abs(weight);
            if (!p) return false;
            const auto coherent = (p>0)==(weight>0);
            if (coherent) m_coherent_c4_l1.m_local += std::abs(weight);
            if (fac==0.0) {
                // with no fac set, assume we allow walkers of the correct predicted sign to become initiators
                return coherent;
            }
            else {
                // if the absolute predicted weight is larger than fac, let walker be initiator
                return std::abs(p) > fac;
            }
        }

        void update() {
            m_coherent_c4_l1.all_sum();
            m_total_c4_l1.all_sum();
            GlobalAccumulation<Row>::update();
        }

    };

}


#endif //M7_HFC2ACCUMULATION_H
