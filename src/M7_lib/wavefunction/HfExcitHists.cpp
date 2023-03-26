//
// Created by rja on 19/02/23.
//

#include <M7_lib/util/Math.h>
#include "HfExcitHists.h"
#include "Wavefunction.h"

hf_excit_hist::IndVals::IndVals(const hdf5::NodeReader &parent, str_t name, double geo_mean_power_thresh) :
        m_inds(hdf5::GroupReader(parent, name), "indices", mpi::on_node_i_am_root(), true),
        m_vals(hdf5::GroupReader(parent, name), "values", mpi::on_node_i_am_root(), true) {
    REQUIRE_EQ(m_inds.nrow(), m_vals.nelement(),
               "number of index arrays is not the same as the number of values");
    uintv_t order;
    // intermediate normalization
    m_vals /= hdf5::DatasetLoader::load_vector<wf_t>(parent, "norm")[0];
    if (mpi::on_node_i_am_root()) m_vals.sort_inds(order, false, true);
    m_inds.reorder_rows(order);
    m_vals.reorder(order);
    mpi::barrier_on_node();
    m_geo_mean = math::geo_mean(m_vals.tbegin(), m_vals.nelement());
    m_thresh = std::pow(m_geo_mean, geo_mean_power_thresh);
    while(m_nelement < m_vals.nelement() && std::abs(m_vals[m_nelement]) >= m_thresh) ++m_nelement;
}

hf_excit_hist::Initializer::Initializer(
        wf::Vectors &wf, const Mbf &hf, str_t fname, double geo_mean_power_thresh, bool cancellation) :
        m_wf(wf), m_hf(hf), m_c2(hdf5::FileReader(fname), "2", geo_mean_power_thresh),
        m_cancellation(cancellation), m_work_mbf(wf.m_sector), m_work_conn(m_work_mbf),
        m_ncreated({2*max_power_by_thresh()+1}) {
    logging::info("Maximum-magnitude C2 coefficient: {}", m_c2.m_vals[0]);
    logging::info("Geometric mean of C2 coefficients: {}", m_c2.m_geo_mean);
    logging::info("Specified power of C2 geo mean ({}) implies a threshold of {} for C2^n permanitiators",
                  geo_mean_power_thresh, m_c2.m_thresh);
    logging::info("Threshold implies maximum relevant power of C2 is {}", max_power_by_thresh());
}

uint_t hf_excit_hist::Initializer::max_power_by_thresh() {
    const auto limit = uint_t(m_wf.m_sector.m_frm.m_elecs) / 2;
    wf_comp_t product = 1.0;
    for (uint_t ipower = 0ul; ipower < limit; ++ipower) {
        product *= std::abs(m_c2.m_vals[ipower]);
        product /= (ipower+1);
        if (product < m_c2.m_thresh) return ipower;
    }
    return limit;
}

bool hf_excit_hist::Initializer::apply(Mbf &mbf, uint_t ientry) {
    const auto a = m_c2.m_inds(ientry, 0);
    const auto b = m_c2.m_inds(ientry, 1);
    const auto i = m_c2.m_inds(ientry, 2);
    const auto j = m_c2.m_inds(ientry, 3);
    if (mbf.get(a) || mbf.get(b) || !mbf.get(i) || !mbf.get(j)) return false;
    mbf.set(a);
    mbf.set(b);
    mbf.clr(i);
    mbf.clr(j);
    return true;
}

bool hf_excit_hist::Initializer::undo(Mbf &mbf, uint_t ientry) {
    const auto a = m_c2.m_inds(ientry, 0);
    const auto b = m_c2.m_inds(ientry, 1);
    const auto i = m_c2.m_inds(ientry, 2);
    const auto j = m_c2.m_inds(ientry, 3);
    if (!mbf.get(a) || !mbf.get(b) || mbf.get(i) || mbf.get(j)) return false;
    mbf.clr(a);
    mbf.clr(b);
    mbf.set(i);
    mbf.set(j);
    return true;
}

bool hf_excit_hist::Initializer::phase(Mbf &mbf) {
    m_work_conn.connect(m_hf, mbf);
    return m_work_conn.phase(m_hf);
}

void hf_excit_hist::Initializer::setup(Mbf& mbf, uint_t imax, uint_t ipower, wf_t prev_product) {
    for (uint_t i = 0ul; i < imax; ++i) {
        if (loop_body(mbf, i, ipower, prev_product)) return;
    }
}

void hf_excit_hist::Initializer::setup(Mbf& mbf) {
    // parallelize over the first loop in the recursion
    const uint_t nelement_local = mpi::evenly_shared_count(m_c2.m_nelement);
    /*
     * because there is a synchronization at the end of each iteration of the top-level loop, each rank must execute
     * reach the communication call same number of times.
     */
    const auto nelement_local_max = mpi::all_max(nelement_local);
    const auto stride = mpi::nrank();
    const auto offset = mpi::irank();

    for (uint_t i = 0ul; i < nelement_local_max; ++i) {
        // only call the loop body if this rank still has work to do
        if (i < nelement_local) {
            const auto ielement = offset + i * stride;
            loop_body(mbf, ielement, 1, 1.0);
        }
        communicate_and_insert();
    }
}

void hf_excit_hist::Initializer::communicate_and_insert() {
    m_wf.m_send_recv.communicate();
    auto& store_row = m_wf.m_store.m_row;
    auto fn = [&](const Spawn& recv_row) {
        m_wf.m_store.lookup(recv_row.m_dst_mbf, store_row);
        if (!store_row) {
            // MBF not already added
            m_wf.m_store.insert(recv_row.m_dst_mbf, store_row);
            store_row.m_permanitiator.set();
            // permanitiators cannot be deleted - need to protect
            store_row.protect();
            const auto exlvl = OpCounts(m_hf, recv_row.m_dst_mbf).m_nfrm_cre;
            DEBUG_ASSERT_FALSE(exlvl & 1ul, "excitation level should be even, since currently only C2 is used as source info");
            const auto ipower = exlvl / 2;
            ++m_ncreated.m_local[ipower * 2];
        }
        store_row.m_weight += wf_t(recv_row.m_delta_weight);
    };
    m_wf.recv().foreach_row_in_use(fn);
}

void hf_excit_hist::Initializer::setup() {
    m_work_mbf = m_hf;
    setup(m_work_mbf);
    m_ncreated.all_sum();
    auto& row = m_wf.m_store.m_row;
    strv_t header = {"excitation level", "number in use"};
    v_t<strv_t> logging_table;
    auto max_power = this->max_power_by_thresh();
    if (m_cancellation) {
        header.emplace_back("number revoked by cancellation");
        /*
         * the present state of the wavefunction m_weight stores the cumulative value of the CI products
         * loop over the WF store table and delete any permanitiators which have fallen beneath thresh due to cancellation
         */
        auto after_cancellation = m_ncreated.m_reduced;
        for (row.restart(); row; ++row) {
            if (row.m_permanitiator.get(0)) {
                const auto ipower = OpCounts(m_hf, row.m_mbf).m_nfrm_cre / 2;
                if (std::abs(row.m_weight[0]) < m_c2.m_thresh) {
                    // the reference is not eligible for deletion
                    if (ipower) {
                        --after_cancellation[ipower];
                        row.unprotect();
                        m_wf.m_store.erase(row.m_mbf);
                    }
                }
                else {
                    if (ipower) row.m_weight = 0.0;
                }
            }
        }

        logging_table.push_back(header);
        for (uint_t i = 1ul; i<=max_power*2; ++i) {
            logging_table.push_back({
                                            convert::to_string(i),
                                            convert::to_string(after_cancellation[i]),
                                            convert::to_string(m_ncreated.m_reduced[i] - after_cancellation[i])
                                    });
        }
    }
    else {
        logging_table.push_back(header);
        for (uint_t i = 1ul; i <= max_power * 2; ++i)
            logging_table.push_back({convert::to_string(i), convert::to_string(m_ncreated.m_reduced[i])});
    }
    logging::info_table("Permanitiator breakdown", logging_table, true);
}

bool hf_excit_hist::Initializer::loop_body(Mbf& mbf, uint_t ielement, uint_t ipower, wf_t prev_product) {
    auto product = prev_product * m_c2.m_vals[ielement];
    // the form of the product is (1/n!)*(c2)^n
    product /= ipower;
    /*
     * do not continue with this branch of contributions if the prev_product is already smaller than the thresh,
     * since the values are sorted in descending order
     */
    if (std::abs(product) < m_c2.m_thresh) return true;
    if (apply(mbf, ielement)) {
        // i-th indices were successfully applied
        auto irank = m_wf.m_dist.irank(mbf);
        auto& send_table = m_wf.send(irank);
        auto& send_row = send_table.m_row;
        send_row.push_back_jump();
        send_row.m_dst_mbf = mbf;
        if (m_cancellation) send_row.m_delta_weight = (phase(mbf) ? -1 : 1) * product;
        setup(mbf, ielement, ipower + 1, product);
        undo(mbf, ielement);
    }
    return false;
}

void hf_excit_hist::initialize(wf::Vectors &wf, const Mbf &hf, str_t fname, double geo_mean_power_thresh, bool cancellation) {
    Initializer(wf, hf, fname, geo_mean_power_thresh, cancellation).setup();
}

void hf_excit_hist::initialize(wf::Vectors &wf, const Mbf &hf, const conf::CiPermanitiator &opts) {
    logging::info("Setting permanitiators based on CI data from \"{}\"", opts.m_path.m_value);
    initialize(wf, hf, opts.m_path, opts.m_geo_mean_power_thresh, opts.m_cancellation);
}

uintv_t hf_excit_hist::Accumulators::make_nexcit_is_accumulated(const uintv_t &nexcits) {
    if (nexcits.empty()) return {};
    auto nmax = *std::max_element(nexcits.cbegin(), nexcits.cend());
    uintv_t out(nmax+1, ~0ul);
    uint_t i = 0;
    for (auto& n: nexcits) out[n] = i++;
    return out;
}

uint_t hf_excit_hist::Accumulators::ind(uint_t nexcit) const {
    return nexcit < m_accumulated_nexcit_inds.size() ? m_accumulated_nexcit_inds[nexcit] : ~0ul;
}

hf_excit_hist::Accumulators::Accumulators(
        const shared_rows::Walker *hf, uintv_t nexcits, wf_comp_t thresh,
        const conf::OptionalFile& save_file, const conf::OptionalFileSeries& chkpt_files) :
        m_hf(hf), m_thresh(thresh), m_work_conn(hf->mbf()),
        m_nexcits(std::move(nexcits)), m_accumulated_nexcit_inds(make_nexcit_is_accumulated(m_nexcits)),
        m_accum_epoch("HF excitation accumulation"),
        m_save_file_path(save_file.path_if_enabled()), m_chkpt_files(m_accum_epoch, chkpt_files){
    if (m_nexcits.empty()) return;
    REQUIRE_TRUE(m_hf, "HF state must be defined for HF excitation accumulation");
    m_tables.reserve(m_nexcits.size());
    m_lookup_keys.reserve(m_nexcits.size());
    if (m_save_file_path.empty() && !m_chkpt_files)
        logging::warn("Accumulating HF excitations of rank {}, but saves and checkpointing are both disabled",
                      convert::to_string(m_nexcits));
    for (auto& nexcit: m_nexcits){
        const OpSig exsig({nexcit, nexcit}, {0, 0});
        const auto name = exsig.to_string()+" excitations of HF state";
        m_tables.emplace_back(name, RdmRow(exsig, 1), false);
        m_tables.back().resize(500);
        m_tables.back().set_expansion_factor(2.0);
        m_lookup_keys.emplace_back(exsig);
    }
}

hf_excit_hist::Accumulators::Accumulators(const conf::HfExcits &opts, const shared_rows::Walker *hf) :
        Accumulators(hf, opts.m_nexcits, opts.m_thresh, opts.m_save, opts.m_chkpt){
}

void hf_excit_hist::Accumulators::add(const Mbf &mbf, wf_t weight) {
    if (!*this) return;
    if (!m_accum_epoch) return;
    m_work_conn.connect(m_hf->mbf(), mbf);
    const auto nexcit = conn::nfrm_cre(m_work_conn);
    if (!nexcit) {
        m_norm += weight;
        return;
    }
    const auto i = ind(nexcit);
    if (i == ~0ul) return;
    auto& table = m_tables[i];
    auto& key = m_lookup_keys[i];
    key = m_work_conn;
    auto& lookup = table.lookup(key);
    // if the excitation is already in the table, it is added regardless of current weight
    if (lookup) lookup.m_values[0] += weight;
    // if it's not already histogrammed, it can be added only if the instantaneous weight is sufficient
    else if (std::abs(weight) >= std::abs(m_hf->weight(0)*m_thresh)) {
        auto& insert_row = table.insert(key);
        insert_row.m_values += weight;
    }
}

void hf_excit_hist::Accumulators::save(const hdf5::NodeWriter &nw) const {
    if (!m_accum_epoch) return;
    auto table = m_tables.begin();
    auto norm = mpi::all_sum(m_norm);
    hdf5::DatasetSaver::save_scalar(nw, "norm", norm);
    for (auto nexcit: m_nexcits) {
        /*
         * HF-related quantities are only applicable to fermionic excitations, so we may dispense with the full wxyz
         * string naming convention designed to support fermion-boson excitations
         */
        table->save(nw, convert::to_string(nexcit), true);
        ++table;
    }
}

void hf_excit_hist::Accumulators::save() const {
    if (m_save_file_path.empty()) return;
    hdf5::FileWriter fw(m_save_file_path);
    save(fw);
}

void hf_excit_hist::Accumulators::attempt_chkpt(uint_t icycle) {
    if (!m_chkpt_files) return;
    auto path = m_chkpt_files.get_file_path(icycle);
    if (path.empty()) return;
    logging::info("Saving HF excitation accumulators to checkpoint file {} on cycle {}", path, icycle);
    hdf5::FileWriter fw(path);
    fw.save_attr("icycle", icycle);
    save(fw);
}
