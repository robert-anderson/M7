//
// Created by rja on 19/02/23.
//

#include "HfExcitHists.h"

hf_excit_hist::IndVals::IndVals(const hdf5::NodeReader &parent, str_t name, wf_comp_t thresh) :
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
    while(m_nelement < m_vals.nelement() && std::abs(m_vals[m_nelement]) >= thresh) ++m_nelement;
}

hf_excit_hist::Initializer::Initializer(
        wf::Fci &wf, const Mbf &hf, str_t fname, wf_comp_t thresh, uint_t max_power, bool cancellation) :
        m_wf(wf), m_hf(hf), m_c2(hdf5::FileReader(fname), "2200", thresh), m_cancellation(cancellation),
        m_work_mbf(wf.m_sector), m_work_conn(m_work_mbf), m_thresh(max_power ? find_first(max_power) : thresh),
        m_ncreated({2*max_power_by_thresh()+1}) {
    logging::info("Maximum-magnitude C2 coefficient: {}", m_c2.m_vals[0]);
    if (max_power) {
        DEBUG_ASSERT_EQ(max_power, max_power_by_thresh(), "max power specified should match match power computed");
        logging::info("Specified maximum relevant power of C2 ({}) implies a threshold of {}", max_power, m_thresh);
    }
    else {
        logging::info("Threshold of {} implies maximum relevant power of C2 is {}", m_thresh, max_power_by_thresh());
    }
}

uint_t hf_excit_hist::Initializer::max_power_by_thresh() {
    const auto limit = uint_t(m_wf.m_sector.m_frm.m_elecs) / 2;
    wf_comp_t product = 1.0;
    for (uint_t ipower = 0ul; ipower < limit; ++ipower) {
        product *= std::abs(m_c2.m_vals[ipower]);
        if (product < m_thresh) return ipower;
    }
    return limit;
}

wf_comp_t hf_excit_hist::Initializer::find_first(uint_t npower) {
    wf_t product = 1.0;
    m_work_mbf = m_hf;
    uint_t ientry = 0ul;

    auto mag = [&](uint_t ientry) {return std::abs(m_c2.m_vals[ientry]);};

    for (uint_t ipower = 1ul; ipower <= npower; ++ipower){
        for (; ientry < m_c2.nelement(); ++ientry) {
            if (apply(m_work_mbf, ientry)) {
                // the entry was successfully applied
                const bool geo_mean_with_next = (ipower == npower && (ientry + 1) < m_c2.nelement());
                /*
                 * to avoid any inconsistencies due to floating-point inequalities, let the last contribution be between
                 * the current entry and the next-largest one
                 */
                product *= geo_mean_with_next ? std::sqrt(mag(ientry) * mag(ientry+1)) : mag(ientry);
                // exit the loop over entries
                break;
            }
        }
    }
    return std::abs(product);
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

void hf_excit_hist::Initializer::setup(Mbf &mbf, uint_t imax, uint_t ipower, wf_t prev_product) {
    for (uint_t i = 0ul; i < imax; ++i) {
        auto product = prev_product * m_c2.m_vals[i];
        /*
         * do not continue with this branch of contributions if the prev_product is already smaller than the thresh,
         * since the values are sorted in descending order
         */
        if (std::abs(product) < m_thresh) return;
        if (apply(mbf, i)) {
            // i-th indices were successfully applied
            if (mpi::i_am(m_wf.m_dist.irank(mbf))) {
                // generated MBF belongs on this MPI rank
                auto& row = m_wf.m_store.lookup_or_insert(mbf);
                row.m_permanitiator.set();
                if (ipower==1) row.m_ref_conn.set();
                // permanitiators cannot be deleted - need to protect if not already
                if (!row.is_protected()) row.protect();
                ++m_ncreated.m_local[ipower * 2];
                if (m_cancellation) row.m_weight += (phase(mbf) ? -1 : 1) * product;
                // if not observing cancellation, no need to store the product
            }
            setup(mbf, i, ipower+1, product);
            undo(mbf, i);
        }
    }
}

void hf_excit_hist::Initializer::setup() {
    m_work_mbf = m_hf;
    setup(m_work_mbf, m_c2.nelement(), 1, 1.0);
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
                if (std::abs(row.m_weight[0]) < m_thresh) {
                    m_work_conn.connect(m_hf, row.m_mbf);
                    --after_cancellation[m_work_conn.m_cre.size()];
                    row.unprotect();
                    m_wf.m_store.erase(row.m_mbf);
                }
                else {
                    row.m_weight = 0.0;
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

void hf_excit_hist::initialize(wf::Fci &wf, const Mbf &hf, str_t fname, wf_comp_t thresh, uint_t max_power, bool cancellation) {
    Initializer(wf, hf, fname, thresh, max_power, cancellation).setup();
}

void hf_excit_hist::initialize(wf::Fci &wf, const Mbf &hf, const conf::CiPermanitiator &opts) {
    logging::info("Setting permanitiators based on CI data from \"{}\"", opts.m_path.m_value);
    initialize(wf, hf, opts.m_path, opts.m_thresh, opts.m_max_power, opts.m_cancellation);
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
    if (m_save_file_name.empty() && !m_chkpt_files)
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
        Accumulators(hf, opts.m_nexcits, opts.m_thresh, opts.m_save, opts.m_chkpt){}

void hf_excit_hist::Accumulators::add(const Mbf &mbf, wf_t weight) {
    if (!*this) return;
    if (!m_accum_epoch) return;
    m_work_conn.connect(m_hf->mbf(), mbf);
    const auto nexcit = m_work_conn.exsig().nfrm_cre();
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
        OpSig exsig({nexcit, nexcit}, {0, 0});
        table->save(nw, exsig.to_string(), true);
        ++table;
    }
}

void hf_excit_hist::Accumulators::save() const {
    if (m_save_file_name.empty()) return;
    hdf5::FileWriter fw(m_save_file_name);
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
