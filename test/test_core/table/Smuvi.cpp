//
// Created by Robert John Anderson on 13/09/2023.
//

#include "test_core/defs.h"
#include "M7_lib/table/Smuvi.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/communication/Communicator.h"

TEST(Smuvi, LookupKeysIndices) {
    const uint_t nsite = 6;
    using smuvi_t = Smuvi<field::FrmOnvSpinChannel, field::Number<uint_t>>;
    smuvi_t smuvi("test smuvi", field::FrmOnvSpinChannel(nullptr, nsite), field::Number<uint_t>(nullptr));
    buffered::FrmOnvSpinChannel tmp_key(nsite);
    const v_t<std::pair<uintv_t, uint_t>> insertions = {
        {{0, 3, 5}, 4},
        {{0, 2, 5}, 3},
        {{1, 3, 5}, 6},
        {{0, 3, 4}, 7},
        {{0, 3, 5}, 9},
        {{0, 1, 5}, 8},
        {{0, 3, 5}, 2},
    };

    buffered::Number<uint_t> tmp_val;
    for (auto& insertion: insertions) {
        tmp_key = insertion.first;
        tmp_val = insertion.second;
        smuvi.insert(tmp_key, tmp_val);
    }

    const auto order_fn = [&](const field::Number<uint_t>& i, const field::Number<uint_t>& j) -> bool {
        return i < j;
    };
    smuvi.collate(order_fn);

    auto fn = [&smuvi](const field::FrmOnvSpinChannel& key, smuvi_t::AccessResult values) -> void {
        auto lookup_values = smuvi.access(key);
        ASSERT_EQ(lookup_values, values);
    };
    smuvi.foreach_key(fn);
}


TEST(Smuvi, Intersection) {
    const uint_t nsite = 9;
    using smuvi_t = Smuvi<field::FrmOnvSpinChannel, field::FrmOnvSpinChannel>;
    smuvi_t smuvi1("test smuvi one", field::FrmOnvSpinChannel(nullptr, nsite), field::FrmOnvSpinChannel(nullptr, nsite));
    smuvi_t smuvi2("test smuvi two", field::FrmOnvSpinChannel(nullptr, nsite), field::FrmOnvSpinChannel(nullptr, nsite));

    const v_t<std::pair<uintv_t, uintv_t>> insertions1 = {
            {{1,4,5,6,7,8}, {1,2,4,5,7,8}},
            {{1,4,5,6,7,8}, {0,1,2,4,5,7}}
    };
    const v_t<std::pair<uintv_t, uintv_t>> insertions2 = {
            {{1,2,4,6,7,8}, {1,2,3,5,6,7}},
            {{1,2,4,5,7,8}, {1,2,4,6,7,8}},
            {{0,1,4,5,7,8}, {1,2,4,5,7,8}}
    };

    buffered::FrmOnvSpinChannel tmp_key(nsite);
    buffered::FrmOnvSpinChannel tmp_val(nsite);
    for (auto& insertion: insertions1) {
        tmp_key = insertion.first;
        tmp_val = insertion.second;
        smuvi1.insert(tmp_key, tmp_val);
    }
    for (auto& insertion: insertions2) {
        tmp_key = insertion.first;
        tmp_val = insertion.second;
        smuvi2.insert(tmp_key, tmp_val);
    }

    const auto order_fn = [&](const field::FrmOnvSpinChannel& i, const field::FrmOnvSpinChannel& j) -> bool {
        return i.lexical_order(j);
    };
    smuvi1.collate(order_fn);
    smuvi2.collate(order_fn);

    /*
     * Each element is a pair with:
     *  first element the pair of keys:
     *   first element the key to access smuvi1
     *   second element the key to access smuvi2
     *  second element the number of items in the intersection of the items between the two smuvi item vectors
     */
    const v_t<std::pair<std::pair<uintv_t, uintv_t>, uint_t>> smuvi_keypairs = {
            {{{1,4,5,6,7,8}, {1,2,4,6,7,8}}, 0},
            {{{1,4,5,6,7,8}, {1,2,4,5,7,8}}, 0},
            {{{1,4,5,6,7,8}, {0,1,4,5,7,8}}, 1} // {1,2,4,5,7,8}
    };
    buffered::FrmOnvSpinChannel tmp_key1(nsite);
    buffered::FrmOnvSpinChannel tmp_key2(nsite);
    for (auto& pair: smuvi_keypairs) {
        const auto key_pair = pair.first;
        tmp_key1 = key_pair.first;
        tmp_key2 = key_pair.second;
        uint_t counter = 0ul;
        smuvi1.foreach_common_value(tmp_key1, smuvi2, tmp_key2,
            [&](const field::FrmOnvSpinChannel &common_string){++counter;}, order_fn);
        ASSERT_EQ(counter, pair.second);
    }
}

TEST(Smuvi, BitsetToBitset) {
    const uint_t nsite = 6;
    using smuvi_t = Smuvi<field::FrmOnvSpinChannel, field::FrmOnvSpinChannel>;
    smuvi_t smuvi("test smuvi", field::FrmOnvSpinChannel(nullptr, nsite), field::FrmOnvSpinChannel(nullptr, nsite));
    buffered::FrmOnvSpinChannel tmp_key(nsite);
    buffered::FrmOnvSpinChannel tmp_val(nsite);

    tmp_key = uintv_t{0, 3, 5};
    tmp_val = uintv_t{0, 3, 4};

    if (mpi::irank() == 0) smuvi.insert(tmp_key, tmp_val);

    const auto order_fn = [&](const field::FrmOnvSpinChannel& i, const field::FrmOnvSpinChannel& j) -> bool {
        return i < j;
    };
    smuvi.collate(order_fn);

    auto fn = [&smuvi](const field::FrmOnvSpinChannel& key, smuvi_t::AccessResult values) -> void {
        auto lookup_values = smuvi.access(key);
        ASSERT_EQ(lookup_values, values);
    };
    smuvi.foreach_key(fn);
}

/**
 * @param input_data
 *  rank-private input data that is in general different on each rank
 * @return
 *  A copy of the all-gathered input data as a map that is the same on every rank
 */
std::map<uintp_t, uintv_t> make_global_data(const v_t<std::pair<uintp_t, uint_t>>& input_data) {
    using namespace integer;
    // pick a large enough number of "columns" in the rectangular map
    const size_t ncol = 1ul << (CHAR_BIT*sizeof(uint_t)/2);
    uintv_t local_keys;
    uintv_t local_values;
    for (auto& pair: input_data) {
        // encode the pair key as a single integer since MPI requires simple datatypes
        local_keys.push_back(rectmap(pair.first.first, pair.first.second, ncol));
        local_values.push_back(pair.second);
    }
    uintv_t global_keys;
    uintv_t global_values;
    mpi::all_gatherv(local_keys, global_keys);
    mpi::all_gatherv(local_values, global_values);
    REQUIRE_EQ(global_values.size(), global_keys.size(), "global keys and value data should be the same length");
    // assemble this global data as a map from pairs to sorted entries
    std::map<uintp_t, std::set<uint_t>> global_map;
    for (uint_t irow = 0; irow < global_keys.size(); ++irow) {
        uintp_t inds;
        // decode to a pair
        inv_rectmap(inds.first, inds.second, ncol, global_keys[irow]);
        auto it = global_map.find(inds);
        if (it == global_map.end()) it = global_map.insert({inds, {}}).first;
        it->second.insert(global_values[irow]);
    }
    // finally, assemble this global data as a map from pairs to sorted vector entries
    std::map<uintp_t, uintv_t> out;
    for (auto& pair: global_map) out.insert({pair.first, uintv_t(pair.second.cbegin(), pair.second.cend())});
    return out;
}

TEST(Smuvi, Comms) {
    const sys::frm::Basis basis(4);
    const uint_t nelec_per_channel = 2;

    /*
     * generate and store all the single-channel configurations for spin-conserving ONV enumeration on all ranks
     */
    v_t<uintv_t> all_channel_setbits;
    {
        auto fn = [&](const uintv_t& inds){all_channel_setbits.push_back(inds);};
        basic_foreach::rtnd::Ordered<true, true> foreach(basis.m_nsite, nelec_per_channel);
        foreach.loop(fn);
    }
    /*
     * the number of entries is the square of the spin channel enumeration length
     */
    const auto nelem = all_channel_setbits.size() * all_channel_setbits.size();
    /*
     * number of elements to draw as a fraction of the total number of elements
     */
    const double nelem_select_fraction = 0.3;
    /*
     * max number of uint_t entries to generate for each ONV
     */
    const uint_t max_nentry = 6;
    /*
     * max value of each uint_t entry
     */
    const uint_t max_entry = 10;

    /*
     * select a rank-specific subset of the elements
     */
    const auto select_elems = hash::unique_in_range(mpi::irank(), uint_t(nelem * nelem_select_fraction), 0, nelem, true);

    const auto nelem_select = select_elems.size();

    /*
     * now generate the rank-specific numbers of integer entries
     */
    const auto elem_entry_counts = hash::in_range<uint_t>(mpi::irank(), nelem_select, 0, max_nentry);
    v_t<uintv_t> elem_entries;
    elem_entries.reserve(nelem_select);
    /*
     * loop over the elements and generate the required number of entries
     */
    for (uint_t ielem = 0; ielem < nelem_select; ++ielem) {
        const auto nentry = elem_entry_counts[ielem];
        elem_entries.emplace_back(hash::unique_in_range<uint_t>(ielem, nentry, 0, max_entry));
    }
    /*
     * bring all the input together into a single object where each key is given as a pair of spin channel indices and
     * each value is an associated entry
     */
    v_t<std::pair<uintp_t, uint_t>> input_data;

    auto select_it = select_elems.cbegin();
    uint_t ielem = 0;
    for (uint_t ialpha = 0; ialpha < all_channel_setbits.size(); ++ialpha) {
        for (uint_t ibeta = 0; ibeta < all_channel_setbits.size(); ++ibeta) {
            // skip over unselected string pairs
            if (select_it != select_elems.cend() && ielem == *select_it) {
                const auto ielem_select = uint_t(std::distance(select_elems.cbegin(), select_it));
                for (const auto &entry: elem_entries[ielem_select]) {
                    input_data.emplace_back(uintp_t{ialpha, ibeta}, entry);
                }
                ++select_it;
            }
            ++ielem;
        }
    }
    ASSERT_EQ(ielem, nelem);
    ASSERT_EQ(uint_t(std::distance(select_elems.cbegin(), select_it)), nelem_select);

    const auto ordered_input_data = input_data;
    /*
     * shuffle the input data to arbitrary order
     */
    const auto nshuffle = input_data.size();
    for (uint_t ishuffle = 0; ishuffle < nshuffle; ++ishuffle) {
        auto shuffle_pair = hash::unique_in_range(ishuffle, 2, 0, input_data.size());
        std::swap(input_data[shuffle_pair[0]], input_data[shuffle_pair[1]]);
    }

    Smuvi<field::FrmOnv, field::Number<uint_t>> smuvi("my_smuvi", field::FrmOnv(nullptr, basis), field::Number<uint_t>(nullptr));
    buffered::FrmOnv mbf(basis);

    buffered::Number<uint_t> val;

    for (const auto& data: input_data) {
        const auto& alpha_string = all_channel_setbits[data.first.first];
        const auto& beta_string = all_channel_setbits[data.first.second];
        val = data.second;
        mbf = {alpha_string, beta_string};
        smuvi.insert(mbf, val);
    }

    const auto order_fn = [&](const field::Number<uint_t>& i, const field::Number<uint_t>& j) -> bool {
        return i < j;
    };
    smuvi.collate(order_fn);

    auto global_data = make_global_data(input_data);
    // look up all keys in the input
    for (uint_t irank=0; irank < mpi::nrank(); ++irank) {
        if (!mpi::i_am(irank)) continue;
        for (const auto &pair: global_data) {
            const auto &alpha_string = all_channel_setbits[pair.first.first];
            const auto &beta_string = all_channel_setbits[pair.first.second];
            mbf = {alpha_string, beta_string};
            auto access_result = smuvi.access(mbf);
            // firstly, length of found data should match that of the global data value-vector
            ASSERT_EQ(access_result.nremain(), pair.second.size());
            // then ensure that all values are the same and occur in the same order
            ASSERT_EQ(access_result.to_vector(), pair.second);
        }
    }
}