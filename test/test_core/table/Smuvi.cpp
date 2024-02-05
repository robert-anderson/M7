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

    const v_t<std::pair<std::pair<uintv_t, uintv_t>, uint_t>> smuvi_keypairs = {
            {{{1,4,5,6,7,8}, {1,2,4,6,7,8}}, 0},
            {{{1,4,5,6,7,8}, {1,2,4,5,7,8}}, 0},
            {{{1,4,5,6,7,8}, {0,1,4,5,7,8}}, 1}
    };
    buffered::FrmOnvSpinChannel tmp_key1(nsite);
    buffered::FrmOnvSpinChannel tmp_key2(nsite);
    for (auto& pair: smuvi_keypairs) {
        const auto key_pair = pair.first;
        tmp_key1 = key_pair.first;
        tmp_key2 = key_pair.second;
        uint_t counter = 0ul;
        smuvi1.foreach_value(tmp_key1, [&](const field::FrmOnvSpinChannel& val){std::cout << val << " ";});
        std::cout << std::endl;
        smuvi2.foreach_value(tmp_key2, [&](const field::FrmOnvSpinChannel& val){std::cout << val << " ";});
        std::cout << std::endl;
        smuvi1.foreach_common_value(tmp_key1, smuvi2, tmp_key2,
                                    [&](const field::FrmOnvSpinChannel &common_string){counter += 1ul;}, order_fn);
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

    smuvi.insert(tmp_key, tmp_val);

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
    const double nelem_select_fraction = 0.8;
    /*
     * max number of uint_t entries to generate for each ONV
     */
    const uint_t max_nentry = 5;
    /*
     * max value of each uint_t entry
     */
    const uint_t max_entry = 20;

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

    uintv_t values {};
    uintv_t alpha_channels {};
    uintv_t beta_channels {};
    for (const auto& data: input_data) {
        const auto& alpha_string = all_channel_setbits[data.first.first];
        const auto& beta_string = all_channel_setbits[data.first.second];
        logging::info_("val {} data: {}", mbf, val);
        val = data.second;
        mbf = {alpha_string, beta_string};
        smuvi.insert(mbf, val);

        alpha_channels.emplace_back(data.first.first);
        beta_channels.emplace_back(data.first.second);
        values.emplace_back(val);
    }

    uintv_t values_global {};
    uintv_t alpha_channels_global {};
    uintv_t beta_channels_global {};
    mpi::all_gatherv(values, values_global);
    mpi::all_gatherv(alpha_channels, alpha_channels_global);
    mpi::all_gatherv(beta_channels, beta_channels_global);
    v_t<buffered::FrmOnv> mbfs_global {};
    for (uint_t i = 0; i < alpha_channels_global.size(); ++i) {
        const auto& alpha_string = all_channel_setbits[alpha_channels_global[i]];
        const auto& beta_string  = all_channel_setbits[beta_channels_global[i]];
        mbf = {alpha_string, beta_string};
        mbfs_global.emplace_back(mbf);
    }

    const auto order_fn = [&](const field::Number<uint_t>& i, const field::Number<uint_t>& j) -> bool {
        return i < j;
    };
    smuvi.collate(order_fn);


    // look up all keys in the input
    for (const auto& data: ordered_input_data) {
        const auto& alpha_string = all_channel_setbits[data.first.first];
        const auto& beta_string = all_channel_setbits[data.first.second];
        mbf = {alpha_string, beta_string};
        auto access_result = smuvi.access(mbf);
    }

#if 0
    std::cout << mbf << std::endl;
    std::cout << smuvi.nitem(mbf) << std::endl;
    std::cout << "++++++++++" << std::endl;
    for (uint_t iitem = 0; iitem < smuvi.nitem(mbf); ++iitem){
        std::cout << smuvi.item_cbegin(mbf)[iitem] << std::endl;
    }

    if (mpi::i_am_root()) {
//        std::cout << smuvi.m_inserter.m_send_recv.recv().to_string() << std::endl;
        std::cout << convert::to_string(mpi::g_nrank_in_shmem_realms) << std::endl;
        std::cout << convert::to_string(mpi::g_irank_root_in_shmem_realms) << std::endl;
    }
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) {
        if (mpi::i_am(irank)) {
            std::cout << mpi::irank() << " " << mpi::irank(mpi::SharedMemory) << std::endl;
        }
        mpi::barrier();
    }
#endif
}