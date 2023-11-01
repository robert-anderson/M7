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

    smuvi.collate();

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
            {{1,4,5,6,7,8}, {1,2,3,5,7,8}}
            // {{0,1,4,5,6,8}, {0,1,4,5,7,8}},
            // {{0,1,5,6,7,8}, {1,2,4,6,7,8}},
            // {{1,4,5,6,7,8}, {1,2,4,5,7,8}},
    };
    const v_t<std::pair<uintv_t, uintv_t>> insertions2 = {
            {{1,2,4,6,7,8}, {1,2,4,5,7,8}},
            {{1,2,4,6,7,8}, {1,2,3,5,7,8}}
            // {{1,2,4,5,7,8}, {0,3,4,6,7,8}},
            // {{1,2,4,5,7,8}, {1,2,4,6,7,8}},
            // {{0,1,4,5,7,8}, {0,1,2,3,4,5}},
    };
    const v_t<std::pair<uintv_t, uintv_t>> smuvi_keypairs = {
            {{1,4,5,6,7,8}, {1,2,4,6,7,8}}
            // {{0,1,4,5,6,8}, {1,2,4,5,7,8}},
            // {{0,1,5,6,7,8}, {1,2,4,5,7,8}},
            // {{1,4,5,6,7,8}, {0,1,4,5,7,8}},
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
    smuvi1.collate();
    smuvi2.collate();

    buffered::FrmOnvSpinChannel tmp_key1(nsite);
    buffered::FrmOnvSpinChannel tmp_key2(nsite);
    for (auto& key_pair: smuvi_keypairs) {
        tmp_key1 = key_pair.first;
        tmp_key2 = key_pair.second;
        uint_t counter = 0ul;
        smuvi1.foreach_common_value(tmp_key1, smuvi2, tmp_key2,
                                    [&](const field::FrmOnvSpinChannel &common_string){
                                        counter += 1ul;
                                    });
        ASSERT_EQ(counter, 2ul);
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

    smuvi.collate();

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
     * generate and store all the single-channel configurations for spin-conserving ONV enumeration
     */
    v_t<uintv_t> all_channel_setbits;
    {
        auto fn = [&](const uintv_t& inds){all_channel_setbits.push_back(inds);};
        basic_foreach::rtnd::Ordered<true, true> foreach(basis.m_nsite, nelec_per_channel);
        foreach.loop(fn);
    }

    const auto nelem = all_channel_setbits.size() * all_channel_setbits.size();
    const uint_t max_nentry = 5;
    const uint_t max_entry = 20;

    const auto elem_entry_counts = hash::in_range<uint_t>(0, nelem, 0, max_nentry);
    v_t<uintv_t> elem_entries;
    elem_entries.reserve(nelem);
    for (uint_t ielem = 0; ielem < nelem; ++ielem) {
        const auto nentry = elem_entry_counts[ielem];
        elem_entries.emplace_back(hash::unique_in_range<uint_t>(ielem, nentry, 0, max_entry));
    }

    v_t<std::pair<uintp_t, uint_t>> input_data;

    uint_t ielem = 0;
    for (uint_t ialpha = 0; ialpha < all_channel_setbits.size(); ++ialpha) {
        for (uint_t ibeta = 0; ibeta < all_channel_setbits.size(); ++ibeta) {
            for (const auto& entry: elem_entries[ielem]) {
                input_data.push_back({{ialpha, ibeta}, entry});
            }
            ++ielem;
        }
    }

    /*
     * shuffle the input data to arbitrary order
     */
    const auto nshuffle = input_data.size();
    for (uint_t ishuffle = 0; ishuffle < nshuffle; ++ishuffle) {
        auto shuffle_pair = hash::in_range(ishuffle, 2, 0, input_data.size());
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

    smuvi.collate();

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