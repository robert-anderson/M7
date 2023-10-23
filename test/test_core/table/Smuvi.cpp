//
// Created by Robert John Anderson on 13/09/2023.
//

#include "test_core/defs.h"
#include "M7_lib/table/Smuvi.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/communication/Communicator.h"

TEST(Smuvi, LookupKeysIndices) {
    const uint_t nsite = 6;
    Smuvi<field::FrmOnvSpinChannel> smuvi("test smuvi", nsite);
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

    for (auto& insertion: insertions) {
        tmp_key = insertion.first;
        smuvi.insert(tmp_key, insertion.second);
    }

    smuvi.collate();

    for (auto& insertion: insertions) {
        tmp_key = insertion.first;
        auto key_lookup_result = smuvi.access_by_key(tmp_key);
        auto index_lookup_result = smuvi.access_by_index(key_lookup_result.m_key_index);
        ASSERT_EQ(index_lookup_result.m_key_row.m_key, tmp_key);
    }

    auto fn = [&](uint_t index, const Smuvi<field::FrmOnvSpinChannel>::SmuviEntriesWithKey& entries) {
        auto index_chk = smuvi.access_by_key(entries.m_key_row.m_key).m_key_index;
        ASSERT_EQ(index, index_chk);
    };
    smuvi.foreach(fn);

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

    Smuvi<field::FrmOnv> smuvi("my_smuvi", basis);
    buffered::FrmOnv mbf(basis);

    for (const auto& data: input_data) {
        const auto& alpha_string = all_channel_setbits[data.first.first];
        const auto& beta_string = all_channel_setbits[data.first.second];
        const auto& entry = data.second;
        mbf = {alpha_string, beta_string};
        smuvi.insert(mbf, entry);
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