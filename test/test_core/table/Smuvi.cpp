//
// Created by Robert John Anderson on 13/09/2023.
//

#include "test_core/defs.h"
#include "M7_lib/table/Smuvi.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/communication/Communicator.h"
#include "M7_lib/foreach/BasicForeach.h"

TEST(Smuvi, LocalTable) {
    typedef SingleFieldRow<field::Number<uint_t>> row_t;
    buffered::smuvi::LocalTable<row_t> table("test", row_t(), {});
    buffered::Number<uint_t> work;
    v_t<uintp_t> pairs = {
        {5, 11}, {3, 15}, {4, 19}, {3, 18}, {4, 18}, {5, 15}, {5, 14}, {3, 17}, {4, 17}
    };
    std::map<uint_t, uintv_t> combined_pairs;
    for (const auto& pair: pairs) {
        auto it = combined_pairs.find(pair.first);
        if (it==combined_pairs.end()) it = combined_pairs.insert({pair.first, {}}).first;
        it->second.push_back(pair.second);
    }

    for (const auto& pair : pairs) {
        work = pair.first;
        table.append(work, pair.second);
    }

    auto row = table.m_row;
    for (row.restart(); row; ++row) {
        const uint_t key = row.m_field;
        auto inds_ptr = table.inds(row);
        ASSERT_TRUE(inds_ptr);
        ASSERT_EQ(*inds_ptr, combined_pairs[key]);
    }
}


template<typename key_t>
struct Smuvi {

    struct InsertRow : Row {
        key_t m_key;
        field::Number<uint_t> m_entry;

        template<typename ...Args>
        InsertRow(const Args&... key_ctor_args): m_key(this, key_ctor_args...), m_entry(this){}

        key_t &key_field() {
            return m_key;
        };
    };

    communicator::BasicSend<InsertRow, InsertRow> m_inserter;
    v_t<buffered::MappedTable<InsertRow>> m_accessors;

    template<typename ...Args>
    Smuvi(str_t name, const Args&... key_ctor_args):
        m_inserter(
            name,
            InsertRow(key_ctor_args...),
            DistribOptions(),
            {100, 1.0},
            InsertRow(key_ctor_args...),
            {100, 1.0}){

//        for (uint_t irank=0ul; irank < mpi::nrank(); ++irank) {
//
//        }
//        m_accessors(name, InsertRow(key_ctor_args...), true) {}
//
//        m_accessors(name, InsertRow(key_ctor_args...), true)

    }



    void insert(const key_t& key, const uint_t& entry) {
        const auto irank_dst = m_inserter.m_dist.irank(key);
        Table<InsertRow>& send = m_inserter.m_send_recv.send(irank_dst);
        send.m_row.push_back_jump();
        send.m_row.m_key = key;
        send.m_row.m_entry = entry;
    }

    void collate() {
        m_inserter.communicate();
//        for (auto&)
    }



};
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

    std::cout << smuvi.m_inserter.m_send_recv.recv().to_string() << std::endl;
}