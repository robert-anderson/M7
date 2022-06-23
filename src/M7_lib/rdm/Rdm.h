//
// Created by Robert J. Anderson on 10/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include <M7_lib/io/Archivable.h>
#include <M7_lib/table/Communicator.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/bilinear/FermionPromoter.h>


#if 0
struct Rdm {};

struct FermionRdm : Communicator<MaeRow<defs::wf_t>, MaeRow<defs::wf_t>, true>, Archivable {
    typedef Communicator<MaeRow<defs::wf_t>, MaeRow<defs::wf_t>, true> base_t;
    const size_t m_nann, m_ncre, m_nelec;
    /**
     * working indices for building promotions and looking up the MEV tables
     */
    buffered::FermionMevInds m_lookup_inds;
    /**
     * all promoters required for an RDM of this rank. Index refers to the number of SQ ops
     * being inserted
     */
    std::vector<FermionPromoter> m_promoters;

    conn::FrmOnv m_conn;
    const bool m_mixed_estimator;
    mutable FrmOps m_com;

    const size_t &nop() const;

    static size_t nrow_estimate(size_t nann, size_t ncre, size_t nsite);

    FermionRdm(const conf::FermionRdm &opts, size_t nrow_crude_est, size_t nsite, size_t nelec);
    FermionRdm(const conf::FermionRdm &opts, size_t nsite, size_t nelec):
    FermionRdm(opts, nrow_estimate(opts.m_rank, opts.m_rank, nsite), nsite, nelec){}

    void make_contribs(const fields::FrmOnv &src_onv, const conn::FrmOnv &conn, const FrmOps &com,
                       const defs::wf_t &src_weight, const defs::wf_t &dst_weight);

    void make_contribs(const fields::FrmOnv &src_onv, const defs::wf_t &src_weight,
                       const fields::FrmOnv &dst_onv, const defs::wf_t &dst_weight) {
        m_conn.connect(src_onv, dst_onv, m_com);
        make_contribs(src_onv, m_conn, m_com, src_weight, dst_weight);
    }

    void make_contribs(const fields::FrmBosOnv &src_onv, const defs::wf_t &src_weight,
                       const fields::FrmBosOnv &dst_onv, const defs::wf_t &dst_weight) {
        make_contribs(src_onv.m_frm, src_weight, dst_onv.m_frm, dst_weight);
    }

    void make_contribs(const fields::FrmOnv &src_onv, const defs::wf_t &src_weight,
                       const fields::FrmOnv &dst_onv, const defs::wf_t &dst_weight, const size_t &nop_conn) {
        m_conn.connect(src_onv, dst_onv, m_com);
        if (m_conn.size() != nop_conn) return;
        make_contribs(src_onv, m_conn, m_com, src_weight, dst_weight);
    }

    void end_cycle() {
        if (!send().buffer_size()) return;
        communicate();
        auto &row = m_comm.recv().m_row;
        if (!m_comm.recv().m_hwm) return;
        for (row.restart(); row.in_range(); row.step()) {
            auto irow_store = *m_store[row.m_inds];
            if (irow_store == ~0ul) irow_store = m_store.insert(row.m_inds);
            m_store.m_row.jump(irow_store);
            m_store.m_row.m_values += row.m_values;
        }
        m_comm.recv().clear();
    }

    void h5_read(hdf5::GroupReader &parent) {
        m_store.clear();
        BufferedTable<MevRow<defs::wf_t>> m_buffer("", {{m_nann, m_ncre}});
        m_buffer.push_back();
        RowHdf5Reader<MevRow<defs::wf_t>> row_reader(m_buffer.m_row, parent, std::to_string(nop()), h5_field_names());

        row_reader.restart();
        for (size_t iitem = 0ul; iitem < row_reader.m_nitem; ++iitem) {
            row_reader.read(iitem);
            auto &send_table = send(m_ra.get_rank(row_reader.m_inds));
            // should never read in the same inds_t twice
            ASSERT(!send_table[row_reader.m_inds]);
            auto irow = send_table.insert(row_reader.m_inds);
            send_table.m_row.jump(irow);
            send_table.m_row.m_values = row_reader.m_values;
        }
    }

    void h5_write(hdf5::GroupWriter &parent) {
        m_store.save(parent, std::to_string(nop()), h5_field_names());
    }

    std::vector<std::string> h5_field_names() {
        return {m_store.m_row.m_inds.m_ann.m_name,
                m_store.m_row.m_inds.m_cre.m_name,
                m_store.m_row.m_values.m_name};
    }

protected:
    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;

public:


};

struct RdmGroup {
    // particle number-conserving, fermion RDMs
    std::array<std::unique_ptr<FermionRdm>, defs::c_nop_mask_frm> m_frm_rdms;

    RdmGroup(const conf::FermionRdm &opts, size_t nsite, size_t nelec){
        for (auto rank: opts.m_ranks.get())
            m_frm_rdms[rank] = smart_ptr::make_unique<FermionRdm>(opts, nsite, nelec);
    }

    /**
     * compute the 2-RDM energy
     * @return
     *  MPI-reduced sum of 0, 1, and 2 body parts of the RDM energy
     *
     *  E_RDM = h0 + h1[i,j] * rdm1[i,j] + <ij|kl> * rdm2[i,j,k,l]
     *
     *  rdm1[i,j] = sum_k rdm2[i,k,j,k] / (n_elec - 1)
     */
    defs::ham_comp_t get_energy(const FermionHamiltonian& ham) const {
        REQUIRE_TRUE_ALL(m_frm_rdms[2]!=nullptr, "cannot compute energy without the 2RDM");
        defs::ham_t e1 = 0.0;
        defs::ham_t e2 = 0.0;
        defs::wf_t trace = 0.0;
        auto row = m_frm_rdms[2]->m_store.m_row;

        for (row.restart(); row.in_range(); row.step()){
            const size_t i=row.m_inds.m_cre[0];
            const size_t j=row.m_inds.m_cre[1];
            ASSERT(i<j);
            const size_t k=row.m_inds.m_ann[0];
            const size_t l=row.m_inds.m_ann[1];
            ASSERT(k<l);
            const auto rdm_element = row.m_values[0];
            e2 += rdm_element*ham.get_element_2(i, j, k, l);
            /*
             * signs of the following contributions come from the Fermi phase of bringing the like-valued creation and
             * annihilation operators together to act on the ket, leading to a factor of (nelec - 1), accounted for
             * after the loop.
             */
            if (i == k) e1 += rdm_element*ham.m_int_1(j,l);
            if (j == l) e1 += rdm_element*ham.m_int_1(i,k);
            if (i == l) e1 -= rdm_element*ham.m_int_1(j,k);
            if (j == k) e1 -= rdm_element*ham.m_int_1(i,l);
            if ((i==k) && (j==l)) trace+=rdm_element;
        }
        // scale the one-body contribution by the number of two-body contributions
        e1 /= ham.nelec()-1;
        e1 = mpi::all_sum(e1);
        e2 = mpi::all_sum(e2);
        trace = mpi::all_sum(trace);
        ASSERT(!datatype::nearly_zero(std::abs(trace), 1e-14));
        const auto norm = datatype::real(trace) / utils::integer::combinatorial(ham.nelec(), 2);
        return datatype::real(ham.m_e_core) + (datatype::real(e1) + datatype::real(e2))/norm;
    }

};

#endif //M7_RDM_H
#endif //M7_RDM_H
