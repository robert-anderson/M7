//
// Created by rja on 25/08/2021.
//

#ifndef M7_EXCITITER_H
#define M7_EXCITITER_H

#include <M7_lib/basis/Suites.h>
#include <M7_lib/caches/CachedOrbsOld.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>

#include "BodyFnTypes.h"

using namespace body_fn_types;

/**
 * deterministic equivalent of the excitation generators, provides customisable loops which can be injected with
 * arbitrary code in the form of a body function argument
 */
struct ExcitIter {
    const size_t m_exsig;
    const Hamiltonian &m_ham;
protected:
    const BasisData m_bd;
    suite::Conns m_work_conn;
    suite::Mbfs m_work_dst;
    CachedOrbs m_work_orbs;
    defs::ham_t m_work_helement;
private:
    /**
     * value is set when one of the non-virtual foreach methods is called so that the canonical definition is informed
     * of whether not the helement value is required, and if so whether to discard connections on the basis that they
     * have a zero H-matrix element
     */
    bool m_need_helement = true;
    bool m_nonzero_helement_only = true;

    /**
     * clear the cached orbs, and set up the above "matrix element requirements" flags
     * @param need
     *  true if currently invoked foreach overload requires the matrix element
     * @param nonzero
     *  true if currently invoked foreach overload requests the discounting of connections with zero matrix element
     */
    void reset(bool need, bool nonzero);

protected:
    /**
     * this method should be called in canonical foreach overrides to decide whether to call the body function
     * @tparam mbf_t
     *  many-body basis function type
     * @param src
     *  source MBF
     * @param conn
     *  connection defining matrix element
     * @return
     *  false if the connection is to be discounted on the grounds that it corresponds to a zero matrix element
     */
    template<typename mbf_t>
    bool set_helement(const mbf_t &src, const conn::from_field_t<mbf_t> &conn) {
        if (!m_need_helement) return true;
        m_work_helement = m_ham.get_element(src, conn);
        return !m_nonzero_helement_only || !consts::nearly_zero(m_work_helement);
    }

private:
    /*
     * The following adaptors transform the body functions of various types for use with the canonical virtual methods.
     *
     * For now, these are kept private since the user can just call one of the foreach overloads and the adapt method
     * will be automatically called before delegating to the canonical (virtual) foreach definition. There is a case for
     * encouraging the use of the adapt method in the calling scope of foreach, since the same wrapper std::function
     * object is otherwise being constructed on every call to the foreach overload. However, this could create more
     * object lifetime problems if used incorrectly, and in any case, the creation of said wrappers should be a small
     * overhead compared to the time taken to execute the foreach method itself.
     */
    /**
     * adapts "cd" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts connection and connected MBF args
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_cd_t<mbf_t> &body_fn) {
        auto &dst_mbf = m_work_dst[mbf];
        return [&](const conn::from_field_t<mbf_t> &conn) {
            conn.apply(mbf, dst_mbf);
            body_fn(conn, dst_mbf);
        };
    }

    /**
     * adapts "ch" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts connection and H matrix element args
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_ch_t<mbf_t> &body_fn) {
        return [&](const conn::from_field_t<mbf_t> &conn) {
            DEBUG_ASSERT_EQ(m_ham.get_element(mbf, conn), m_work_helement,
                            "canonical foreach definition did not compute the correct matrix element");
            body_fn(conn, m_work_helement);
        };
    }

    /**
     * adapts "cdh" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts connection, connected MBF args, and H matrix element as args
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_cdh_t<mbf_t> &body_fn) {
        auto &dst_mbf = m_work_dst[mbf];
        return [&](const conn::from_field_t<mbf_t> &conn) {
            DEBUG_ASSERT_EQ(m_ham.get_element(mbf, conn), m_work_helement,
                            "canonical foreach definition did not compute the correct matrix element");
            conn.apply(mbf, dst_mbf);
            body_fn(conn, dst_mbf, m_work_helement);
        };
    }

    /**
     * adapts "d" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts connected MBF as its only argument
     * @param nonzero_h_only
     *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_d_t<mbf_t> &body_fn) {
        auto &dst_mbf = m_work_dst[mbf];
        return [&](const conn::from_field_t<mbf_t> &conn) {
            conn.apply(mbf, dst_mbf);
            body_fn(dst_mbf);
        };
    }

    /**
     * adapts "h" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts H matrix element as its only argument
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_h_t &body_fn) {
        return [&](const conn::from_field_t<mbf_t> &conn) {
            DEBUG_ASSERT_EQ(m_ham.get_element(mbf, conn), m_work_helement,
                            "canonical foreach definition did not compute the correct matrix element");
            body_fn(m_work_helement);
        };
    }

    /**
     * adapts "dh" closure to "c" closure
     * @tparam mbf_t
     *  many-body basis function type
     * @param mbf
     *  source many-body basis function from which to enumerate all connections and call the given function
     * @param body_fn
     *  void function which accepts connection, connected MBF args, and H matrix element as args
     * @param nonzero_h_only
     *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
     */
    template<typename mbf_t>
    fn_c_t<mbf_t> adapt(const mbf_t &mbf, const fn_dh_t<mbf_t> &body_fn) {
        auto &dst_mbf = m_work_dst[mbf];
        return [&](const conn::from_field_t<mbf_t> &conn) {
            DEBUG_ASSERT_EQ(m_ham.get_element(mbf, conn), m_work_helement,
                            "canonical foreach definition did not compute the correct matrix element");
            conn.apply(mbf, dst_mbf);
            body_fn(dst_mbf, m_work_helement);
        };
    }


public:

    ExcitIter(const Hamiltonian &ham, size_t exsig);

    virtual ~ExcitIter() {}

    /*
     * these are the "canonical" foreach definitions - the only methods which need to be overridden in derived classes
     * to change the behaviour of the loop
     */
    /**
     * @param src
     *  fermion ONV
     * @param conn
     *  workspace for setting the connection
     * @param body
     *  canonical function body object
     */
    virtual void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<FrmOnv> &body) = 0;

    /**
     * @param src
     *  fermion-boson ONV
     * @param conn
     *  workspace for setting the connection
     * @param body
     *  canonical function body object
     */
    virtual void foreach(const field::FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) = 0;

    /**
     * @param src
     *  boson ONV
     * @param conn
     *  workspace for setting the connection
     * @param body
     *  canonical function body object
     */
    virtual void foreach(const field::BosOnv &src, conn::BosOnv &conn, const fn_c_t<BosOnv> &body) = 0;


    /*
     * these are all the other foreach definitions - they do not need to be virtual because their behavior is agnostic
     * to that of the canonical definition they delegate to
     */
    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_c_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(false, nonzero_h_only);
        this->foreach(mbf, m_work_conn[mbf], body_fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_cd_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(false, nonzero_h_only);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_ch_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(true, nonzero_h_only);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_cdh_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(true, nonzero_h_only);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_d_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(false, nonzero_h_only);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_h_t &body_fn) {
        reset(true, true);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_dh_t<mbf_t> &body_fn, bool nonzero_h_only) {
        reset(true, nonzero_h_only);
        auto fn = adapt(mbf, body_fn);
        this->foreach(mbf, m_work_conn[mbf], fn);
    }
};

#endif //M7_EXCITITER_H
