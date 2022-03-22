//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONVCONNECTION_H
#define M7_BOSONVCONNECTION_H

#include <M7_lib/field/BosOnvField.h>
#include <M7_lib/parallel/MPIAssert.h>

struct BosOpPair {
    const size_t m_imode;
    const size_t m_nop;
    BosOpPair(size_t imode, size_t nop);
};

class BosOps {
    std::vector<BosOpPair> m_pairs;
    size_t m_nop = 0ul;
public:
    BosOps(size_t nmode);

    const std::vector<BosOpPair>& pairs() const;

    defs::inds to_vector() const;

    void from_vector(const defs::inds& imodes);

    void clear();

    size_t size() const;

    void add(BosOpPair&& pair);

    void set(BosOpPair&& pair);

    void set(const size_t& imode);

    void set(const size_t& imode, const size_t& jmode);

    const BosOpPair& operator[](const size_t& ipair) const;

    /**
     * @param iop
     *  operator index in the expanded product
     * @return
     *  mode index corresponding to the operator index
     */
    size_t get_imode(size_t iop) const;
};

struct BosOnvConnection {
    BosOps m_ann, m_cre;

    explicit BosOnvConnection(size_t nmode);

    explicit BosOnvConnection(BasisData bd);

    explicit BosOnvConnection(const BosOnvField& mbf);

    void clear();

    size_t size() const;

    void connect(const BosOnvField& src, const BosOnvField& dst);

    void connect(const BosOnvField& src, const BosOnvField& dst, BosOps& com);

    void apply(const BosOnvField &src, BosOnvField &dst) const;

    void apply(const BosOnvField &src, BosOps &com) const;

    void apply(const BosOnvField &src, BosOnvField &dst, BosOps &com) const;

    size_t exsig() const;

    bool respects_occ_range(const BosOnvField &src, size_t nboson_max) const;

    /**
     * compute the "occupation factor" required to keep the boson ONV basis orthonormal.
     * b |n>  = sqrt(n) |n-1>
     * b+ |n> = sqrt(n+1) |n+1>
     *
     * b^m |n>  = sqrt(n(n-1)...(n-m+1)) |n-m>
     * b+^m |n>  = sqrt((n+1)(n+2)...(n+m)) |n+m>
     * @param src
     * @return
     */
    size_t occ_fac_square(const BosOnvField& src) const;

    double occ_fac(const BosOnvField& src) const;

private:
    static size_t occ_fac_square_com(const size_t& occ, const size_t& nop_com);

public:

    /**
     * if the matrix element in question is actually contracted, with common indices among the creation and annihilation
     * operators, the occupation factors change.
     * the ONV can be rearranged to group operators in mode-wise normal order. for a given mode, common indices result
     * in the following factor:
     * b+b |n> = n |n>
     * b+^2b^2 |n> = n(n-1) |n>
     * b+^mb^m |n> = n(n-1)...(n-m+1) |n>
     * @param src
     * @param com
     * @return
     */
    size_t occ_fac_square(const BosOnvField& src, const BosOps& com) const;

    double occ_fac(const BosOnvField& src, const BosOps& com) const;

};


#endif //M7_BOSONVCONNECTION_H