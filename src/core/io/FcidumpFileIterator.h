//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_FCIDUMPFILEITERATOR_H
#define M7_FCIDUMPFILEITERATOR_H

#include "TensorFileIterator.h"

static const std::regex header_terminator_regex{R"(\&END)"};

static constexpr std::array<std::array<size_t, 4>, 8> orderings{
        {
                {0, 1, 2, 3},
                {1, 0, 2, 3},
                {0, 1, 3, 2},
                {1, 0, 3, 2},
                {2, 3, 0, 1},
                {2, 3, 1, 0},
                {3, 2, 0, 1},
                {3, 2, 1, 0}
        }
};

template<typename T>
class FcidumpFileIterator : public TensorFileIterator<T> {
public:
    const size_t m_norb;
    const size_t m_isymm;
    const size_t m_nelec;
    const defs::inds m_orbsym;
    const bool m_spin_resolved;

    FcidumpFileIterator(const std::string &filename) :
            TensorFileIterator<T>(filename, 4, false),
            m_norb(read_header_int(filename, "NORB")), m_isymm(m_norb>3?isymm(filename):1),
            m_nelec(read_header_int(filename, "NELEC")),
            m_orbsym(read_header_array(filename, "ORBSYM")),
            m_spin_resolved(read_header_bool(filename, "UHF") || read_header_bool(filename, "TREL")) {
        /*
         * a valid FCIDUMP file will correspond to a tensor with all extents in the shape
         * array equal to one another
         */
        ASSERT(std::adjacent_find(
                TensorFileIterator<T>::m_shape.begin(),
                TensorFileIterator<T>::m_shape.end(),
                std::not_equal_to<size_t>()) == TensorFileIterator<T>::m_shape.end())
        /*
         * the number of orbitals read from the header should be equal to the elements
         * of the tensor shape array
         */
        ASSERT(m_norb == TensorFileIterator<T>::m_shape[0]);
    }

    size_t nspinorb() const{
        return m_spin_resolved?m_norb:m_norb*2;
    }

    size_t nsite() const{
        return m_spin_resolved?m_norb/2:m_norb;
    }
private:

    /*
     *  Chemical integral symmetries
     *
     *  [ij|kl] = i*(1)j(1) 1/r k*(2)l(2)
     *
     *  if shape >= 4
     *  [01|23] [10|23] [01|32] [10|32]
     *  [23|01] [23|10] [32|01] [32|10]
     *
     *  if shape = 3
     *  [01|22] [10|22]
     *  [22|01] [22|10]
     *
     *  if shape = 2
     *  [00|11]
     *  [11|00]
     *
     *  if shape = 1
     *  [00|00]
     *
     */

    static size_t isymm(const std::string &filename) {
        auto index_is_defined = [](size_t i) { return i < (size_t) (-1); };
        TensorFileIterator<T> tensorFileIterator(filename, 4, false);
        defs::inds inds(4);
        T value;
        // this will eventually hold all orderings of the first example of an
        // integral with 4 distinct indices
        std::array<defs::inds, 8> inds_distinct{};
        size_t isymm{};
        while (tensorFileIterator.next(inds, value)) {
            if (std::all_of(inds.begin(), inds.end(), index_is_defined)) {
                // we have a two body integral
                if (!isymm) {
                    inds_distinct[0].assign(inds.begin(), inds.end());
                    std::sort(inds.begin(), inds.end());
                    // still looking for an example of four distinct indices
                    if (std::adjacent_find(inds.begin(), inds.end()) == inds.end()) {
                        for (size_t i = 1ul; i < orderings.size(); ++i) {
                            for (size_t j = 0ul; j < 4; ++j)
                                inds_distinct[i].push_back(inds_distinct[0][orderings[i][j]]);
                        }
                        isymm++;
                    }
                } else {
                    for (auto tmp : inds_distinct) {
                        if (std::equal(inds.begin(), inds.end(), tmp.begin())){
                            isymm++;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT(isymm);
        // 8->1; 4->2; 2->4; 1->8
        return 8 / isymm;
    }

    static size_t read_header_int(const std::string &fname, const std::string &label, size_t default_ = 0) {
        const auto regex = std::regex(label + R"(\s*\=\s*[0-9]+)");
        FileIterator iterator(fname);
        std::string line;
        std::smatch match;
        std::string match_string;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) {
                std::string match_string(match.str());
                std::regex_search(match_string, match, std::regex(R"([0-9]+)"));
                return std::atol(match.str().c_str());
            }
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    static defs::inds read_header_array(const std::string &fname, const std::string &label) {
        return defs::inds{0, 1};
    }

    static size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_ = false) {
        std::regex regex;
        if (default_) regex = std::regex(label + R"(\s?\=\s?\.FALSE\.)");
        else regex = std::regex(label + R"(\s?\=\s?\.TRUE\.)");
        FileIterator iterator(fname);
        std::string line;
        std::smatch match;
        std::string match_string;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) return !default_;
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }
};

#endif //M7_FCIDUMPFILEITERATOR_H
