/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 * The behaviour of the algorithms implemented in this code, or the on/off
 * status of certain subalgorithms commonly changes when some condition is
 * fulfilled.
 * e.g. in an FCIQMC calculation, a crucial idea is the transition into
 * variable shift behaviour, the intiation of the semistochastic adaptation,
 * or the onset of RDM accumulation. These three concepts can be represented
 * by an instance of the Epoch class.
 *
 * When an FciqmcCalculation enters an Epoch, all MPI ranks must be informed,
 * and the class should expose methods to report the MPI-synchronized:
 * 1. boolean state of the Epoch
 * 2. cycle number on which the Epoch began
 */

#ifndef M7_EPOCH_H
#define M7_EPOCH_H

#include <utility>
#include <string>
#include <vector>

#include "Reduction.h"
#include "MPIAssert.h"

class Epoch {
    str_t m_name;
    Reduction<uint_t> m_icycle_start;

public:
    explicit Epoch(str_t name);

    /**
     * Update the state of the Epoch, no need to check state before
     * calling update unless condition is expensive to compute.
     *
     * @param icycle the current cycle of the main loop
     * @param condition local value of the condition for the Epoch to begin
     * @return true if Epoch begins on icycle, else false
     */
    bool update(uint_t icycle, bool condition);

    void terminate(uint_t icycle);

    /**
     * @return the MPI-synchronized cycle number on which the Epoch began
     */
    const uint_t& icycle_start() const;

    /**
     * @return the MPI-synchronized state of the Epoch
     */
    operator bool() const;

    bool started_last_cycle(uint_t icycle) const;

    bool started_this_cycle(uint_t icycle) const;

};

class Epochs {
    std::vector<Epoch> m_epochs;
public:
    Epochs(str_t name, uint_t n, str_t element_identifier="element"){
        m_epochs.reserve(n);
        for (uint_t i=0ul; i<n; ++i) m_epochs.emplace_back(name+" (" +element_identifier + " " + std::to_string(i) + ")");
    }

    uint_t nelement() const {
        return m_epochs.size();
    }

    Epoch& operator[](const uint_t i){
        DEBUG_ASSERT_LT(i, m_epochs.size(), "Epoch index OOB");
        return m_epochs[i];
    }

    const Epoch& operator[](const uint_t i) const{
        DEBUG_ASSERT_LT(i, m_epochs.size(), "Epoch index OOB");
        return m_epochs[i];
    }

    void terminate(const uint_t& icycle){
        for (auto& epoch:m_epochs) epoch.terminate(icycle);
    }

    /**
     * @return
     * the cycle on which the last epoch to start started, or ~0ul if any epochs have yet to start
     */
    uint_t icycle_start_last() const;

    /**
     * @return
     *  true only if all epochs have begun
     */
    operator bool() const;

};


#endif //M7_EPOCH_H
