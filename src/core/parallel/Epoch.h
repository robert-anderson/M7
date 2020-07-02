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

#include "Distributed.h"

class Epoch {
    Distributed<size_t> m_icycle_start;

    Epoch();

    /**
     * Update the state of the Epoch, no need to check state before
     * calling update unless condition is expensive to compute.
     *
     * @param icycle the current cycle of the main loop
     * @param condition local value of the condition for the Epoch to begin
     */
    void update(size_t icycle, bool condition);

    /**
     * @return the MPI-synchronized cycle number on which the Epoch began
     */
    const size_t& start();

    /**
     * @return the MPI-synchronized state of the Epoch
     */
    bool has_begun();

};


#endif //M7_EPOCH_H
