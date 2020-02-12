//
// Created by Robert John Anderson on 2020-02-08.
//

#include "InputOptions.h"

const std::string InputOptions::description =
        "\nM7: Many-body Stochastic Expectation Value Estimation Networks\n"
        "Command line interface\n";

InputOptions::InputOptions(CLI::App &app) :m_app(app){

    add_option("-f,--fcidump-path", fcidump_path,
               "Path to the FCIDUMP file, this is read-only");

    add_option("-i,--nwalker-initial", nwalker_initial,
               "Number of walkers with which the populations will be initialized");

    add_option("-n,--nwalker-target", nwalker_target,
               "Number of walkers at which the population will be stabilised after the growth phase", true);

    add_option("-a,--nadd-initiator", nadd_initiator,
               "Number of walkers defining the initiator threshold");

    add_option("-s,--ndet-semistoch", ndet_semistoch,
               "Number of determinants selected to comprise the deterministic subspace");
}
