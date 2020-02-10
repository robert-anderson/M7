//
// Created by Robert John Anderson on 2020-02-08.
//

#ifndef M7_INPUTOPTIONS_H
#define M7_INPUTOPTIONS_H

#include <string>
#include <external/CLI11/include/CLI/App.hpp>

class InputOptions {

    CLI::App &m_app;

public:
    std::string fcidump_path = "FCIDUMP";
    double nwalker_initial = 1.0;
    double nwalker_target = 0.0;
    double nadd_initiator = 3.0;
    size_t ndet_semistoch = 0;

    InputOptions(CLI::App &app):m_app(app){

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

    template <typename T>
    void add_option(const std::string cli_options,
            T &variable_to_bind, const std::string description, bool required=false){
        auto opt = m_app.add_option(cli_options, variable_to_bind, description);
        if (required) opt->required();
        else opt->capture_default_str();
    }

    static const std::string description;

};


#endif //M7_INPUTOPTIONS_H
