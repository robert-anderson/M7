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

    InputOptions(CLI::App &app);

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
