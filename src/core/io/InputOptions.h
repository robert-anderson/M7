//
// Created by rja on 24/05/2020.
//

#ifndef M7_INPUTOPTIONS_H
#define M7_INPUTOPTIONS_H

#include "Options.h"
#include <external/CLI11/include/CLI/App.hpp>
#include <utility>

class InputError : public std::exception {
    const std::string m_msg;

    InputError(std::string msg) : m_msg(std::move(msg)) {}

    virtual const char *what() const throw() {
        return m_msg.c_str();
    }
};

class InputOptions : public Options {
    CLI::App &m_app;
public:

    explicit InputOptions(CLI::App &app);

    template<typename T>
    void add_option(const std::string cli_options,
                    T &variable_to_bind, const std::string description, bool required = false) {
        auto opt = m_app.add_option(cli_options, variable_to_bind, description);
        if (required) opt->required();
        else opt->capture_default_str();
    }

    void add_flag(const std::string cli_options, bool &variable_to_bind, const std::string description) {
        m_app.add_flag(cli_options, variable_to_bind, description);
    }

    static const std::string program_description;


};


#endif //M7_INPUTOPTIONS_H