//
// Created by Robert John Anderson on 2020-02-08.
//

#ifndef M7_INPUTOPTIONS_H
#define M7_INPUTOPTIONS_H

#include <string>
#include <external/CLI11/include/CLI/App.hpp>
#include <utility>


class InputError : public std::exception {
    const std::string m_msg;

    InputError(std::string msg) : m_msg(std::move(msg)) {}

    virtual const char *what() const throw() {
        return m_msg.c_str();
    }
};

struct InputOptions {
private:
    CLI::App &m_app;
public:
    std::string fcidump_path = "FCIDUMP";
    std::string stats_path = "M7.stats";
    bool exact_propagation = false;
    double nwalker_initial = 1.0;
    double nwalker_target = 0.0;
    double nadd_initiator = 3.0;
    double max_bloom = 1.0;
    size_t prng_seed = 0;
    size_t prng_ngen = 1000;
    size_t ndet_semistoch = 0;
    size_t spin_restrict = 0;
    double walker_factor_initial = 1.0;
    double buffer_factor_initial = 10.0;
    double min_spawn_mag = 0.0;
    size_t nload_balance_block = 10;
    double tau_initial = 0.05;
    bool dynamic_tau = false;
    size_t nenough_spawns_for_dynamic_tau = 100;
    double shift_initial = 0.0;
    double shift_damp = 1.0;
    size_t shift_update_period = 1;
    size_t ncycle = ~0ul;

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

    static const std::string description;

    bool validate() const;

};


#endif //M7_INPUTOPTIONS_H
