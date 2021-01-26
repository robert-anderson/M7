//
// Created by rja on 26/01/2021.
//

#include "Logging.h"

std::shared_ptr<spdlog::logger> g_reduced_stdout_logger = nullptr;
std::shared_ptr<spdlog::logger> g_reduced_file_logger = nullptr;
std::shared_ptr<spdlog::logger> g_local_file_logger = nullptr;