#ifndef UTIL_CFG_HPP
#define UTIL_CFG_HPP

#include "nlohmann/json.hpp"

namespace cfg {
    typedef nlohmann::json Config;

    // GetConfig returns new Config that was parsed from file in JSON format
    const Config GetConfig(const std::string file);

    // kConfig is a config in JSON format, that must be initialized with
    // GetConfig() function.
    inline Config kConfig;
}

#endif // UTIL_CFG_HPP
