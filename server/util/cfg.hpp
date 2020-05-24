#ifndef CFG_HPP
#define CFG_HPP

#include "nlohmann/json.hpp"

namespace cfg {
    typedef nlohmann::json Config;
    const Config GetConfig(const std::string);

    inline Config CONFIG;
}

#endif // CFG_HPP
