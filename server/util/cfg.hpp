#ifndef CFG_HPP
#define CFG_HPP

#include "nlohmann/json.hpp"

namespace cfg {
    typedef nlohmann::json Config;
    const Config GetConfig(std::string);
}

#endif // CFG_HPP
