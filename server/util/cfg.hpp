#ifndef CFG_HPP
#define CFG_HPP

#include "nlohmann/json.hpp"

namespace cfg {
    const nlohmann::json GetConfig(std::string);
}

#endif // CFG_HPP
