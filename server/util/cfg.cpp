#include "cfg.hpp"

#include <fstream>

#include "nlohmann/json.hpp"

const cfg::Config cfg::GetConfig(const std::string path) {
    std::ifstream ifs(path);
    return nlohmann::json::parse(ifs);
}
