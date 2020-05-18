#include "cfg.hpp"

#include <fstream>

const cfg::Config cfg::GetConfig(std::string path) {
    std::ifstream ifs(path);
    return nlohmann::json::parse(ifs);
}
