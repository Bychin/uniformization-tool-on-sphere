#include "cfg.hpp"

#include <fstream>

const nlohmann::json cfg::GetConfig(std::string path) {
    std::ifstream ifs(path);
    return nlohmann::json::parse(ifs);
}
