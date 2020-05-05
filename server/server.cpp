#include <boost/algorithm/string.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include "api/isoline.hpp"                 // IsolineAPI
#include "distributions/angular_gauss.hpp" // AngularGauss
#include "grids/classic_grid.hpp" // ClassicGrid
#include "grids/spiral_grid.hpp" // SpiralGrid
#include "nlohmann/json.hpp"
#include "util/cfg.hpp" // cfg::GetConfig()
#include "yhirose/httplib.h"

namespace ublas = boost::numeric::ublas;
using json = nlohmann::json;

const auto config = cfg::GetConfig("./cfg.json"); // TODO command line argument

// TODO std::tuple<std::array<double, 3>, std::array<double, 6>> -> AngularGaussParams
// database stores distributions with raw mean and cov parameters.
std::unordered_map<std::tuple<std::array<double, 3>, std::array<double, 6>>,
                   std::tuple<ClassicGrid*, SpiralGrid*>,
                   boost::hash<std::tuple<std::array<double, 3>, std::array<double, 6>>>>
    database;

std::tuple<ClassicGrid*, SpiralGrid*> GetGrids(std::array<double, 3>& mean, std::array<double, 6>& cov) {
    auto key = std::make_tuple(mean, cov);
    if (database.find(key) != database.end())
        return database[key];

    ublas::matrix<double> cov_mat(3, 3); // TODO better appr and func
    cov_mat(0,0) = cov[0];
    cov_mat(0,1) = cov[1];
    cov_mat(0,2) = cov[2];
    cov_mat(1,0) = cov[1];
    cov_mat(1,1) = cov[3];
    cov_mat(1,2) = cov[4];
    cov_mat(2,0) = cov[2];
    cov_mat(2,1) = cov[4];
    cov_mat(2,2) = cov[5];

    std::cout << "Will calculate grids with params: mean=" << mean.data() << ", cov=" << cov.data() << std::endl;

    AngularGauss* distribution = new AngularGauss(mean, cov_mat);
    std::cout << "Will calculate ClassicGrid: mean=" << config["classic_grid_div"].get<int>() << std::endl;
    ClassicGrid* classic_grid = new ClassicGrid(config["classic_grid_div"].get<int>(), distribution);
    std::cout << "Will calculate SpiralGrid: mean=" << config["spiral_grid_points"].get<int>() << std::endl;
    SpiralGrid* spiral_grid = new SpiralGrid(config["spiral_grid_points"].get<int>(), distribution);
    std::cout << "Done with grids" << std::endl;

    auto result = std::make_tuple(classic_grid, spiral_grid);
    database[key] = result;

    return result;
}

IsolineAPI* GetIsolineAPI(const httplib::Request& req) {
    if (!req.has_param("mean"))
        throw "no 'mean' URL param";
    if (!req.has_param("cov"))
        throw "no 'cov' URL param";
    if (!req.has_param("ratio"))
        throw "no 'ratio' URL param";

    std::vector<std::string> raw_mean_vec;
    boost::algorithm::split(raw_mean_vec, req.get_param_value("mean"), boost::is_any_of(","));
    if (raw_mean_vec.size() != 3) {
        throw "wrong 'mean' URL param";
    }

    std::array<double, 3> mean_vec;
    for (size_t i = 0; i < raw_mean_vec.size(); ++i)
        mean_vec[i] = boost::lexical_cast<double>(raw_mean_vec[i]);

    std::vector<std::string> raw_cov_vec;
    boost::algorithm::split(raw_cov_vec, req.get_param_value("cov"), boost::is_any_of(","));
    if (raw_cov_vec.size() != 6) {
        throw "wrong 'cov' URL param";
    }

    std::array<double, 6> cov_vec;
    for (size_t i = 0; i < raw_cov_vec.size(); ++i)
        cov_vec[i] = boost::lexical_cast<double>(raw_cov_vec[i]);

    std::vector<std::string> raw_ratios;
    boost::algorithm::split(raw_ratios, req.get_param_value("ratio"), boost::is_any_of(","));
    if (raw_ratios.size() == 0) {
        throw "wrong 'ratio' URL param";
    }

    std::vector<double> ratios_vec(raw_ratios.size());
    for (size_t i = 0; i < raw_ratios.size(); ++i)
        ratios_vec[i] = boost::lexical_cast<double>(raw_ratios[i]);

    for (auto& i: ratios_vec)
        std::cout << i << ' ';

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;
    std::tie(classic_grid, spiral_grid) = GetGrids(mean_vec, cov_vec);

    IsolineAPI* api = new IsolineAPI(mean_vec, cov_vec, ratios_vec, classic_grid, spiral_grid); // TODO clever pointer + destructors (for everything lol)
    return api;
}

void InitServer(void) {
    using namespace httplib;

    Server server;

    server.Get("/api/isoline", [&](const Request& req, Response& res) {
        res.status = 200;
        res.set_header("Access-Control-Allow-Origin", "*"); // TODO
        json response;

        try {
            auto api =GetIsolineAPI(req);
            auto isolines = api->GetIsolines();
            response["body"] = isolines;
        } catch (std::exception& e) {
            response["code"] = 500;
            response["body"] = e.what();
            res.set_content(response.dump(), "application/json");
            return;
        }

        response["code"] = 200;
        res.set_content(response.dump(), "application/json");
    });

    server.Get("/stop", [&](const Request& req, Response& res) { server.stop(); });

    const auto host = config["host"].get<std::string>();
    const auto port = config["port"].get<int>();

    server.listen(host.c_str(), port);
}

int main(void) {
    InitServer();
}
