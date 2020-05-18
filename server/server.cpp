#include <boost/algorithm/string.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <exception>

#include "psvt/statistics.cpp" // Yeah this source is VERY messy
#include "nlohmann/json.hpp"
#include "yhirose/httplib.h"

#include "api/isoline.hpp"                 // IsolineAPI
#include "api/stats.hpp"                   // IsolineAPI
#include "distributions/angular_gauss.hpp" // AngularGauss
#include "grids/classic_grid.hpp"          // ClassicGrid
#include "grids/spiral_grid.hpp"           // SpiralGrid
#include "util/cfg.hpp"                    // cfg::Config

namespace ublas = boost::numeric::ublas;
using json = nlohmann::json;

cfg::Config CONFIG;

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
    cov_mat(0,0) = cov[0], cov_mat(0,1) = cov[1], cov_mat(0,2) = cov[2];
    cov_mat(1,0) = cov[1], cov_mat(1,1) = cov[3], cov_mat(1,2) = cov[4];
    cov_mat(2,0) = cov[2], cov_mat(2,1) = cov[4], cov_mat(2,2) = cov[5];

    std::cout << "Will calculate grids with params: mean=" << mean.data() << ", cov=" << cov.data() << std::endl;

    AngularGauss* distribution = new AngularGauss(mean, cov_mat);
    std::cout << "Will calculate ClassicGrid: mean=" << CONFIG["classic_grid_div"].get<int>() << std::endl;
    ClassicGrid* classic_grid = new ClassicGrid(CONFIG["classic_grid_div"].get<int>(), distribution);
    std::cout << "Will calculate SpiralGrid: mean=" << CONFIG["spiral_grid_points"].get<int>() << std::endl;
    SpiralGrid* spiral_grid = new SpiralGrid(CONFIG["spiral_grid_points"].get<int>(), distribution);
    std::cout << "Done with grids" << std::endl;

    auto result = std::make_tuple(classic_grid, spiral_grid);
    database[key] = result;

    return result;
}

IsolineAPI* GetIsolineAPI(const httplib::Request& req) {
    if (!req.has_param("mean"))
        throw std::invalid_argument("no 'mean' URL param");
    if (!req.has_param("cov"))
        throw std::invalid_argument("no 'cov' URL param");
    if (!req.has_param("ratio"))
        throw std::invalid_argument("no 'ratio' URL param");

    std::vector<std::string> raw_mean_vec;
    boost::algorithm::split(raw_mean_vec, req.get_param_value("mean"), boost::is_any_of(","));
    if (raw_mean_vec.size() != 3) {
        throw std::invalid_argument("wrong 'mean' URL param");
    }

    std::array<double, 3> mean_vec;
    for (size_t i = 0; i < raw_mean_vec.size(); ++i)
        mean_vec[i] = boost::lexical_cast<double>(raw_mean_vec[i]);

    std::vector<std::string> raw_cov_vec;
    boost::algorithm::split(raw_cov_vec, req.get_param_value("cov"), boost::is_any_of(","));
    if (raw_cov_vec.size() != 6) {
        throw std::invalid_argument("wrong 'cov' URL param");
    }

    std::array<double, 6> cov_vec;
    for (size_t i = 0; i < raw_cov_vec.size(); ++i)
        cov_vec[i] = boost::lexical_cast<double>(raw_cov_vec[i]);

    std::vector<std::string> raw_ratios;
    boost::algorithm::split(raw_ratios, req.get_param_value("ratio"), boost::is_any_of(","));
    if (raw_ratios.size() == 0) {
        throw std::invalid_argument("wrong 'ratio' URL param");
    }

    std::vector<double> ratios_vec(raw_ratios.size());
    for (size_t i = 0; i < raw_ratios.size(); ++i)
        ratios_vec[i] = boost::lexical_cast<double>(raw_ratios[i]);

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;
    std::tie(classic_grid, spiral_grid) = GetGrids(mean_vec, cov_vec);

    IsolineAPI* api = new IsolineAPI(ratios_vec, classic_grid, spiral_grid); // TODO clever pointer + destructors (for everything lol)
    return api;
}

StatsAPI* GetStatsAPI(const httplib::Request& req) {
    json params = json::parse(req.body);

    if (params.find("mean") == params.end())
        throw std::invalid_argument("no 'mean' data");
    if (params.find("cov") == params.end())
        throw std::invalid_argument("no 'cov' data");
    if (params.find("points") == params.end())
        throw std::invalid_argument("no 'points' data");

    std::vector<std::string> raw_mean_vec;
    boost::algorithm::split(raw_mean_vec, params["mean"].get<std::string>(), boost::is_any_of(","));
    if (raw_mean_vec.size() != 3) {
        throw std::invalid_argument("wrong 'mean' data");
    }

    std::array<double, 3> mean_vec;
    for (size_t i = 0; i < raw_mean_vec.size(); ++i)
        mean_vec[i] = boost::lexical_cast<double>(raw_mean_vec[i]);

    std::vector<std::string> raw_cov_vec;
    boost::algorithm::split(raw_cov_vec, params["cov"].get<std::string>(), boost::is_any_of(","));
    if (raw_cov_vec.size() != 6) {
        throw std::invalid_argument("wrong 'cov' data");
    }

    std::array<double, 6> cov_vec;
    for (size_t i = 0; i < raw_cov_vec.size(); ++i)
        cov_vec[i] = boost::lexical_cast<double>(raw_cov_vec[i]);

    std::vector<std::string> raw_points;
    boost::algorithm::split(raw_points, params["points"].get<std::string>(), boost::is_any_of(","));
    if (raw_points.size() == 0 || raw_points.size()%3 != 0) {
        throw std::invalid_argument("wrong 'points' data");
    }

    std::vector<std::array<double, 3>> points_vec(raw_points.size() / 3);
    for (size_t i = 0; i < raw_points.size(); i += 3) {
        std::array<double, 3> point = {
            boost::lexical_cast<double>(raw_points[i]),
            boost::lexical_cast<double>(raw_points[i+1]),
            boost::lexical_cast<double>(raw_points[i+2])
        };
        points_vec[i/3] = point;
    }

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;
    std::tie(classic_grid, spiral_grid) = GetGrids(mean_vec, cov_vec);

    StatsAPI* api = new StatsAPI(points_vec, classic_grid, spiral_grid); // TODO clever pointer + destructors (for everything lol)
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
            auto api = GetIsolineAPI(req);
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

    server.Post("/api/stats", [&](const Request& req, Response& res) {
        res.status = 200;
        res.set_header("Access-Control-Allow-Origin", "*"); // TODO
        res.set_header("Access-Control-Request-Method", "POST");
        res.set_header("Access-Control-Allow-Headers", "*");
        json response;

        try {
            auto api = GetStatsAPI(req);
            auto points = api->points;

            std::vector<double> s_stats;
            s_stats.reserve(points.size());
            int s_skipped_counter = 0;

            std::vector<double> t_stats;
            t_stats.reserve(points.size());
            int t_skipped_counter = 0;

            for (int i = 0; i < points.size(); ++i) {
                double value = api->spiral_grid->distr->Calc(points[i]);
                try {
                    t_stats.push_back(api->CalculateTStat(value));
                } catch(std::exception& e) {
                    std::cerr << "error: could not calculate t stat for point #" << i << " (" <<
                    points[i][0] << ", " << points[i][1] << ", " << points[i][2] << "): " << e.what() << std::endl;
                    std::cerr << "warning: t stat for point #" << i << " will be skipped" << std::endl;
                    ++t_skipped_counter;
                }
                try {
                    s_stats.push_back(api->CalculateSStat(points[i], value));
                } catch(std::exception& e) {
                    std::cerr << "error: could not calculate s stat for point #" << i << " (" <<
                    points[i][0] << ", " << points[i][1] << ", " << points[i][2] << "): " << e.what() << std::endl;
                    std::cerr << "warning: s stat for point #" << i << " will be skipped" << std::endl;
                    ++s_skipped_counter;
                }
            }

            if (t_skipped_counter)
                std::cerr << "warning: " << t_skipped_counter << " points were skipped for t stat calculation" << std::endl;
            if (s_skipped_counter)
                std::cerr << "warning: " << s_skipped_counter << " points were skipped for s stat calculation" << std::endl;

            if (t_stats.size() != 0) {
                double *OrdUnif = new double[t_stats.size()];
                double KS_measure, KS_estim, AD_measure, AD_estim;
                StTests1D(OrdUnif, KS_measure, KS_estim, AD_measure, AD_estim, t_stats.data(), t_stats.size());
                std::cout << "t statistics results:" << std::endl;
                std::cout << "KS_measure=" << KS_measure << " KS_estim=" << KS_estim << std::endl;
                std::cout << "AD_measure=" << AD_measure << " AD_estim=" << AD_estim << std::endl;
                delete[] OrdUnif;
            }

            if (s_stats.size() != 0) {
                double *OrdUnif = new double[s_stats.size()];
                double KS_measure, KS_estim, AD_measure, AD_estim;
                StTests1D(OrdUnif, KS_measure, KS_estim, AD_measure, AD_estim, s_stats.data(), s_stats.size());
                std::cout << "s statistics results:" << std::endl;
                std::cout << "KS_measure=" << KS_measure << " KS_estim=" << KS_estim << std::endl;
                std::cout << "AD_measure=" << AD_measure << " AD_estim=" << AD_estim << std::endl;
                delete[] OrdUnif;
            }

            response["body"]["t"] = t_stats;
            response["body"]["s"] = s_stats;

        } catch (std::exception& e) {
            response["code"] = 500;
            response["body"] = e.what();
            res.set_content(response.dump(), "application/json");
            return;
        }

        response["code"] = 200;
        res.set_content(response.dump(), "application/json");
    });

    server.Options("/api/stats", [](const Request& req, Response& res) {
        res.status = 204;
        res.set_header("Access-Control-Allow-Origin", "*");
        res.set_header("Access-Control-Request-Method", "POST");
        res.set_header("Access-Control-Allow-Headers", "*");
    });

    server.Get("/stop", [&](const Request& req, Response& res) { server.stop(); });

    const auto host = CONFIG["host"].get<std::string>();
    const auto port = CONFIG["port"].get<int>();

    server.listen(host.c_str(), port);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "please specify config path as second argument" << std::endl;
        exit(1);
    }

    CONFIG = cfg::GetConfig(argv[1]);
    InitServer();
}
