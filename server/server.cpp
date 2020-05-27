#include <chrono>
#include <exception>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "nlohmann/json.hpp"
#include "psvt/statistics.hpp"
#include "yhirose/httplib.h"

#include "api/isoline.hpp"                 // IsolineAPI
#include "api/stats.hpp"                   // StatsAPI
#include "distributions/angular_gauss.hpp" // AngularGauss
#include "grids/classic_grid.hpp"          // ClassicGrid
#include "grids/spiral_grid.hpp"           // SpiralGrid
#include "util/cfg.hpp"                    // cfg::kConfig

namespace ublas = boost::numeric::ublas;
using json = nlohmann::json;

// GridsCache stores grids by raw mean and cov parameters.
typedef std::unordered_map<std::tuple<std::array<double, 3>, std::array<double, 6>>,
                           std::tuple<ClassicGrid*, SpiralGrid*>,
                           boost::hash<std::tuple<std::array<double, 3>, std::array<double, 6>>>>
    GridsCache;

GridsCache* cache;

std::tuple<ClassicGrid*, SpiralGrid*> GetGrids(std::array<double, 3>& mean, std::array<double, 6>& cov) {
    auto key = std::make_tuple(mean, cov);
    if (cache->find(key) != cache->end()) {
        std::cout << "debug: got grids from cache" << std::endl;
        return cache->at(key);
    }

    ublas::matrix<double> cov_mat(3, 3);
    cov_mat(0,0) = cov[0], cov_mat(0,1) = cov[1], cov_mat(0,2) = cov[2];
    cov_mat(1,0) = cov[1], cov_mat(1,1) = cov[3], cov_mat(1,2) = cov[4];
    cov_mat(2,0) = cov[2], cov_mat(2,1) = cov[4], cov_mat(2,2) = cov[5];

    std::cout << "debug: going to calculate grids with params: mean=[" << mean[0] << "," << mean[1] << "," << mean[2] << "], cov=["
        << cov[0] << "," << cov[1] << "," << cov[2] << "," << cov[3] << "," << cov[4] << "," << cov[5] << "]\n";

    AngularGauss* distribution = new AngularGauss(mean, cov_mat);

    std::cout << "debug: going to calculate ClassicGrid with bw=" << cfg::kConfig["classic_grid_div"].get<int>() << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    ClassicGrid* classic_grid = new ClassicGrid(cfg::kConfig["classic_grid_div"].get<int>(), distribution);
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "debug: ClassicGrid is ready in " << duration.count() << "ms" << std::endl;

    std::cout << "debug: going to calculate SpiralGrid with n=" << cfg::kConfig["spiral_grid_points"].get<int>() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    SpiralGrid* spiral_grid = new SpiralGrid(cfg::kConfig["spiral_grid_points"].get<int>(), distribution);
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "debug: SpiralGrid is ready in " << duration.count() << "ms" << std::endl;

    auto result = std::make_tuple(classic_grid, spiral_grid);
    (*cache)[key] = result;

    return result;
}

IsolineAPI* GetIsolineAPI(const httplib::Request& req) {
    if (!req.has_param("mean"))
        throw std::invalid_argument("invalid request: no 'mean' URL param");
    if (!req.has_param("cov"))
        throw std::invalid_argument("invalid request: no 'cov' URL param");
    if (!req.has_param("ratio"))
        throw std::invalid_argument("invalid request: no 'ratio' URL param");

    std::vector<std::string> raw_mean_vec;
    boost::algorithm::split(raw_mean_vec, req.get_param_value("mean"), boost::is_any_of(","));
    if (raw_mean_vec.size() != 3) {
        throw std::invalid_argument("invalid request: wrong 'mean' URL param");
    }

    std::array<double, 3> mean_vec;
    for (size_t i = 0; i < raw_mean_vec.size(); ++i)
        mean_vec[i] = boost::lexical_cast<double>(raw_mean_vec[i]);

    std::vector<std::string> raw_cov_vec;
    boost::algorithm::split(raw_cov_vec, req.get_param_value("cov"), boost::is_any_of(","));
    if (raw_cov_vec.size() != 6) {
        throw std::invalid_argument("invalid request: wrong 'cov' URL param");
    }

    std::array<double, 6> cov_vec;
    for (size_t i = 0; i < raw_cov_vec.size(); ++i)
        cov_vec[i] = boost::lexical_cast<double>(raw_cov_vec[i]);

    std::vector<std::string> raw_ratios;
    boost::algorithm::split(raw_ratios, req.get_param_value("ratio"), boost::is_any_of(","));
    if (raw_ratios.size() == 0) {
        throw std::invalid_argument("invalid request: wrong 'ratio' URL param");
    }

    std::vector<double> ratios_vec(raw_ratios.size());
    for (size_t i = 0; i < raw_ratios.size(); ++i)
        ratios_vec[i] = boost::lexical_cast<double>(raw_ratios[i]);

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;
    std::tie(classic_grid, spiral_grid) = GetGrids(mean_vec, cov_vec);

    IsolineAPI* api = new IsolineAPI(ratios_vec, classic_grid, spiral_grid);
    return api;
}

StatsAPI* GetStatsAPI(const httplib::Request& req) {
    json params = json::parse(req.body);

    if (params.find("mean") == params.end())
        throw std::invalid_argument("invalid request: no 'mean' data");
    if (params.find("cov") == params.end())
        throw std::invalid_argument("invalid request: no 'cov' data");
    if (params.find("points") == params.end())
        throw std::invalid_argument("invalid request: no 'points' data");

    std::vector<std::string> raw_mean_vec;
    boost::algorithm::split(raw_mean_vec, params["mean"].get<std::string>(), boost::is_any_of(","));
    if (raw_mean_vec.size() != 3) {
        throw std::invalid_argument("invalid request: wrong 'mean' data");
    }

    std::array<double, 3> mean_vec;
    for (size_t i = 0; i < raw_mean_vec.size(); ++i)
        mean_vec[i] = boost::lexical_cast<double>(raw_mean_vec[i]);

    std::vector<std::string> raw_cov_vec;
    boost::algorithm::split(raw_cov_vec, params["cov"].get<std::string>(), boost::is_any_of(","));
    if (raw_cov_vec.size() != 6) {
        throw std::invalid_argument("invalid request: wrong 'cov' data");
    }

    std::array<double, 6> cov_vec;
    for (size_t i = 0; i < raw_cov_vec.size(); ++i)
        cov_vec[i] = boost::lexical_cast<double>(raw_cov_vec[i]);

    std::vector<std::string> raw_points;
    boost::algorithm::split(raw_points, params["points"].get<std::string>(), boost::is_any_of(","));
    if (raw_points.size() == 0 || raw_points.size()%3 != 0) {
        throw std::invalid_argument("invalid request: wrong 'points' data");
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

    StatsAPI* api = new StatsAPI(points_vec, classic_grid, spiral_grid);
    return api;
}

void InitServer(void) {
    using namespace httplib;

    Server server;

    server.Get("/api/isoline", [&](const Request& req, Response& res) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "debug: /api/isoline: got new GET request\n";

        res.status = 200;
        res.set_header("Access-Control-Allow-Origin", "*");
        json response;

        try {
            auto api = GetIsolineAPI(req);
            api->Validate();
            auto isolines = api->GetIsolines();
            response["body"] = isolines;
            delete api;
        } catch (std::exception& e) {
            response["code"] = 500;
            response["body"] = e.what();
            res.set_content(response.dump(), "application/json");
            return;
        }

        response["code"] = 200;
        res.set_content(response.dump(), "application/json");

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << "debug: /api/isoline: response is ready in " << duration.count() << "ms" << std::endl;
    });

    server.Post("/api/stats", [&](const Request& req, Response& res) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "debug: /api/stats: got new POST request\n";

        res.status = 200;
        res.set_header("Access-Control-Allow-Origin", "*");
        res.set_header("Access-Control-Request-Method", "POST");
        res.set_header("Access-Control-Allow-Headers", "*");
        json response;

        try {
            auto api = GetStatsAPI(req);
            api->Validate();
            auto points = api->Points();

            std::vector<double> s_stats;
            s_stats.reserve(points.size());
            int s_skipped_counter = 0;

            std::vector<double> t_stats;
            t_stats.reserve(points.size());
            int t_skipped_counter = 0;

            bool s_stat_with_debug_info = cfg::kConfig["s_stat_with_debug_info"].get<bool>();
            std::vector<json> s_stats_debug_info;
            s_stats_debug_info.reserve(points.size());

            for (int i = 0; i < points.size(); ++i) {
                double value = api->CalculateFunc(points[i]);
                try {
                    t_stats.push_back(api->CalculateTStat(value));
                } catch(std::exception& e) {
                    std::cerr << "error: could not calculate t stat for point #" << i << " (" <<
                        points[i][0] << ", " << points[i][1] << ", " << points[i][2] << "): " << e.what() << std::endl;
                    std::cout << "warning: t stat for point #" << i << " will be skipped" << std::endl;
                    ++t_skipped_counter;
                }
                try {
                    if (!s_stat_with_debug_info) {
                        s_stats.push_back(api->CalculateSStat(points[i], value));
                        continue;
                    }

                    auto result = api->CalculateSStatWithDebugInfo(points[i], value);
                    s_stats.push_back(result["s"].get<double>());
                    s_stats_debug_info.push_back(result);

                } catch(std::exception& e) {
                    std::cerr << "error: could not calculate s stat for point #" << i << " (" <<
                        points[i][0] << ", " << points[i][1] << ", " << points[i][2] << "): " << e.what() << std::endl;
                    std::cout << "warning: s stat for point #" << i << " will be skipped" << std::endl;
                    ++s_skipped_counter;
                }
            }

            if (t_skipped_counter)
                std::cout << "warning: " << t_skipped_counter << " points were skipped for t stat calculation" << std::endl;
            if (s_skipped_counter)
                std::cout << "warning: " << s_skipped_counter << " points were skipped for s stat calculation" << std::endl;

            if (t_stats.size() != 0) {
                auto tests_result = UniformityTests(t_stats);
                response["body"]["t_results"] = json({});
                response["body"]["t_results"]["ks_measure"] = tests_result.KS_measure, response["body"]["t_results"]["ks_est"] = tests_result.KS_estimate;
                response["body"]["t_results"]["ad_measure"] = tests_result.AD_measure, response["body"]["t_results"]["ad_est"] = tests_result.AD_estimate;
                std::cout << "info: t statistics results:" << std::endl;
                std::cout << "info: KS_measure=" << tests_result.KS_measure << " KS_estimate=" << tests_result.KS_estimate << std::endl;
                std::cout << "info: AD_measure=" << tests_result.AD_measure << " AD_estimate=" << tests_result.AD_estimate << std::endl;
            }

            if (s_stats.size() != 0) {
                auto tests_result = UniformityTests(s_stats);
                response["body"]["s_results"] = json({});
                response["body"]["s_results"]["ks_measure"] = tests_result.KS_measure, response["body"]["s_results"]["ks_est"] = tests_result.KS_estimate;
                response["body"]["s_results"]["ad_measure"] = tests_result.AD_measure, response["body"]["s_results"]["ad_est"] = tests_result.AD_estimate;
                std::cout << "info: s statistics results:" << std::endl;
                std::cout << "info: KS_measure=" << tests_result.KS_measure << " KS_estimate=" << tests_result.KS_estimate << std::endl;
                std::cout << "info: AD_measure=" << tests_result.AD_measure << " AD_estimate=" << tests_result.AD_estimate << std::endl;
            }

            response["body"]["t"] = t_stats;
            response["body"]["s"] = s_stats;

            if (s_stat_with_debug_info)
                response["body"]["s_stats_debug_info"] = s_stats_debug_info;

            delete api;

        } catch (std::exception& e) {
            response["code"] = 500;
            response["body"] = e.what();
            res.set_content(response.dump(), "application/json");
            return;
        }

        response["code"] = 200;
        res.set_content(response.dump(), "application/json");

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << "debug: /api/stats: response is ready in " << duration.count() << "ms" << std::endl;
    });

    server.Options("/api/stats", [](const Request& req, Response& res) {
        res.status = 204;
        // WARNING: for local use only
        res.set_header("Access-Control-Allow-Origin", "*");
        res.set_header("Access-Control-Request-Method", "POST");
        res.set_header("Access-Control-Allow-Headers", "*");
    });

    server.Get("/stop", [&](const Request& req, Response& res) {server.stop();});

    const auto host = cfg::kConfig["host"].get<std::string>();
    const auto port = cfg::kConfig["port"].get<int>();

    server.listen(host.c_str(), port);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "please specify config path as second argument" << std::endl;
        exit(1);
    }

    cache = new GridsCache;
    cfg::kConfig = cfg::GetConfig(argv[1]);
    InitServer();
}
