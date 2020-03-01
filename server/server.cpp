#include "yhirose/httplib.h"

#include "cfg.hpp" // cfg::GetConfig()

const auto config = cfg::GetConfig("./cfg.json");

void InitServer(void) {
    using namespace httplib;

    Server server;

    server.Get("/hi", [](const Request& req, Response& res) {
        auto text_body = "Hello World!";
        res.set_content(text_body, "text/plain");
    });

    server.Get(R"(/numbers/(\d+))", [&](const Request& req, Response& res) {
        auto numbers = req.matches[1];
        res.set_content(numbers, "text/plain");
    });

    // server.Get("/stop", [&](const Request& req, Response& res) { server.stop(); });

    const auto host = config["host"].get<std::string>().c_str();
    const auto port = config["port"].get<int>();

    server.listen(host, port);
}

int main(void) {
    InitServer();
}
