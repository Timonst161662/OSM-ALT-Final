#include <crow/app.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <thread>
#include <chrono>

#include "Dijkstra.h"
#include "GraphUtils.h"
#include "AStar.h"
#include "ALT.h"

using namespace std;

// Globale Strukturen
unordered_map<int, Node> nodeMap;
vector<vector<pair<int, double>>> adjListVec;
WeightedAdjList adjListW;

// Landmark-Daten
std::vector<int> landmarkIds;
std::vector<std::vector<double>> tables_forward;
std::vector<std::vector<double>> tables_backward;

int main_all5Algos(int argc, char* argv[]) {
    std::cout << "[DEBUG] CWD: " << std::filesystem::current_path() << std::endl;

    // Graph laden
    auto [parsedNodeMap, parsedAdjListVec] = parseFMIFile("graph_output_fixed.fmi");
    if (parsedNodeMap.empty() || parsedAdjListVec.empty()) {
        std::cerr << "[ERROR] Graph konnte nicht geladen werden.\n";
        return 1;
    }
    nodeMap = std::move(parsedNodeMap);
    adjListVec = std::move(parsedAdjListVec);

    // Reverse-Graph bauen
    std::vector<std::vector<std::pair<int, double>>> adjListVecRev(adjListVec.size());
    for (int u = 0; u < (int)adjListVec.size(); ++u) {
        for (const auto& [v, w] : adjListVec[u]) {
            if (v >= 0 && v < (int)adjListVecRev.size())
                adjListVecRev[v].emplace_back(u, w);
        }
    }
    // Ungerichtet: adjListVecRev = adjListVec

    // Landmark-IDs (nur forward)
    landmarkIds.clear();
    std::string tableFolder = "tables_best_K220";

    for (auto& de : std::filesystem::directory_iterator(tableFolder)) {
        if (!de.is_regular_file()) continue;
        const std::string fname = de.path().filename().string();

        if (fname.rfind("L_rev_", 0) == 0) continue;

        int id = -1;
        auto dot = fname.find_last_of('.');
        if (dot == std::string::npos) continue;

        if (fname.rfind("L_fwd_", 0) == 0) {
            id = std::stoi(fname.substr(6, dot - 6));
        }
        else if (fname.rfind("L_", 0) == 0) {
            id = std::stoi(fname.substr(2, dot - 2));
        }
        else {
            continue;
        }
        landmarkIds.push_back(id);
    }

    std::sort(landmarkIds.begin(), landmarkIds.end());
    landmarkIds.erase(std::unique(landmarkIds.begin(), landmarkIds.end()), landmarkIds.end());
    std::cout << "[INFO] Landmark-IDs im Ordner (forward-only): " << landmarkIds.size() << "\n";

    // Tabellen laden (forward); backward = forward
    tables_forward = loadAllLandmarkTables(tableFolder, landmarkIds, /*reverse=*/false);
    const auto& tables_backward_ref = tables_forward;

    // ALT-Globals
    size_t maxIndex = adjListVec.size();
    for (size_t u = 0; u < adjListVec.size(); ++u) {
        for (const auto& e : adjListVec[u]) {
            int v = e.first;
            if (v >= 0) maxIndex = std::max(maxIndex, (size_t)v + 1);
        }
    }
    for (const auto& kv : nodeMap) {
        maxIndex = std::max(maxIndex, (size_t)kv.first + 1);
    }
    initAltGlobals(maxIndex);

    // Warmup
    auto warm_tables = [](const std::vector<std::vector<double>>& T) {
        volatile double sink = 0;
        const size_t step = 4096 / sizeof(double);
        for (const auto& col : T) for (size_t i = 0; i < col.size(); i += step) sink += col[i];
        (void)sink;
        };
    warm_tables(tables_forward);

    auto warm_reverse_adj = [](const auto& rev) {
        volatile double sink = 0;
        for (const auto& bucket : rev) if (!bucket.empty()) sink += bucket[0].second;
        (void)sink;
        };
    auto warm_tls = [&](size_t N) {
        const double INF = std::numeric_limits<double>::infinity();
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        // thread_local ALT-Buffer
        if (g_bi_dist_b.size() < N) g_bi_dist_b.resize(N, INF);
        if (g_bi_next_b.size() < N) g_bi_next_b.resize(N, -1);
        g_bi_touched_b.clear();

        if (bi_forwardHeuristic.size() < N) bi_forwardHeuristic.resize(N, NaN);
        if (bi_forwardBestLandmark.size() < N) bi_forwardBestLandmark.resize(N, -1);
        if (bi_backwardHeuristic.size() < N) bi_backwardHeuristic.resize(N, NaN);
        if (bi_backwardBestLandmark.size() < N) bi_backwardBestLandmark.resize(N, -1);
        bi_forwardTouched.clear();
        bi_backwardTouched.clear();
        };
    auto warm_bi_once = [&]() {
        int s = -1, t = -1;
        for (size_t i = 0; i < adjListVec.size(); ++i) {
            if (!adjListVec[i].empty()) { s = (int)i; t = adjListVec[i][0].first; break; }
        }
        if (s >= 0 && t >= 0 && nodeMap.count(s) && nodeMap.count(t)) {
            std::vector<double> dist; std::vector<int> prev;
            (void)alt_bidirectional(
                s, t,
                adjListVec, adjListVecRev,
                tables_forward, tables_backward_ref,
                landmarkIds, nodeMap,
                dist, prev
            );
        }
        };

    // Per-Thread Warmup
    const unsigned THREADS = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> warms; warms.reserve(THREADS);
    const size_t Nmax = std::max(adjListVec.size(), adjListVecRev.size());

    for (unsigned i = 0; i < THREADS; ++i) {
        warms.emplace_back([&] {
            warm_tables(tables_forward);
            warm_reverse_adj(adjListVecRev);
            warm_tls(Nmax);
            warm_bi_once();
            });
    }
    for (auto& th : warms) th.join();

    // Crow starten
    crow::SimpleApp app;
    using namespace crow;

    app.route_dynamic("/")
        .methods("GET"_method)
        ([](const crow::request&, crow::response& res) {
        std::ifstream file("static/index.html");
        if (!file.is_open()) {
            std::cerr << "Failed to open static/index.html. CWD=" << std::filesystem::current_path() << std::endl;
        }
        std::stringstream buffer; buffer << file.rdbuf();
        res.set_header("Content-Type", "text/html");
        res.write(buffer.str());
        res.end();
            });

    CROW_ROUTE(app, "/route")
        .methods("GET"_method)
        ([&](const crow::request& req) {
        auto source_str = req.url_params.get("source");
        auto target_str = req.url_params.get("target");
        crow::json::wvalue result;

        if (!source_str || !target_str) {
            result["found"] = false;
            return crow::response{ result };
        }

        int source = std::stoi(source_str);
        int target = std::stoi(target_str);

        if (nodeMap.find(source) == nodeMap.end() || nodeMap.find(target) == nodeMap.end()) {
            result["found"] = false;
            return crow::response{ result };
        }

        using Clock = std::chrono::high_resolution_clock;
        auto ms = [](auto dt) { return std::chrono::duration_cast<std::chrono::milliseconds>(dt).count(); };

        // ALT OR (Pfad, Metriken, used landmarks)
        {
            std::vector<double> dist_or;
            std::vector<int>    prev_or;

            auto t0 = Clock::now();
            bool found_or = alt_uniOrBi(
                source, target,
                adjListVec, adjListVecRev,
                tables_forward, tables_backward_ref,
                landmarkIds, nodeMap,
                dist_or, prev_or,
                QueryBucket::Unknown, /*log_decision=*/true
            );
            auto t1 = Clock::now();
            const auto elapsed = ms(t1 - t0);

            result["alt_or_time_ms"] = elapsed;
            result["alt_or_found"] = found_or;

            if (found_or && target < (int)dist_or.size() && std::isfinite(dist_or[target])) {
                const double d = dist_or[target];
                result["alt_or_dist"] = d;

                auto path = reconstruct_path(source, target, prev_or);
                std::vector<crow::json::wvalue> nodes; nodes.reserve(path.size());
                for (int id : path) {
                    const auto& n = nodeMap.at(id);
                    crow::json::wvalue jn;
                    jn["id"] = id;
                    jn["lat"] = n.lat;
                    jn["lon"] = n.lon;
                    nodes.push_back(std::move(jn));
                }
                result["path_nodes"] = std::move(nodes);
                result["found"] = true;
            }
            else {
                result["found"] = false;
            }

            // used landmarks
            std::vector<crow::json::wvalue> lm_arr;
            lm_arr.reserve(g_alt_lastUsedLandmarks.size());
            for (int lid : g_alt_lastUsedLandmarks) {
                auto it = nodeMap.find(lid);
                if (it == nodeMap.end()) continue;
                const auto& nd = it->second;
                crow::json::wvalue jlm;
                jlm["id"] = lid;
                jlm["lat"] = nd.lat;
                jlm["lon"] = nd.lon;
                lm_arr.push_back(std::move(jlm));
            }
            result["used_landmarks"] = std::move(lm_arr);
        }

        // A* (Metriken)
        {
            std::vector<double> dist_astar; std::vector<int> prev_astar;
            auto t0 = Clock::now();
            bool ok = astar(source, target, adjListVec, dist_astar, prev_astar, nodeMap);
            auto t1 = Clock::now();
            result["astar_found"] = ok;
            result["astar_time_ms"] = ms(t1 - t0);
            if (ok && target < (int)dist_astar.size() && std::isfinite(dist_astar[target])) {
                result["astar_dist"] = dist_astar[target];
            }
        }

        // ALT uni (Metriken)
        {
            std::vector<double> dist_uni; std::vector<int> prev_uni;
            auto t0 = Clock::now();
            bool ok = alt(source, target, adjListVec, tables_forward, tables_backward_ref, landmarkIds, nodeMap, dist_uni, prev_uni);
            auto t1 = Clock::now();
            result["alt_uni_found"] = ok;
            result["alt_uni_time_ms"] = ms(t1 - t0);
            if (ok && target < (int)dist_uni.size() && std::isfinite(dist_uni[target])) {
                result["alt_uni_dist"] = dist_uni[target];
            }
        }

        // ALT bi (Metriken)
        {
            std::vector<double> dist_bi; std::vector<int> prev_bi;
            auto t0 = Clock::now();
            bool ok = alt_bidirectional(source, target, adjListVec, adjListVecRev,
                tables_forward, tables_backward_ref,
                landmarkIds, nodeMap, dist_bi, prev_bi);
            auto t1 = Clock::now();
            result["alt_bi_found"] = ok;
            result["alt_bi_time_ms"] = ms(t1 - t0);
            if (ok && target < (int)dist_bi.size() && std::isfinite(dist_bi[target])) {
                result["alt_bi_dist"] = dist_bi[target];
            }
        }

        return crow::response{ result };
            });

    CROW_ROUTE(app, "/<string>")
        ([](const crow::request&, crow::response& res, std::string filename) {
        std::string path = "static/" + filename;
        std::ifstream file(path);
        if (!file) { res.code = 404; res.end(); return; }
        std::stringstream buffer; buffer << file.rdbuf();
        if (filename.find(".css") != std::string::npos) {
            res.set_header("Content-Type", "text/css");
        }
        else if (filename.find(".js") != std::string::npos) {
            res.set_header("Content-Type", "application/javascript");
        }
        res.write(buffer.str());
        res.end();
            });

    // App-Start mit fester Threadzahl
    app.port(18080).concurrency(THREADS).run();
    return 0;
}

int main(int argc, char* argv[]) {
    std::cout << "[DEBUG] CWD: " << std::filesystem::current_path() << std::endl;

    // Graph laden
    auto [parsedNodeMap, parsedAdjListVec] = parseFMIFile("graph_output_fixed.fmi");
    if (parsedNodeMap.empty() || parsedAdjListVec.empty()) {
        std::cerr << "[ERROR] Graph konnte nicht geladen werden.\n";
        return 1;
    }
    nodeMap = std::move(parsedNodeMap);
    adjListVec = std::move(parsedAdjListVec);

    // Reverse-Graph bauen
    std::vector<std::vector<std::pair<int, double>>> adjListVecRev(adjListVec.size());
    for (int u = 0; u < (int)adjListVec.size(); ++u) {
        for (const auto& [v, w] : adjListVec[u]) {
            if (v >= 0 && v < (int)adjListVecRev.size())
                adjListVecRev[v].emplace_back(u, w);
        }
    }
    // Ungerichtet: adjListVecRev = adjListVec

    // Landmark-IDs (nur forward)
    landmarkIds.clear();
    std::string tableFolder = "tables_best_K220";

    for (auto& de : std::filesystem::directory_iterator(tableFolder)) {
        if (!de.is_regular_file()) continue;
        const std::string fname = de.path().filename().string();

        if (fname.rfind("L_rev_", 0) == 0) continue;

        int id = -1;
        auto dot = fname.find_last_of('.');
        if (dot == std::string::npos) continue;

        if (fname.rfind("L_fwd_", 0) == 0) {
            id = std::stoi(fname.substr(6, dot - 6));
        }
        else if (fname.rfind("L_", 0) == 0) {
            id = std::stoi(fname.substr(2, dot - 2));
        }
        else {
            continue;
        }
        landmarkIds.push_back(id);
    }

    std::sort(landmarkIds.begin(), landmarkIds.end());
    landmarkIds.erase(std::unique(landmarkIds.begin(), landmarkIds.end()), landmarkIds.end());
    std::cout << "[INFO] Landmark-IDs im Ordner (forward-only): " << landmarkIds.size() << "\n";

    // Tabellen laden (forward); backward = forward
    tables_forward = loadAllLandmarkTables(tableFolder, landmarkIds, /*reverse=*/false);
    const auto& tables_backward_ref = tables_forward;

    // ALT-Globals
    size_t maxIndex = adjListVec.size();
    for (size_t u = 0; u < adjListVec.size(); ++u) {
        for (const auto& e : adjListVec[u]) {
            int v = e.first;
            if (v >= 0) maxIndex = std::max(maxIndex, (size_t)v + 1);
        }
    }
    for (const auto& kv : nodeMap) {
        maxIndex = std::max(maxIndex, (size_t)kv.first + 1);
    }
    initAltGlobals(maxIndex);

    // Warmup
    auto warm_tables = [](const std::vector<std::vector<double>>& T) {
        volatile double sink = 0;
        const size_t step = 4096 / sizeof(double);
        for (const auto& col : T) for (size_t i = 0; i < col.size(); i += step) sink += col[i];
        (void)sink;
        };
    warm_tables(tables_forward);

    auto warm_reverse_adj = [](const auto& rev) {
        volatile double sink = 0;
        for (const auto& bucket : rev) if (!bucket.empty()) sink += bucket[0].second;
        (void)sink;
        };
    auto warm_tls = [&](size_t N) {
        const double INF = std::numeric_limits<double>::infinity();
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        // thread_local ALT-Buffer
        if (g_bi_dist_b.size() < N) g_bi_dist_b.resize(N, INF);
        if (g_bi_next_b.size() < N) g_bi_next_b.resize(N, -1);
        g_bi_touched_b.clear();

        if (bi_forwardHeuristic.size() < N) bi_forwardHeuristic.resize(N, NaN);
        if (bi_forwardBestLandmark.size() < N) bi_forwardBestLandmark.resize(N, -1);
        if (bi_backwardHeuristic.size() < N) bi_backwardHeuristic.resize(N, NaN);
        if (bi_backwardBestLandmark.size() < N) bi_backwardBestLandmark.resize(N, -1);
        bi_forwardTouched.clear();
        bi_backwardTouched.clear();
        };
    auto warm_bi_once = [&]() {
        int s = -1, t = -1;
        for (size_t i = 0; i < adjListVec.size(); ++i) {
            if (!adjListVec[i].empty()) { s = (int)i; t = adjListVec[i][0].first; break; }
        }
        if (s >= 0 && t >= 0 && nodeMap.count(s) && nodeMap.count(t)) {
            std::vector<double> dist; std::vector<int> prev;
            (void)alt_bidirectional(
                s, t,
                adjListVec, adjListVecRev,
                tables_forward, tables_backward_ref,
                landmarkIds, nodeMap,
                dist, prev
            );
        }
        };

    // Per-Thread Warmup
    const unsigned THREADS = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> warms; warms.reserve(THREADS);
    const size_t Nmax = std::max(adjListVec.size(), adjListVecRev.size());

    for (unsigned i = 0; i < THREADS; ++i) {
        warms.emplace_back([&] {
            warm_tables(tables_forward);
            warm_reverse_adj(adjListVecRev);
            warm_tls(Nmax);
            warm_bi_once();
            });
    }
    for (auto& th : warms) th.join();

    // Crow starten
    crow::SimpleApp app;
    using namespace crow;

    app.route_dynamic("/")
        .methods("GET"_method)
        ([](const crow::request&, crow::response& res) {
        std::ifstream file("static/index.html");
        if (!file.is_open()) {
            std::cerr << "Failed to open static/index.html. CWD=" << std::filesystem::current_path() << std::endl;
        }
        std::stringstream buffer; buffer << file.rdbuf();
        res.set_header("Content-Type", "text/html");
        res.write(buffer.str());
        res.end();
            });

    CROW_ROUTE(app, "/route")
        .methods("GET"_method)
        ([&](const crow::request& req) {
        auto source_str = req.url_params.get("source");
        auto target_str = req.url_params.get("target");
        crow::json::wvalue result;

        if (!source_str || !target_str) {
            result["found"] = false;
            return crow::response{ result };
        }

        int source = std::stoi(source_str);
        int target = std::stoi(target_str);

        if (nodeMap.find(source) == nodeMap.end() || nodeMap.find(target) == nodeMap.end()) {
            result["found"] = false;
            return crow::response{ result };
        }

        using Clock = std::chrono::high_resolution_clock;
        auto ms = [](auto dt) { return std::chrono::duration_cast<std::chrono::milliseconds>(dt).count(); };

        // ALT OR (Pfad, Metriken, used landmarks)
        {
            std::vector<double> dist_or;
            std::vector<int>    prev_or;

            auto t0 = Clock::now();
            bool found_or = alt_uniOrBi(
                source, target,
                adjListVec, adjListVecRev,
                tables_forward, tables_backward_ref,
                landmarkIds, nodeMap,
                dist_or, prev_or,
                QueryBucket::Unknown, /*log_decision=*/true
            );
            auto t1 = Clock::now();
            const auto elapsed = ms(t1 - t0);

            result["alt_or_time_ms"] = elapsed;
            result["alt_or_found"] = found_or;

            if (found_or && target < (int)dist_or.size() && std::isfinite(dist_or[target])) {
                const double d = dist_or[target];
                result["alt_or_dist"] = d;

                auto path = reconstruct_path(source, target, prev_or);
                std::vector<crow::json::wvalue> nodes; nodes.reserve(path.size());
                for (int id : path) {
                    const auto& n = nodeMap.at(id);
                    crow::json::wvalue jn;
                    jn["id"] = id;
                    jn["lat"] = n.lat;
                    jn["lon"] = n.lon;
                    nodes.push_back(std::move(jn));
                }
                result["path_nodes"] = std::move(nodes);
                result["found"] = true;
            }
            else {
                result["found"] = false;
            }

            // used landmarks
            std::vector<crow::json::wvalue> lm_arr;
            lm_arr.reserve(g_alt_lastUsedLandmarks.size());
            for (int lid : g_alt_lastUsedLandmarks) {
                auto it = nodeMap.find(lid);
                if (it == nodeMap.end()) continue;
                const auto& nd = it->second;
                crow::json::wvalue jlm;
                jlm["id"] = lid;
                jlm["lat"] = nd.lat;
                jlm["lon"] = nd.lon;
                lm_arr.push_back(std::move(jlm));
            }
            result["used_landmarks"] = std::move(lm_arr);
        }

        return crow::response{ result };
            });

    CROW_ROUTE(app, "/<string>")
        ([](const crow::request&, crow::response& res, std::string filename) {
        std::string path = "static/" + filename;
        std::ifstream file(path);
        if (!file) { res.code = 404; res.end(); return; }
        std::stringstream buffer; buffer << file.rdbuf();
        if (filename.find(".css") != std::string::npos) {
            res.set_header("Content-Type", "text/css");
        }
        else if (filename.find(".js") != std::string::npos) {
            res.set_header("Content-Type", "application/javascript");
        }
        res.write(buffer.str());
        res.end();
            });

    // App-Start mit fester Threadzahl
    app.port(18080).concurrency(THREADS).run();
    return 0;
}
