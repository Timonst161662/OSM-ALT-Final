#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>      // std::reverse
#include <random>         // std::mt19937
#include <filesystem>     // debug output / dirs
#include <unordered_set>  // std::unordered_set
#include <chrono>         // timing
#include <thread>
#include <atomic>
#include <nlohmann/json.hpp>

#include "Geo.h"          // << neu: Geo-Header (hält Point/Hilfsfunktionen im Namespace geo)
#include "GraphUtils.h"
#include "ALT.h"
#include "AStar.h"
#include "Dijkstra.h"

using json = nlohmann::json;
using namespace std;

// ===================== OPTION B: Stratifizierte Query-Erzeugung =====================
#include <array>
#include <random>

// Für Option B benutzen wir eine Proxy-Distanz (Haversine) auf Basis von nodeMap:
struct QueryPair {
    int s;
    int t;
    double d_proxy; // Luftlinie in Metern (nur zum Binning/Sortieren)
};

static inline double dist_proxy(int a, int b, const std::unordered_map<int, Node>& nm) {
    const auto& A = nm.at(a);
    const auto& B = nm.at(b);
    // Falls geo::haversine nicht im Scope ist, ggf. qualifizieren:
    return geo::haversine(A.lat, A.lon, B.lat, B.lon);
}

// --- Hilfsfunktion: Quantile-Schätzer für Query-Distanzen ---
static std::array<double, 4> estimate_quantiles(
    const std::vector<int>& cand,
    const std::unordered_map<int, Node>& nm,
    std::mt19937& rng,
    int sample_count)
{
    std::uniform_int_distribution<size_t> pick(0, cand.size() - 1);
    std::vector<double> ds;
    ds.reserve(sample_count);

    for (int i = 0; i < sample_count; ++i) {
        int s = (int)cand[pick(rng)];
        int t = (int)cand[pick(rng)];
        if (s == t) continue;
        double d = haversine_source_target_m(s, t, nm);
        if (std::isfinite(d)) ds.push_back(d);
    }

    if (ds.empty()) return { 0,0,0,0 };

    std::sort(ds.begin(), ds.end());
    auto q = [&](double p) -> double {
        size_t k = (size_t)std::clamp<size_t>(std::llround(p * (ds.size() - 1)), 0, ds.size() - 1);
        return ds[k];
        };

    return { q(0.20), q(0.40), q(0.60), q(0.80) };
}



// Baut numQueries Queries gemäß Quoten je Distanz-Bucket (XS,S,M,L,XL = 20/30/30/15/5)
static std::vector<QueryPair> build_queries_stratified(
    int numQueries,
    const std::vector<int>& cand,
    const std::unordered_map<int, Node>& nm,
    std::mt19937& rng)
{
    if (cand.size() < 2 || numQueries <= 0) return {};

    // 1) Quantile
    auto Q = estimate_quantiles(cand, nm, rng, std::max(1000, numQueries * 10));

    // 2) Zielverteilung
    const double frac[5] = { 0.20, 0.30, 0.30, 0.15, 0.05 };
    int target[5]; int sum = 0;
    for (int i = 0; i < 5; ++i) { target[i] = (int)std::round(frac[i] * numQueries); sum += target[i]; }
    while (sum > numQueries) { for (int i = 4; i >= 0 && sum > numQueries; --i) if (target[i] > 0) { --target[i]; --sum; } }
    while (sum < numQueries) { for (int i = 0; i < 5 && sum < numQueries; ++i) { ++target[i]; ++sum; } }

    auto bin_of = [&](double d)->int {
        if (d < Q[0]) return 0; // XS
        if (d < Q[1]) return 1; // S
        if (d < Q[2]) return 2; // M
        if (d < Q[3]) return 3; // L
        return 4;              // XL
        };

    std::uniform_int_distribution<size_t> pick(0, cand.size() - 1);
    std::vector<QueryPair> out; out.reserve(numQueries);
    int filled[5]{ 0,0,0,0,0 };

    // 3) Erstpass: fülle Bins gemäß Ziel
    const int MAX_TRIES = std::max(numQueries * 50, 10000);
    for (int tries = 0; (int)out.size() < numQueries && tries < MAX_TRIES; ++tries) {
        int s = (int)cand[pick(rng)], t = (int)cand[pick(rng)];
        if (s == t) continue;
        double d = dist_proxy(s, t, nm);
        int b = bin_of(d);
        if (filled[b] < target[b]) {
            out.push_back({ s, t, d });
            ++filled[b];
        }
    }
    // 4) Rest auffüllen (falls ein Bin „zu schwer“ war)
    while ((int)out.size() < numQueries) {
        int s = (int)cand[pick(rng)], t = (int)cand[pick(rng)];
        if (s == t) continue;
        out.push_back({ s, t, dist_proxy(s, t, nm) });
    }

    std::cout << "[QS] q20=" << (long long)Q[0]
        << " q40=" << (long long)Q[1]
        << " q60=" << (long long)Q[2]
        << " q80=" << (long long)Q[3]
        << " | target XS/S/M/L/XL="
        << target[0] << "/" << target[1] << "/" << target[2] << "/"
        << target[3] << "/" << target[4] << "\n";
    std::cout << "[QS] filled XS/S/M/L/XL="
        << filled[0] << "/" << filled[1] << "/" << filled[2] << "/"
        << filled[3] << "/" << filled[4] << "\n";
    return out;
}


int main_Benchmark1K() {
    // === Konfiguration ===
    const int MAX_QUERIES = 1000;
    std::string tableFolder = "tables_best_K220";
    namespace fs = std::filesystem;

    // Helper
    auto safeAt = [](const std::vector<double>& v, int idx) -> double {
        return (idx >= 0 && idx < (int)v.size() && std::isfinite(v[idx]))
            ? v[idx]
            : std::numeric_limits<double>::quiet_NaN();
        };

    std::cout << "[DEBUG] CWD: " << fs::current_path() << "\n";

    // === Schritt 1: Graph laden ===
    auto [nodeMap, adjListVec] = parseFMIFile("graph_output_fixed.fmi");
    if (nodeMap.empty() || adjListVec.empty()) {
        std::cerr << "[ERROR] Graph konnte nicht geladen werden.\n";
        return 1;
    }

    // === Schritt 2: Reverse-Graph (Vektor) aus forward ableiten (für ALT_bi) ===
    std::vector<std::vector<std::pair<int, double>>> adjListVecRev(adjListVec.size());
    for (int u = 0; u < (int)adjListVec.size(); ++u) {
        for (const auto& [v, w] : adjListVec[u]) {
            if (v >= 0 && v < (int)adjListVecRev.size())
                adjListVecRev[v].emplace_back(u, w);
        }
    }

    // === Schritt 3: Landmark-IDs & Tabellen laden (wie in main_web) ===
    if (!fs::exists(tableFolder) || !fs::is_directory(tableFolder)) {
        std::cerr << "[ERROR] LM-Ordner nicht gefunden: " << tableFolder
            << " (CWD=" << fs::current_path() << ")\n";
        return 1;
    }

    std::vector<int> landmarkIds;
    landmarkIds.reserve(512);
    int files_seen = 0;
    for (const auto& de : fs::directory_iterator(tableFolder)) {
        if (!de.is_regular_file()) continue;
        const std::string fname = de.path().filename().string();
        ++files_seen;

        if (fname.rfind("L_rev_", 0) == 0) continue; // Reverse explizit skippen

        int id = -1;
        auto dot = fname.find_last_of('.');
        if (dot == std::string::npos) continue;

        if (fname.rfind("L_fwd_", 0) == 0) {                 // "L_fwd_<id>.tbl"
            id = std::stoi(fname.substr(6, dot - 6));
        }
        else if (fname.rfind("L_", 0) == 0) {              // "L_<id>.tbl"
            id = std::stoi(fname.substr(2, dot - 2));
        }
        else {
            continue;
        }
        landmarkIds.push_back(id);
    }
    std::sort(landmarkIds.begin(), landmarkIds.end());
    landmarkIds.erase(std::unique(landmarkIds.begin(), landmarkIds.end()), landmarkIds.end());

    std::cout << "[INFO] LM-Dateien gesehen: " << files_seen
        << " | akzeptierte Landmark-IDs: " << landmarkIds.size() << "\n";
    if (landmarkIds.empty()) {
        std::cerr << "[ERROR] Keine Landmark-IDs erkannt. Prüfe Dateinamen (L_fwd_*.tbl / L_*.tbl) "
            "und den Ordnerpfad relativ zum CWD.\n";
        return 1;
    }

    auto tables_forward = loadAllLandmarkTables(tableFolder, landmarkIds, /*reverse=*/false);
    if (tables_forward.size() != landmarkIds.size()) {
        std::cerr << "[ERROR] Tabellen/ID-Mismatch: tables=" << tables_forward.size()
            << " vs ids=" << landmarkIds.size() << "\n";
        return 1;
    }
    const auto& tables_backward = tables_forward; // Alias, kein Extra-RAM
    std::cout << "[INFO] LM-Tables geladen: " << tables_forward.size() << " (FWD) | REV=Alias\n";

    // Sanity-Check: Spaltenlängen >= maxNodeId+1 ?
    size_t maxNodeId = 0;
    for (const auto& kv : nodeMap) maxNodeId = std::max(maxNodeId, (size_t)kv.first);
    const size_t expectedLen = maxNodeId + 1;
    for (size_t i = 0; i < tables_forward.size(); ++i) {
        const auto& col = tables_forward[i];
        if (col.size() < expectedLen) {
            std::cerr << "[ERROR] LM col " << i << " length=" << col.size()
                << " < expected " << expectedLen
                << " (Node-IDs größer als Tabellendimension?)\n";
            return 1;
        }
    }
    std::cout << "[CHECK] LM tables ok: " << tables_forward.size()
        << " Landmarken, column length >= " << expectedLen << "\n";

    // === Warmup (identisch zu main_all5Algos) ===================================

// Unidirectional: Landmark-Tabellen seitenweise berühren
    auto warm_tables = [](const std::vector<std::vector<double>>& T) {
        volatile double sink = 0;
        const size_t step = 4096 / sizeof(double);
        for (const auto& col : T)
            for (size_t i = 0; i < col.size(); i += step) sink += col[i];
        (void)sink;
        };
    warm_tables(tables_forward);

    // Bidirectional: ein paar Reverse-Buckets anfassen (Cache-Touch)
    auto warm_reverse_adj = [](const auto& rev) {
        volatile double sink = 0;
        for (const auto& bucket : rev)
            if (!bucket.empty()) sink += bucket[0].second;
        (void)sink;
        };

    // TLS/Heuristik-Puffer auf Zielgröße bringen + touched-Listen leeren
    auto warm_tls = [&](size_t N) {
        const double INF = std::numeric_limits<double>::infinity();
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        if (g_bi_dist_b.size() < N) g_bi_dist_b.resize(N, INF);
        if (g_bi_next_b.size() < N) g_bi_next_b.resize(N, -1);
        g_bi_touched_b.clear();

        if (bi_forwardHeuristic.size() < N)     bi_forwardHeuristic.resize(N, NaN);
        if (bi_forwardBestLandmark.size() < N)  bi_forwardBestLandmark.resize(N, -1);
        if (bi_backwardHeuristic.size() < N)    bi_backwardHeuristic.resize(N, NaN);
        if (bi_backwardBestLandmark.size() < N) bi_backwardBestLandmark.resize(N, -1);
        bi_forwardTouched.clear();
        bi_backwardTouched.clear();
        };

    // Ein superkurzer ALT_bi-Aufruf (erste Kante) um PQ/Heuristik zu „berühren“
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
                tables_forward, tables_backward,
                landmarkIds, nodeMap,
                dist, prev
            );
        }
        };

    // Pro CPU-Thread einmal warm machen (wie in main_all5Algos)
    const unsigned THREADS = std::max(1u, std::thread::hardware_concurrency());
    const size_t Nmax = std::max(adjListVec.size(), adjListVecRev.size());

    std::vector<std::thread> warms; warms.reserve(THREADS);
    for (unsigned i = 0; i < THREADS; ++i) {
        warms.emplace_back([&] {
            warm_tables(tables_forward);
            warm_reverse_adj(adjListVecRev);
            warm_tls(Nmax);
            warm_bi_once();
            });
    }
    for (auto& th : warms) th.join();


    // === ALT-Globals initialisieren ===
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

    // === Schritt 4: Queries zufällig erzeugen (echtes Uniform-Random) ===
    const bool USE_RANDOM = true;            // nur zur Doku; hier immer true
    const int  MAX_TRIES = MAX_QUERIES * 20; // Obergrenze gegen Endlosschleifen

    // optional reproduzierbar machen (Seed via Umgebungsvariable OSM_SEED oder fallback)
    unsigned long seed = [] {
#ifdef _WIN32
        // Windows-sicheres getenv
        char* buf = nullptr;
        size_t sz = 0;
        if (_dupenv_s(&buf, &sz, "OSM_SEED") == 0 && buf) {
            unsigned long val = std::strtoul(buf, nullptr, 10);
            free(buf);
            return val;
        }
        return (unsigned long)std::chrono::high_resolution_clock::now().time_since_epoch().count();
#else
        if (const char* s = std::getenv("OSM_SEED")) return std::strtoul(s, nullptr, 10);
        return (unsigned long)std::chrono::high_resolution_clock::now().time_since_epoch().count();
#endif
        }();

    // --- Hier wirklich den rng definieren ---
    std::mt19937 rng(static_cast<uint32_t>(seed));
    std::cout << "[INFO] Random seed = " << seed << "\n";



    // Kandidatenliste aus den tatsächlichen Node-IDs aufbauen
    std::vector<int> nodes; nodes.reserve(nodeMap.size());
    for (const auto& kv : nodeMap) nodes.push_back(kv.first);
    if (nodes.empty()) { std::cerr << "[ERROR] Keine Nodes vorhanden.\n"; return 1; }

    std::uniform_int_distribution<size_t> pick(0, nodes.size() - 1);

    // Doppelte (s,t) vermeiden: 64-bit Key (s<<32 | t)
    auto pack_pair = [](int s, int t) -> uint64_t {
        return (uint64_t)(uint32_t)s << 32 | (uint32_t)t;
        };
    std::unordered_set<uint64_t> seen; seen.reserve(MAX_QUERIES * 2);

    // optional: Mindestluftlinie in Metern, um triviale Nachbarn zu vermeiden (0 = aus)
    const double MIN_HAVERSINE_M = 0.0; // z.B. 1000.0 für >1km

    std::vector<std::pair<int, int>> queries;
    queries.reserve(MAX_QUERIES);

    int tries = 0;
    while ((int)queries.size() < MAX_QUERIES && tries < MAX_TRIES) {
        ++tries;
        int s = nodes[pick(rng)];
        int t = nodes[pick(rng)];
        if (s == t) continue;

        // Distanz-Filter (optional)
        if (MIN_HAVERSINE_M > 0.0) {
            double d = haversine_source_target_m(s, t, nodeMap);
            if (!(std::isfinite(d) && d >= MIN_HAVERSINE_M)) continue;
        }

        // unique halten
        if (!seen.insert(pack_pair(s, t)).second) continue;

        queries.emplace_back(s, t);
    }

    if (queries.empty()) {
        std::cerr << "[ERROR] Konnte keine Random-Queries erzeugen (Filter zu streng?).\n";
        return 1;
    }
    std::cout << "[INFO] Random Queries erzeugt: " << queries.size()
        << " (von max " << MAX_QUERIES << ", tries=" << tries << ")\n";



    // === CSV-Ausgabe vorbereiten ===
    std::ofstream csv("routing_results.csv");
    if (!csv) {
        std::cerr << "[ERROR] Konnte 'routing_results.csv' nicht zum Schreiben öffnen.\n";
        return 1;
    }
    csv << "QueryID,Algorithm,SourceID,SourceLat,SourceLon,TargetID,TargetLat,TargetLon,Distance,TimeMs,PathNodeIDs,PathCoords\n";

    // --- Diagnose-Potenzial: max |d(L,t) - d(L,s)| über 64 random Paare ---
    auto randId = [&]() {
        auto it = nodeMap.begin();
        std::advance(it, std::min<size_t>(rand() % nodeMap.size(), nodeMap.size() - 1));
        return it->first;
        };
    double avgMaxPhi = 0.0; int samples = 0;
    for (int k = 0; k < 64; ++k) {
        int s = randId(), t = randId();
        if (s == t) continue;
        double maxPhi = 0.0;
        for (size_t j = 0; j < landmarkIds.size(); ++j) {
            const auto& col = tables_forward[j];
            double dst = col[t], dss = col[s];
            if (!std::isfinite(dst) || !std::isfinite(dss)) continue;
            maxPhi = std::max(maxPhi, std::abs(dst - dss));
        }
        if (maxPhi > 0) { avgMaxPhi += maxPhi; ++samples; }
    }
    if (samples == 0) {
        std::cerr << "[ERROR] Keine sinnvollen Landmark-Potenziale (alle 0/NaN?). ALT ~ A*\n";
    }
    else {
        avgMaxPhi /= samples;
        std::cout << "[CHECK] avg max-|phi| über 64 Paare: " << avgMaxPhi << " (Meter)\n";
        if (avgMaxPhi < 1e3) {
            std::cerr << "[WARN] Potenziale extrem klein – ALT-Heuristik schwach/ausgeschaltet.\n";
        }
    }

    fs::create_directories("debug/routes");
    using Clock = std::chrono::high_resolution_clock;
    auto ms = [](auto dt) { return std::chrono::duration_cast<std::chrono::milliseconds>(dt).count(); };

    // ===== Sammelstruktur für Summary =====
    struct Row {
        double d_m;
        bool ok_dij;    double t_dij_ms;
        bool ok_astar;  double t_astar_ms;
        bool ok_alt_uni;double t_alt_uni_ms;
        bool ok_alt_bi; double t_alt_bi_ms;
        bool ok_alt_or; double t_alt_or_ms;
    };
    std::vector<Row> rows; rows.reserve(std::min<int>(MAX_QUERIES, (int)queries.size()));

    // ===== Hauptschleife =====
    int qi = 0;
    for (const auto& [source, target] : queries) {
        if (qi >= MAX_QUERIES) break;
        ++qi;

        // ===== Messungen =====
        const auto& srcNode = nodeMap.at(source);
        const auto& tgtNode = nodeMap.at(target);
        const double d_m = haversine_source_target_m(source, target, nodeMap);

        std::vector<double> dist_dij;  std::vector<int> prev_dij;
        auto t0 = Clock::now();
        bool found_dij = dijkstra_early_exit(source, target, adjListVec, dist_dij, prev_dij);
        auto t1 = Clock::now();
        auto time_dij = ms(t1 - t0);
        if (found_dij)
            csv << qi << ",Dijkstra," << source << "," << srcNode.lat << "," << srcNode.lon << ","
            << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
            << safeAt(dist_dij, target) << "," << time_dij << ",," << "\n";

        std::vector<double> dist_astar; std::vector<int> prev_astar;
        auto t2 = Clock::now();
        bool found_astar = astar(source, target, adjListVec, dist_astar, prev_astar, nodeMap);
        auto t3 = Clock::now();
        auto time_astar = ms(t3 - t2);
        if (found_astar)
            csv << qi << ",A*," << source << "," << srcNode.lat << "," << srcNode.lon << ","
            << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
            << safeAt(dist_astar, target) << "," << time_astar << ",," << "\n";

        std::vector<double> dist_alt_uni; std::vector<int> prev_alt_uni;
        auto t4 = Clock::now();
        bool found_alt_uni = alt(source, target,
            adjListVec,
            tables_forward, tables_backward,
            landmarkIds, nodeMap,
            dist_alt_uni, prev_alt_uni);
        auto t5 = Clock::now();
        auto time_alt_uni = ms(t5 - t4);
        if (found_alt_uni) {
            auto path = rebuildPath(source, target, prev_alt_uni);
            csv << qi << ",ALT_uni,"
                << source << "," << srcNode.lat << "," << srcNode.lon << ","
                << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
                << safeAt(dist_alt_uni, target) << "," << time_alt_uni << ","
                << "\"" << pathIdsCsv(path) << "\","
                << "\"" << pathCoordsCsv(path, nodeMap) << "\"\n";
        }

        std::vector<double> dist_alt_bi; std::vector<int> prev_alt_bi;
        auto t6 = Clock::now();
        bool found_alt_bi = alt_bidirectional(source, target,
            adjListVec, adjListVecRev,
            tables_forward, tables_backward,
            landmarkIds, nodeMap,
            dist_alt_bi, prev_alt_bi);
        auto t7 = Clock::now();
        auto time_alt_bi = ms(t7 - t6);

        double dist_bi_val = std::numeric_limits<double>::quiet_NaN();
        if (target >= 0 && target < (int)dist_alt_bi.size() && std::isfinite(dist_alt_bi[target])) {
            dist_bi_val = dist_alt_bi[target];
        }
        else {
            std::cerr << "[ALT_bi] WARN: target=" << target
                << " size=" << dist_alt_bi.size()
                << " finite?=" << (target < (int)dist_alt_bi.size()
                    && std::isfinite(dist_alt_bi[target]))
                << "\n";
        }
        if (found_alt_bi) {
            auto path_bi = rebuildPath(source, target, prev_alt_bi);
            csv << qi << ",ALT_bi,"
                << source << "," << srcNode.lat << "," << srcNode.lon << ","
                << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
                << dist_bi_val << "," << time_alt_bi << ","
                << "\"" << pathIdsCsv(path_bi) << "\","
                << "\"" << pathCoordsCsv(path_bi, nodeMap) << "\"\n";
        }

        std::vector<double> dist_alt_or; std::vector<int> prev_alt_or;
        auto t8 = Clock::now();
        bool found_alt_or = alt_uniOrBi(source, target,
            adjListVec, adjListVecRev,
            tables_forward, tables_backward,
            landmarkIds, nodeMap,
            dist_alt_or, prev_alt_or,
            QueryBucket::Unknown, true);
        auto t9 = Clock::now();
        auto time_alt_or = ms(t9 - t8);

        if (found_alt_or) {
            const bool took_bi = (haversine_source_target_m(source, target, nodeMap) > OR_SWITCH_METERS);
            if (!took_bi) {
                auto path = rebuildPath(source, target, prev_alt_or);
                csv << qi << ",ALT_orUNI,"
                    << source << "," << srcNode.lat << "," << srcNode.lon << ","
                    << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
                    << safeAt(dist_alt_or, target) << "," << time_alt_or << ","
                    << "\"" << pathIdsCsv(path) << "\","
                    << "\"" << pathCoordsCsv(path, nodeMap) << "\"\n";
            }
            else {
                csv << qi << ",ALT_orBI,"
                    << source << "," << srcNode.lat << "," << srcNode.lon << ","
                    << target << "," << tgtNode.lat << "," << tgtNode.lon << ","
                    << safeAt(dist_alt_or, target) << "," << time_alt_or << ",,\n";
            }
        }

        // ===== Kompakter Abgabe-Log (Zeit + Länge in km + Knotenanzahl) =====
        auto count_nodes = [&](const std::vector<int>& prev, int s, int t) -> int {
            auto path = rebuildPath(s, t, prev);
            return (int)path.size();
            };
        auto km = [](double meters) -> double { return meters / 1000.0; };

        // Referenz-Distanz (für Header)
        double ref_dist_m = std::numeric_limits<double>::quiet_NaN();
        if (found_dij && std::isfinite(safeAt(dist_dij, target))) {
            ref_dist_m = safeAt(dist_dij, target);
        }
        else if (found_astar && std::isfinite(safeAt(dist_astar, target))) {
            ref_dist_m = safeAt(dist_astar, target);
        }
        else if (found_alt_uni && std::isfinite(safeAt(dist_alt_uni, target))) {
            ref_dist_m = safeAt(dist_alt_uni, target);
        }
        else if (found_alt_bi && std::isfinite(dist_bi_val)) {
            ref_dist_m = dist_bi_val;
        }
        else if (found_alt_or && std::isfinite(safeAt(dist_alt_or, target))) {
            ref_dist_m = safeAt(dist_alt_or, target);
        }

        // Pfadlängen (Knoten)
        int len_dij = found_dij ? count_nodes(prev_dij, source, target) : -1;
        int len_as = found_astar ? count_nodes(prev_astar, source, target) : -1;
        int len_uni = found_alt_uni ? count_nodes(prev_alt_uni, source, target) : -1;
        int len_bi = found_alt_bi ? count_nodes(prev_alt_bi, source, target) : -1;
        int len_or = found_alt_or ? count_nodes(prev_alt_or, source, target) : -1;

        // Distanzen (m)
        double dij_m = found_dij ? safeAt(dist_dij, target) : std::numeric_limits<double>::quiet_NaN();
        double as_m = found_astar ? safeAt(dist_astar, target) : std::numeric_limits<double>::quiet_NaN();
        double uni_m = found_alt_uni ? safeAt(dist_alt_uni, target) : std::numeric_limits<double>::quiet_NaN();
        double bi_m = found_alt_bi ? dist_bi_val : std::numeric_limits<double>::quiet_NaN();
        double or_m = found_alt_or ? safeAt(dist_alt_or, target) : std::numeric_limits<double>::quiet_NaN();

        // Header + kompakte Zeilen
        std::cout << "=============================================\n";
        std::cout << "[TEST #" << qi << "] Source=" << source << " Target=" << target << "\n";


        auto line = [&](const char* name, long long t_ms, double dist_m, int nodes) {
            if (!std::isfinite(dist_m) || nodes < 0) return;
            std::cout << std::left << std::setw(12) << name << ": "
                << std::right << std::setw(5) << t_ms << " ms | "
                << std::fixed << std::setprecision(3) << std::setw(10) << km(dist_m) << " km | "
                << nodes << " Knoten\n";
            };

        if (found_dij)    line("Dijkstra", (long long)time_dij, dij_m, len_dij);
        if (found_astar)  line("A*", (long long)time_astar, as_m, len_as);
        if (found_alt_uni)line("ALT uni", (long long)time_alt_uni, uni_m, len_uni);
        if (found_alt_bi) line("ALT bi", (long long)time_alt_bi, bi_m, len_bi);
        if (found_alt_or) line("ALT UniOrBi", (long long)time_alt_or, or_m, len_or);


        // === Messreihe für Summary sammeln ===
        rows.push_back({
            d_m,
            found_dij,   double(time_dij),
            found_astar, double(time_astar),
            found_alt_uni, double(time_alt_uni),
            found_alt_bi,  double(time_alt_bi),
            found_alt_or,  double(time_alt_or)
            });
    }

    csv.close();

    // ======= SUMMARY ALS .TXT (Overall + Quantil-Buckets aus den echten Distanzen) =======
    {
        if (rows.empty()) {
            std::cerr << "[WARN] Keine Messzeilen für Summary.\n";
        }
        else {
            // Quantile q20,q40,q60,q80 aus rows.d_m
            std::vector<double> ds; ds.reserve(rows.size());
            for (auto& r : rows) ds.push_back(r.d_m);

            auto nth = [&](double p)->double {
                size_t k = (size_t)std::clamp<size_t>(std::llround(p * (ds.size() - 1)), 0, ds.size() - 1);
                std::vector<double> tmp = ds;
                std::nth_element(tmp.begin(), tmp.begin() + k, tmp.end());
                return tmp[k];
                };
            const double q20 = nth(0.20), q40 = nth(0.40), q60 = nth(0.60), q80 = nth(0.80);

            auto bin_of = [&](double d)->int {
                if (d < q20) return 0;        // XS
                if (d < q40) return 1;        // S
                if (d < q60) return 2;        // M
                if (d < q80) return 3;        // L
                return 4;                     // XL
                };

            struct Agg { double sum = 0.0; int n = 0; };
            auto upd = [](Agg& a, bool ok, double v) {
                if (ok && std::isfinite(v)) { a.sum += v; a.n += 1; }
                };
            auto avg = [](const Agg& a) { return (a.n ? (a.sum / a.n) : std::numeric_limits<double>::quiet_NaN()); };

            // Overall
            Agg o_dij, o_astar, o_uni, o_bi, o_or;
            for (const auto& r : rows) {
                upd(o_dij, r.ok_dij, r.t_dij_ms);
                upd(o_astar, r.ok_astar, r.t_astar_ms);
                upd(o_uni, r.ok_alt_uni, r.t_alt_uni_ms);
                upd(o_bi, r.ok_alt_bi, r.t_alt_bi_ms);
                upd(o_or, r.ok_alt_or, r.t_alt_or_ms);
            }

            // Buckets
            struct BA { Agg dij, astar, uni, bi, orr; int count = 0; };
            std::array<BA, 5> B{};
            for (const auto& r : rows) {
                int b = bin_of(r.d_m);
                B[b].count++;
                upd(B[b].dij, r.ok_dij, r.t_dij_ms);
                upd(B[b].astar, r.ok_astar, r.t_astar_ms);
                upd(B[b].uni, r.ok_alt_uni, r.t_alt_uni_ms);
                upd(B[b].bi, r.ok_alt_bi, r.t_alt_bi_ms);
                upd(B[b].orr, r.ok_alt_or, r.t_alt_or_ms);
            }

            // Ausgabe
            std::ofstream txt("routing_summary.txt");
            if (!txt) {
                std::cerr << "[ERROR] Konnte 'routing_summary.txt' nicht schreiben.\n";
            }
            else {
                txt << std::fixed << std::setprecision(3);
                txt << "Routing Summary\n";
                txt << "================\n\n";
                txt << "Queries (verarbeitet): " << rows.size() << "\n";
                txt << "Quantile (Meter): q20=" << q20 << "  q40=" << q40
                    << "  q60=" << q60 << "  q80=" << q80 << "\n\n";

                txt << "Overall Averages (ms) — nur erfolgreiche Läufe je Algo:\n";
                txt << "  Dijkstra : " << avg(o_dij) << "  (n=" << o_dij.n << ")\n";
                txt << "  A*       : " << avg(o_astar) << "  (n=" << o_astar.n << ")\n";
                txt << "  ALT uni  : " << avg(o_uni) << "  (n=" << o_uni.n << ")\n";
                txt << "  ALT bi   : " << avg(o_bi) << "  (n=" << o_bi.n << ")\n";
                txt << "  ALT OR   : " << avg(o_or) << "  (n=" << o_or.n << ")\n\n";

                auto bucket_name = [&](int b)->const char* {
                    static const char* N[5] = { "XS","S","M","L","XL" };
                    return N[b];
                    };

                txt << "Per-Quantile-Bucket Averages (ms):\n";
                for (int b = 0; b < 5; ++b) {
                    txt << "  [" << bucket_name(b) << "]  Queries: " << B[b].count << "\n";
                    txt << "    Dijkstra : " << avg(B[b].dij) << "  (n=" << B[b].dij.n << ")\n";
                    txt << "    A*       : " << avg(B[b].astar) << "  (n=" << B[b].astar.n << ")\n";
                    txt << "    ALT uni  : " << avg(B[b].uni) << "  (n=" << B[b].uni.n << ")\n";
                    txt << "    ALT bi   : " << avg(B[b].bi) << "  (n=" << B[b].bi.n << ")\n";
                    txt << "    ALT OR   : " << avg(B[b].orr) << "  (n=" << B[b].orr.n << ")\n\n";
                }

                txt << "Aktueller OR_SWITCH_METERS: " << OR_SWITCH_METERS << " m\n";
                txt.close();
                std::cout << "[INFO] TXT geschrieben: routing_summary.txt\n";
            }
        }
    }

    return 0;
}






#include <vector>
#include <unordered_map>
#include <array>
#include <cmath>



// ----- Konfig -----
static constexpr int NUM_QUERIES = 1000;
static const std::string FMI_FILE = "graph_output_fixed.fmi";
static const std::string TABLE_FOLDER = "tables_coast_manual_105LM"; // ggf. anpassen

// ====== Drop-in Implementierungen: dist_proxy, estimate_quantiles, build_queries_stratified ======
#include <random>
#include <numeric>




// ====== Ende Drop-in ======


extern std::vector<std::vector<double>> loadAllLandmarkTables(
    const std::string& folder,
    const std::vector<int>& landmarkIds,
    bool backward);

// ---- LM-IDs aus Ordner lesen (L_fwd_*.tbl / L_rev_*.tbl) ----
static std::vector<int> scanLandmarkIdsFromFolder(const std::string& tableFolder) {
    std::vector<int> landmarkIds;
    for (auto& de : std::filesystem::directory_iterator(tableFolder)) {
        if (!de.is_regular_file()) continue;
        auto fname = de.path().filename().string();
        if (fname.rfind("L_fwd_", 0) == 0 || fname.rfind("L_rev_", 0) == 0) {
            auto us = fname.find_last_of('_');
            auto dot = fname.find_last_of('.');
            if (us != std::string::npos && dot != std::string::npos && dot > us + 1) {
                try { landmarkIds.push_back(std::stoi(fname.substr(us + 1, dot - us - 1))); }
                catch (...) {}
            }
        }
    }
    std::sort(landmarkIds.begin(), landmarkIds.end());
    landmarkIds.erase(std::unique(landmarkIds.begin(), landmarkIds.end()), landmarkIds.end());
    return landmarkIds;
}





