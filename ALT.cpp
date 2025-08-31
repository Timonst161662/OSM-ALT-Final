#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <queue>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "ALT.h"
#include "GraphUtils.h"
#include "Geo.h"

using namespace std;
using json = nlohmann::json;

bool g_tbl_verbose = false;
static bool g_alt_compact = true;
void setAltCompactLog(bool on) { g_alt_compact = on; }

double OR_SWITCH_METERS = 674087.637; //Entscheidungsschwelle für uni- or bidirectional
bool   OR_LOG_DECISION = true;     // Konsolenlog für Benchmarks

//ALT-Bidirectional: Heuristik-Caches

thread_local std::vector<double> g_bi_dist_b;
thread_local std::vector<int>    g_bi_next_b;
thread_local std::vector<int>    g_bi_touched_b;

thread_local std::vector<double> bi_forwardHeuristic;
thread_local std::vector<int>    bi_forwardBestLandmark;
thread_local std::vector<int>    bi_forwardTouched;

thread_local std::vector<double> bi_backwardHeuristic;
thread_local std::vector<int>    bi_backwardBestLandmark;
thread_local std::vector<int>    bi_backwardTouched;


// Bindings aus Geo:
using geo::haversine; 

//Luftlinie in Metern zwischen zwei Node-IDs
double haversine_source_target_m(
    int source, int target,
    const std::unordered_map<int, Node>& nodeMap)
{
    const auto& A = nodeMap.at(source);
    const auto& B = nodeMap.at(target);
    
    return haversine(A.lat, A.lon, B.lat, B.lon);
}

std::vector<std::vector<double>> loadAllLandmarkTables(
    const std::string& folder,
    const std::vector<int>& landmarkIds,
    bool reverse)
{
    std::vector<std::vector<double>> tables;
    tables.reserve(landmarkIds.size());

    size_t expected_len = 0;
    size_t okCount = 0, failCount = 0;
    std::unordered_map<size_t, int> formatCount;

    for (int lid : landmarkIds) {
        std::vector<double> table;

        if (reverse) {
            std::ostringstream p1; p1 << folder << "/L_rev_" << lid << ".tbl";
            table = loadDistancesFromBinary(p1.str());
            if (table.empty()) {
                std::ostringstream p2; p2 << folder << "/L_" << lid << ".tbl";
                table = loadDistancesFromBinary(p2.str());
            }
        }
        else {
            std::ostringstream p1; p1 << folder << "/L_fwd_" << lid << ".tbl";
            table = loadDistancesFromBinary(p1.str());
            if (table.empty()) {
                std::ostringstream p2; p2 << folder << "/L_" << lid << ".tbl";
                table = loadDistancesFromBinary(p2.str());
            }
        }

        if (!table.empty() && expected_len == 0) {
            expected_len = table.size();
        }

        if (table.empty()) {
            std::cerr << "[WARN] Tabelle fehlt/leer (LM " << lid
                << ", " << (reverse ? "rev" : "fwd") << ")\n";
            if (expected_len == 0) {
                tables.push_back({});
                continue;
            }
            table.assign(expected_len, std::numeric_limits<double>::infinity());
            failCount++;
        }
        else {
            okCount++;
            formatCount[table.size()]++;
        }

        tables.push_back(std::move(table));
    }

    std::cout << "[INFO] LM-Tables geladen: ok=" << okCount
        << ", fail=" << failCount
        << " | Formate: "
        << [&] {
        std::ostringstream oss;
        for (auto& kv : formatCount)
            oss << kv.first << "x" << kv.second << " ";
        return oss.str();
        }() << "\n";

    return tables;
}

#ifndef ALT_TOPN
#define ALT_TOPN 4
#endif

//ALT (unidirektional) mit Vektor-Adjazenz (adjListVec)
bool alt(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjListVec, // <— statt WeightedAdjList
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    std::vector<double>& distances,
    std::vector<int>& previous)
{
    //Guards
    if (!nodeMap.count(source) || !nodeMap.count(target)) return false;
    if (source == target) {
        distances = g_alt_distances;
        previous = g_alt_previous;
        return true;
    }

    //Speicher vorbereiten / Reset
    {
        const size_t need = std::max(static_cast<size_t>(std::max(source, target)) + 1,
            g_alt_distances.size());
        if (g_alt_distances.size() < need) {
            g_alt_distances.resize(need, std::numeric_limits<double>::infinity());
            g_alt_previous.resize(need, -1);
        }
    }
    resetAltTouched();
    g_alt_distances[source] = 0.0;
    g_alt_previous[source] = -1;
    g_alt_touched.push_back(source);

    //Aktive Landmarks wählen (Top-N via forward-Diff)
    std::vector<int> activeLandmarks;
    {
        std::vector<std::pair<double, int>> diffs;
        diffs.reserve(landmarkIds.size());
        for (size_t i = 0; i < landmarkIds.size(); ++i) {
            if (i >= tables_forward.size()) continue;
            const auto& tf = tables_forward[i];
            if (source >= (int)tf.size() || target >= (int)tf.size()) continue;
            const double fs = tf[source], ft = tf[target];
            if (!std::isfinite(fs) || !std::isfinite(ft)) continue;
            double d = std::abs(fs - ft);
            if (!std::isfinite(d)) d = 0.0;
            diffs.emplace_back(d, (int)i);
        }
        std::sort(diffs.rbegin(), diffs.rend());
        const int topN = std::min(ALT_TOPN, (int)diffs.size());
        activeLandmarks.reserve(topN);
        for (int j = 0; j < topN; ++j) activeLandmarks.push_back(diffs[j].second);
    }
    //Logging: Vorfilter (aktive LMs)
    /*std::cout << "[ALT] activeLM.size=" << activeLandmarks.size()
        << " (of " << landmarkIds.size() << ")\n";
    std::cout << "[ALT] activeLM IDs: ";
    for (size_t j = 0; j < activeLandmarks.size(); ++j) {
        int idx = activeLandmarks[j];
        if (idx >= 0 && idx < (int)landmarkIds.size()) {
            std::cout << landmarkIds[idx];
            if (j + 1 < activeLandmarks.size()) std::cout << ", ";
        }
    }
    std::cout << "\n";*/

    //Heuristik-Cache (vektorbasiert)
    const size_t N = g_alt_distances.size();
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> hcache_val(N, NaN);
    std::vector<int>    hcache_best(N, -1);

    auto haversine_safe = [&](int a, int b) {
        const auto& A = nodeMap.at(a);
        const auto& B = nodeMap.at(b);
        return haversine(A.lat, A.lon, B.lat, B.lon);
        };

    //verwendete LMs tracken (für UI)
    g_alt_lastUsedLandmarks.clear();
    std::vector<uint8_t> lmUsedIdx(landmarkIds.size(), 0);

    auto computeH = [&](int v) -> std::pair<double, int> {
        if ((size_t)v < hcache_val.size()) {
            double hv = hcache_val[(size_t)v];
            if (hv == hv) return { hv, hcache_best[(size_t)v] };
        }

        double max_diff = 0.0;
        int    best = -1;
        int    valid = 0;

        for (int i : activeLandmarks) {
            if (i < 0 || i >= (int)tables_forward.size() || i >= (int)tables_backward.size())
                continue;
            const auto& tf = tables_forward[i];
            const auto& tb = tables_backward[i];
            if (v >= (int)tf.size() || target >= (int)tf.size()) continue;
            if (v >= (int)tb.size() || target >= (int)tb.size()) continue;

            const double f_t = tf[target], f_v = tf[v];
            const double b_v = tb[v], b_t = tb[target];
            if (!std::isfinite(f_t) || !std::isfinite(f_v) ||
                !std::isfinite(b_v) || !std::isfinite(b_t)) continue;

            double d1 = std::abs(f_t - f_v); // |d(L,target) - d(L,v)|
            double d2 = std::abs(b_v - b_t); // |d(v,L)      - d(target,L)|
            if (!std::isfinite(d1)) d1 = 0.0;
            if (!std::isfinite(d2)) d2 = 0.0;

            const double val = (d1 > d2) ? d1 : d2;
            if (val > max_diff) { max_diff = val; best = i; }
            ++valid;
        }

        if (valid == 0 || max_diff <= 0.0) {
            max_diff = haversine_safe(v, target); // Fallback
        }
        if (best >= 0 && (size_t)best < lmUsedIdx.size()) lmUsedIdx[(size_t)best] = 1;

        if ((size_t)v < hcache_val.size()) {
            hcache_val[(size_t)v] = max_diff;
            hcache_best[(size_t)v] = best;
        }
        return { max_diff, best };
        };

    //PQ (f = g + h)
    using PQItem = std::pair<double, int>;
    std::priority_queue<PQItem, std::vector<PQItem>, std::greater<>> pq;

    { // Start
        auto [h0, _] = computeH(source);
        pq.emplace(h0, source); // g=0 -> f=h(source)
    }

    //Loop
    while (!pq.empty()) {
        const auto [fest, u] = pq.top();
        pq.pop();

        const auto [hu, _b] = computeH(u);
        if (fest > g_alt_distances[u] + hu) continue; // veralteter PQ-Eintrag

        if (u == target) {
            distances = g_alt_distances;
            previous = g_alt_previous;

            g_alt_lastUsedLandmarks.clear();
            for (size_t i = 0; i < lmUsedIdx.size(); ++i)
                if (lmUsedIdx[i] && i < landmarkIds.size())
                    g_alt_lastUsedLandmarks.push_back(landmarkIds[(int)i]);
            return true;
        }

        if ((size_t)u >= adjListVec.size()) continue; // Schutz (sollte nicht auftreten)
        for (const auto& [v, w] : adjListVec[(size_t)u]) {
            if (!std::isfinite(w) || w < 0) continue;
            if ((size_t)v >= g_alt_distances.size()) continue; // sollte bei initAltGlobals nicht passieren

            const double nd = g_alt_distances[u] + w;
            if (!std::isfinite(nd)) continue;

            if (nd < g_alt_distances[v]) {
                if (!std::isfinite(g_alt_distances[v])) g_alt_touched.push_back(v); // first touch
                g_alt_distances[v] = nd;
                g_alt_previous[v] = u;

                const auto [hv, _b2] = computeH(v);
                pq.emplace(nd + hv, v);
            }
        }
    }

    //Kein Pfad
    distances = g_alt_distances;
    previous = g_alt_previous;
    g_alt_lastUsedLandmarks.clear();
    for (size_t i = 0; i < lmUsedIdx.size(); ++i)
        if (lmUsedIdx[i] && i < landmarkIds.size())
            g_alt_lastUsedLandmarks.push_back(landmarkIds[(int)i]);
    //Logging: tatsächlich verwendete LMs
    /*std::cout << "[ALT] usedLMs.count=" << g_alt_lastUsedLandmarks.size() << " → IDs: ";
    for (size_t j = 0; j < g_alt_lastUsedLandmarks.size(); ++j) {
        std::cout << g_alt_lastUsedLandmarks[j];
        if (j + 1 < g_alt_lastUsedLandmarks.size()) std::cout << ", ";
    }
    std::cout << "\n";*/

    return false;
}

//alt_bidirectional
bool alt_bidirectional(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjListVec,     
    const std::vector<std::vector<std::pair<int, double>>>& adjListVecRev,   
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,                 
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    std::vector<double>& out_distances,
    std::vector<int>& out_previous)
{
    //Guards
    if (!nodeMap.count(source) || !nodeMap.count(target)) return false;
    if (source == target) { out_distances = g_alt_distances; out_previous = g_alt_previous; return true; }

    //g_alt_* vorbereiten
    {
        const size_t need = std::max(static_cast<size_t>(std::max(source, target)) + 1, g_alt_distances.size());
        if (g_alt_distances.size() < need) {
            g_alt_distances.resize(need, std::numeric_limits<double>::infinity());
            g_alt_previous.resize(need, -1);
        }
    }
    resetAltTouched();

    //einmalig genug Platz für viele Touches geben
    if (g_alt_touched.capacity() < (1u << 16))  // ~65k
        g_alt_touched.reserve(1u << 16);

    g_alt_distances[(size_t)source] = 0.0;
    g_alt_previous[(size_t)source] = -1;
    g_alt_touched.push_back(source);

    const double INF = std::numeric_limits<double>::infinity();
    const size_t N = g_alt_distances.size();

    //Priority-Queues (f = g + h)
    using PQItem = std::pair<double, int>;  // <fScore, nodeId>

    //Startkapazitäten
    std::vector<PQItem> bufF; bufF.reserve(32768);
    std::vector<PQItem> bufB; bufB.reserve(32768);

   
    std::priority_queue<PQItem, std::vector<PQItem>, std::greater<>>
        pq_f(std::greater<>{}, std::move(bufF)),
        pq_b(std::greater<>{}, std::move(bufB));

    auto fmin = [&](const auto& pq)->double { return pq.empty() ? INF : pq.top().first; };

    //Aktive Landmarks (TopN nach |tf[source]-tf[target]|)
    std::vector<int> activeLandmarks;
    {
        std::vector<std::pair<double, int>> diffs; diffs.reserve(landmarkIds.size());
        for (size_t i = 0; i < landmarkIds.size(); ++i) {
            if (i >= tables_forward.size()) continue;
            const auto& tf = tables_forward[i];
            if (source >= (int)tf.size() || target >= (int)tf.size()) continue;
            const double fs = tf[source], ft = tf[target];
            if (!std::isfinite(fs) || !std::isfinite(ft)) continue;
            double d = std::abs(fs - ft); if (!std::isfinite(d)) d = 0.0;
            diffs.emplace_back(d, (int)i);
        }
        std::sort(diffs.rbegin(), diffs.rend());
        const int topN = std::min(ALT_TOPN, (int)diffs.size());
        activeLandmarks.reserve(topN);
        for (int j = 0; j < topN; ++j) activeLandmarks.push_back(diffs[j].second);
    }

    //Heuristik-Caches
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    if (bi_forwardHeuristic.size() < N) bi_forwardHeuristic.resize(N, NaN);
    if (bi_forwardBestLandmark.size() < N) bi_forwardBestLandmark.resize(N, -1);
    if (bi_backwardHeuristic.size() < N) bi_backwardHeuristic.resize(N, NaN);
    if (bi_backwardBestLandmark.size() < N) bi_backwardBestLandmark.resize(N, -1);

    if (!bi_forwardTouched.empty()) {
        for (int v : bi_forwardTouched) if ((size_t)v < N) { bi_forwardHeuristic[(size_t)v] = NaN; bi_forwardBestLandmark[(size_t)v] = -1; }
        bi_forwardTouched.clear();
    }
    if (!bi_backwardTouched.empty()) {
        for (int v : bi_backwardTouched) if ((size_t)v < N) { bi_backwardHeuristic[(size_t)v] = NaN; bi_backwardBestLandmark[(size_t)v] = -1; }
        bi_backwardTouched.clear();
    }

    std::vector<uint8_t> lmUsedIdx(landmarkIds.size(), 0);

    auto haversine_safe = [&](int a, int b) {
        const auto& A = nodeMap.at(a); const auto& B = nodeMap.at(b);
        return haversine(A.lat, A.lon, B.lat, B.lon);
        };

    auto computeH_fwd = [&](int v) -> double {
        if ((size_t)v < N) {
            double hv = bi_forwardHeuristic[(size_t)v];
            if (hv == hv) return hv; // Cache-Hit
        }
        double bestVal = 0.0; int bestLm = -1; int valid = 0;
        for (int i : activeLandmarks) {
            if (i < 0 || i >= (int)tables_forward.size() || i >= (int)tables_backward.size()) continue;
            const auto& tf = tables_forward[i];
            const auto& tb = tables_backward[i];
            if (v < 0 || target < 0 || v >= (int)tf.size() || target >= (int)tf.size()) continue;
            if (v >= (int)tb.size() || target >= (int)tb.size()) continue;

            double d1 = std::abs(tf[target] - tf[v]);      if (!std::isfinite(d1)) d1 = 0.0;
            double d2 = std::abs(tb[v] - tb[target]); if (!std::isfinite(d2)) d2 = 0.0;
            double val = (d1 > d2) ? d1 : d2;

            if (val > bestVal) { bestVal = val; bestLm = i; }
            ++valid;
        }
        if (valid == 0 || bestVal <= 0.0) bestVal = haversine_safe(v, target);
        if (bestLm >= 0 && (size_t)bestLm < lmUsedIdx.size()) lmUsedIdx[(size_t)bestLm] = 1;

        if ((size_t)v < N) {
            bi_forwardHeuristic[(size_t)v] = bestVal;
            bi_forwardBestLandmark[(size_t)v] = bestLm;
            bi_forwardTouched.push_back(v);
        }
        return bestVal;
        };

    auto computeH_bwd = [&](int v) -> double {
        if ((size_t)v < N) {
            double hv = bi_backwardHeuristic[(size_t)v];
            if (hv == hv) return hv; // Cache-Hit
        }
        double bestVal = 0.0; int bestLm = -1; int valid = 0;
        for (int i : activeLandmarks) {
            if (i < 0 || i >= (int)tables_forward.size() || i >= (int)tables_backward.size()) continue;
            const auto& tf = tables_forward[i];
            const auto& tb = tables_backward[i];
            if (v < 0 || source < 0 || v >= (int)tf.size() || source >= (int)tf.size()) continue;
            if (v >= (int)tb.size() || source >= (int)tb.size()) continue;

            double d1 = std::abs(tf[source] - tf[v]);      if (!std::isfinite(d1)) d1 = 0.0;
            double d2 = std::abs(tb[v] - tb[source]); if (!std::isfinite(d2)) d2 = 0.0;
            double val = (d1 > d2) ? d1 : d2;

            if (val > bestVal) { bestVal = val; bestLm = i; }
            ++valid;
        }
        if (valid == 0 || bestVal <= 0.0) bestVal = haversine_safe(v, source);
        if (bestLm >= 0 && (size_t)bestLm < lmUsedIdx.size()) lmUsedIdx[(size_t)bestLm] = 1;

        if ((size_t)v < N) {
            bi_backwardHeuristic[(size_t)v] = bestVal;
            bi_backwardBestLandmark[(size_t)v] = bestLm;
            bi_backwardTouched.push_back(v);
        }
        return bestVal;
        };

    //Backward-Zustände
    if (g_bi_dist_b.size() < N) {
        g_bi_dist_b.resize(N, INF);
        g_bi_next_b.resize(N, -1);
    }
    else {
        for (int i : g_bi_touched_b) {
            if ((size_t)i < g_bi_dist_b.size()) { g_bi_dist_b[(size_t)i] = INF; g_bi_next_b[(size_t)i] = -1; }
        }
        g_bi_touched_b.clear();
    }

    //Startzustände
    {
        const double h0f = computeH_fwd(source);
        pq_f.emplace(g_alt_distances[(size_t)source] + h0f, source);

        g_bi_dist_b[(size_t)target] = 0.0;
        g_bi_next_b[(size_t)target] = -1;
        g_bi_touched_b.push_back(target);

        const double h0b = computeH_bwd(target);
        pq_b.emplace(g_bi_dist_b[(size_t)target] + h0b, target);
    }

    auto edge_w = [&](int u, int v)->double {
        if (u < 0 || (size_t)u >= adjListVec.size()) return INF;
        for (const auto& e : adjListVec[(size_t)u]) if (e.first == v) return e.second;
        return INF;
        };

    double mu = INF;
    int    meet = -1;

    //Vorwärts-Expansion
    auto relax_forward = [&]() -> bool {
        while (!pq_f.empty()) {
            auto [fest, u] = pq_f.top(); pq_f.pop();
            if ((size_t)u >= N) continue;
            const double hu = computeH_fwd(u);
            if (fest > g_alt_distances[(size_t)u] + hu) continue;

            // meet auf u
            if (std::isfinite(g_bi_dist_b[(size_t)u])) {
                const double cand = g_alt_distances[(size_t)u] + g_bi_dist_b[(size_t)u];
                if (cand < mu) { mu = cand; meet = u; }
            }

            if ((size_t)u >= adjListVec.size()) return true;
            for (const auto& [v, w] : adjListVec[(size_t)u]) {
                if (!std::isfinite(w) || w < 0) continue;
                if ((size_t)v >= N) continue;

                const double nd = g_alt_distances[(size_t)u] + w;
                if (!std::isfinite(nd) || nd >= mu) continue;

                if (nd < g_alt_distances[(size_t)v]) {
                    if (!std::isfinite(g_alt_distances[(size_t)v])) g_alt_touched.push_back(v);
                    g_alt_distances[(size_t)v] = nd;
                    g_alt_previous[(size_t)v] = u;

                    // meet auf v
                    if (std::isfinite(g_bi_dist_b[(size_t)v])) {
                        const double cand = nd + g_bi_dist_b[(size_t)v];
                        if (cand < mu) { mu = cand; meet = v; }
                        if (nd >= mu) continue;
                    }

                    const double hv = computeH_fwd(v);
                    pq_f.emplace(nd + hv, v);
                }
            }
            return true;
        }
        return false;
        };

    //Rückwärts-Expansion
    auto relax_backward = [&]() -> bool {
        while (!pq_b.empty()) {
            auto [fest, u] = pq_b.top(); pq_b.pop();
            if ((size_t)u >= N) continue;
            const double hu = computeH_bwd(u);
            if (fest > g_bi_dist_b[(size_t)u] + hu) continue;

            // meet auf u
            if (std::isfinite(g_alt_distances[(size_t)u])) {
                const double cand = g_alt_distances[(size_t)u] + g_bi_dist_b[(size_t)u];
                if (cand < mu) { mu = cand; meet = u; }
            }

            if ((size_t)u >= adjListVecRev.size()) return true;
            for (const auto& [v, w] : adjListVecRev[(size_t)u]) {
                if (!std::isfinite(w) || w < 0) continue;
                if ((size_t)v >= N) continue;

                const double nd = g_bi_dist_b[(size_t)u] + w;
                if (!std::isfinite(nd) || nd >= mu) continue;

                if (nd < g_bi_dist_b[(size_t)v]) {
                    if (!std::isfinite(g_bi_dist_b[(size_t)v])) g_bi_touched_b.push_back(v);
                    g_bi_dist_b[(size_t)v] = nd;
                    g_bi_next_b[(size_t)v] = u; // Vorwärts Richtung target

                    //meet auf v
                    if (std::isfinite(g_alt_distances[(size_t)v])) {
                        const double cand = g_alt_distances[(size_t)v] + nd;
                        if (cand < mu) { mu = cand; meet = v; }
                        if (nd >= mu) continue;
                    }

                    const double hv = computeH_bwd(v);
                    pq_b.emplace(nd + hv, v);
                }
            }
            return true;
        }
        return false;
        };

    //Hauptschleife
    while (!pq_f.empty() || !pq_b.empty()) {
        double fF = fmin(pq_f);
        double fB = fmin(pq_b);
        if (std::min(fF, fB) >= mu) break;

        bool do_fwd = (fF <= fB); // Gleichstand -> forward
        if (do_fwd) {
            if (!relax_forward()) { if (!relax_backward()) break; }
        }
        else {
            if (!relax_backward()) { if (!relax_forward()) break; }
        }
    }

    if (meet == -1 || !std::isfinite(mu)) {
        // verwendete LMs
        g_alt_lastUsedLandmarks.clear();
        for (size_t i = 0; i < lmUsedIdx.size(); ++i)
            if (lmUsedIdx[i] && i < landmarkIds.size())
                g_alt_lastUsedLandmarks.push_back(landmarkIds[(int)i]);
        out_distances = g_alt_distances; out_previous = g_alt_previous;
        return false;
    }

    //Pfad rekonstruieren
    std::vector<int> left;
    for (int a = meet; a != -1; a = (a >= 0 && (size_t)a < N ? g_alt_previous[(size_t)a] : -1)) left.push_back(a);
    std::reverse(left.begin(), left.end());

    std::vector<int> right;
    for (int a = (meet >= 0 && (size_t)meet < N ? g_bi_next_b[(size_t)meet] : -1);
        a != -1; a = (a >= 0 && (size_t)a < N ? g_bi_next_b[(size_t)a] : -1))
    {
        right.push_back(a);
        if (a == target) break;
    }

    if (left.empty() || left.front() != source || right.empty() || right.back() != target) {
        out_distances = g_alt_distances; out_previous = g_alt_previous;
        return false;
    }

    std::vector<int> full = left;
    full.insert(full.end(), right.begin(), right.end());

    //g_alt_* entlang des finalen Pfads (rechte Hälfte nachtragen)
    double acc = 0.0;
    g_alt_distances[(size_t)source] = 0.0;
    g_alt_previous[(size_t)source] = -1;

    auto edge_w_safe = [&](int u, int v)->double {
        if (u < 0 || (size_t)u >= adjListVec.size()) return INF;
        for (const auto& e : adjListVec[(size_t)u]) if (e.first == v) return e.second;
        return INF;
        };

    for (size_t i = 1; i < full.size(); ++i) {
        const int u = full[i - 1], v = full[i];
        const double w = edge_w_safe(u, v);
        if (!std::isfinite(w)) { out_distances = g_alt_distances; out_previous = g_alt_previous; return false; }
        acc += w;

        if ((size_t)v >= g_alt_distances.size()) {
            g_alt_distances.resize((size_t)v + 1, INF);
            g_alt_previous.resize((size_t)v + 1, -1);
        }
        if (!std::isfinite(g_alt_distances[(size_t)v])) g_alt_touched.push_back(v);
        g_alt_distances[(size_t)v] = acc;
        g_alt_previous[(size_t)v] = u;
    }

    //verwendete LMs publizieren
    g_alt_lastUsedLandmarks.clear();
    for (size_t i = 0; i < lmUsedIdx.size(); ++i)
        if (lmUsedIdx[i] && i < landmarkIds.size())
            g_alt_lastUsedLandmarks.push_back(landmarkIds[(int)i]);

    out_distances = g_alt_distances;
    out_previous = g_alt_previous;
    return true;
}