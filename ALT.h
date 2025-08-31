#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <limits>
#include <queue>
#include <nlohmann/json.hpp>
#include "GraphUtils.h"
#include "Dijkstra.h"

using namespace std;
using json = nlohmann::json;

// === OR-Switch Konfiguration ===
// Schaltet bei Luftlinie > OR_SWITCH_METERS auf ALT-BI, sonst ALT-UNI.
extern double OR_SWITCH_METERS;   // z.B. 300000.0
extern bool   OR_LOG_DECISION;    // z.B. true

extern thread_local std::vector<double> g_bi_dist_b;
extern thread_local std::vector<int>    g_bi_next_b;
extern thread_local std::vector<int>    g_bi_touched_b;

extern thread_local std::vector<double> bi_forwardHeuristic;
extern thread_local std::vector<int>    bi_forwardBestLandmark;
extern thread_local std::vector<int>    bi_forwardTouched;

extern thread_local std::vector<double> bi_backwardHeuristic;
extern thread_local std::vector<int>    bi_backwardBestLandmark;
extern thread_local std::vector<int>    bi_backwardTouched;


#include <atomic>
inline std::atomic<bool> g_alt_abort{ false };
inline void alt_request_abort() { g_alt_abort.store(true, std::memory_order_relaxed); }
inline void alt_clear_abort() { g_alt_abort.store(false, std::memory_order_relaxed); }
inline bool alt_should_abort() { return g_alt_abort.load(std::memory_order_relaxed); } // <— NEU

enum class QueryBucket { XS = 0, S = 1, M = 2, L = 3, XL = 4, Unknown = 99 };


extern bool g_tbl_verbose;  // default: false
void setAltCompactLog(bool on);  // schaltet Einzeiler an/aus

// === Warmup-Pass für Landmark-Tabellen & ALT-Heuristik ===
#include <random>
#include <chrono>

static volatile double g_warmup_sink = 0.0; // verhindert Wegoptimieren

template <typename T>
static void warmupTables(const std::vector<std::vector<T>>& tables, const char* label)
{
    using Clock = std::chrono::high_resolution_clock;
    const auto t0 = Clock::now();

    double acc = 0.0;
    for (const auto& tbl : tables) {
        const T* p = tbl.data();
        const size_t n = tbl.size();
        for (size_t i = 0; i < n; ++i) acc += static_cast<double>(p[i]);
    }

    g_warmup_sink += acc;
    const auto t1 = Clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[INFO] LM-Warmup (" << label << "): " << ms << " ms, sink=" << g_warmup_sink << "\n";
}

// Hilfsfunktion: sammelt bis zu 'limit' Nachbarn aus der Adjazenz für einen Startknoten
static void collectSomeNeighbors(int start,
    const WeightedAdjList& adj,
    std::vector<int>& out, size_t limit = 64)
{
    out.clear();
    auto it = adj.find(start);
    if (it == adj.end()) return;
    for (const auto& [v, w] : it->second) {
        out.push_back(v);
        if (out.size() >= limit) break;
    }
    // falls zu wenige unmittelbare Nachbarn, noch zweite Ebene anstubsen
    if (out.size() < limit) {
        for (int u : out) {
            auto jt = adj.find(u);
            if (jt == adj.end()) continue;
            for (const auto& [v2, w2] : jt->second) {
                if (out.size() >= limit) break;
                out.push_back(v2);
            }
            if (out.size() >= limit) break;
        }
    }
}

// Führt für mehrere (s,t) Paare die ALT-Selektor- und Heuristik-Zugriffe durch (uni + bi)
static void warmupALT_SelectorsAndHeuristic(
    const WeightedAdjList& adjF,
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    const std::vector<std::pair<int, int>>& seedPairs,   // einige (s,t) – z.B. erste 16 Queries
    int topN = 4,
    size_t touchPerSide = 64                            // wie viele Knoten rund um s/t für Heuristik anfassen
) {
    using Clock = std::chrono::high_resolution_clock;
    const auto t0 = Clock::now();

    if (tables_forward.empty()) return;

    const size_t L = std::min(tables_forward.size(), landmarkIds.size());

    std::vector<std::pair<double, int>> diffs;
    diffs.reserve(L);
    std::vector<int> activeLm; activeLm.reserve((size_t)topN);

    std::vector<int> neigh; neigh.reserve(touchPerSide);

    double sink = 0.0;

    for (const auto& st : seedPairs) {
        const int s = st.first;
        const int t = st.second;

        // --- (1) Selektor wie ALT: |F[t] - F[s]| über alle Landmarks, Top-N
        diffs.clear();
        for (size_t i = 0; i < L; ++i) {
            const auto& F = tables_forward[i];
            if (s < 0 || t < 0 || s >= (int)F.size() || t >= (int)F.size()) continue;
            const double fs = F[s], ft = F[t];
            if (!std::isfinite(fs) || !std::isfinite(ft)) continue;
            const double d = std::abs(ft - fs);
            diffs.emplace_back(d, (int)i);
        }
        std::partial_sort(diffs.begin(),
            diffs.begin() + std::min((size_t)topN, diffs.size()),
            diffs.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
        diffs.resize(std::min((size_t)topN, diffs.size()));

        activeLm.clear();
        for (const auto& p : diffs) activeLm.push_back(p.second);

        // --- (2) Heuristikzugriffe (uni): best = max(|F[t]-F[v]|, |B[v]-B[t]|)
        collectSomeNeighbors(s, adjF, neigh, touchPerSide);
        neigh.push_back(s);
        neigh.push_back(t);
        for (int v : neigh) {
            double best = 0.0;
            for (int idx : activeLm) {
                const auto& F = tables_forward[idx];
                const auto& B = (idx < (int)tables_backward.size()) ? tables_backward[idx] : tables_forward[idx];
                if (v >= 0 && t >= 0 &&
                    v < (int)F.size() && t < (int)F.size() &&
                    std::isfinite(F[v]) && std::isfinite(F[t])) {
                    best = std::max(best, std::abs(F[t] - F[v]));
                }
                if (v >= 0 && t >= 0 &&
                    v < (int)B.size() && t < (int)B.size() &&
                    std::isfinite(B[v]) && std::isfinite(B[t])) {
                    best = std::max(best, std::abs(B[v] - B[t]));
                }
            }
            // zusätzlich Haversine antippen (Fallback in alt)
            const auto& nv = nodeMap.at(v);
            const auto& nt = nodeMap.at(t);
            best = std::max(best, geo::haversine(nv.lat, nv.lon, nt.lat, nt.lon));
            sink += best;
        }

        // --- (3) Heuristikzugriffe (bi, rückwärts-Seite): goal=source
        collectSomeNeighbors(t, adjF, neigh, touchPerSide);
        neigh.push_back(t);
        neigh.push_back(s);
        for (int v : neigh) {
            double best = 0.0;
            for (int idx : activeLm) {
                const auto& F = tables_forward[idx];
                const auto& B = (idx < (int)tables_backward.size()) ? tables_backward[idx] : tables_forward[idx];
                if (v >= 0 && s >= 0 &&
                    v < (int)F.size() && s < (int)F.size() &&
                    std::isfinite(F[v]) && std::isfinite(F[s])) {
                    best = std::max(best, std::abs(F[s] - F[v]));
                }
                if (v >= 0 && s >= 0 &&
                    v < (int)B.size() && s < (int)B.size() &&
                    std::isfinite(B[v]) && std::isfinite(B[s])) {
                    best = std::max(best, std::abs(B[v] - B[s]));
                }
            }
            const auto& nv = nodeMap.at(v);
            const auto& ns = nodeMap.at(s);
            best = std::max(best, geo::haversine(nv.lat, nv.lon, ns.lat, ns.lon));
            sink += best;
        }
    }

    g_warmup_sink += sink;
    const auto t1 = Clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[INFO] ALT-Warmup (Selector+Heuristic uni/bi): " << ms << " ms, sink=" << g_warmup_sink << "\n";
}



struct AltLandmarkUsage {
    std::vector<uint32_t> counts;
    std::vector<uint8_t> used_bitset;
    int best_at_source = -1;
    int best_at_target = -1;
    explicit AltLandmarkUsage(size_t K = 0)
        : counts(K, 0), used_bitset((K + 7) / 8, 0) {
    }
    inline void mark_used(size_t l) {
        counts[l] += 1u;
        used_bitset[l >> 3] |= (1u << (l & 7));
    }
    inline bool was_used(size_t l) const {
        return (used_bitset[l >> 3] >> (l & 7)) & 1u;
    }
};

struct PQItem {
    double key;
    int node;
    bool operator>(const PQItem& o) const { return key > o.key; }
};

struct PQMinCmp {                 // Min-Heap
    bool operator()(const PQItem& a, const PQItem& b) const {
        return a.key > b.key;
    }
};

// Rückwärtsgraph direkt aus der kompakten Adjazenzliste bauen
static inline WeightedAdjList buildReverseFromVec(
    const std::vector<std::vector<std::pair<int, double>>>& adjListVec)
{
    WeightedAdjList R;
    R.reserve(adjListVec.size());
    for (int u = 0; u < (int)adjListVec.size(); ++u) {
        for (const auto& [v, w] : adjListVec[u]) {
            R[v].emplace_back(u, w);
        }
    }
    return R;
}

// Heuristik-Ergebnis inkl. Argmax-Landmark (für Debug)
struct HeuRes { double h; int bestL; };

double computeALTHeuristic(int v, int target,
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,
    const std::vector<int>& landmarkIds);

std::vector<std::vector<double>> loadAllLandmarkTables(
    const std::string& folder,
    const std::vector<int>& landmarkIds,
    bool reverse = false);

std::vector<int> rankTopKByPeek(
    const std::string& folder,
    const std::vector<int>& allLmIds,
    int s, int t, int K);

bool alt(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjListVec, // <— statt WeightedAdjList
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    std::vector<double>& distances,
    std::vector<int>& previous);

bool alt_bidirectional(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjListVec,      // forward
    const std::vector<std::vector<std::pair<int, double>>>& adjListVecRev,   // reverse
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,                // bei ungerichtet Alias zu forward ok
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    std::vector<double>& out_distances,   // entlang finalem Pfad (nur Pfadknoten gesetzt)
    std::vector<int>& out_previous);

double haversine_source_target_m(
    int source, int target,
    const std::unordered_map<int, Node>& nodeMap);

inline bool alt_uniOrBi(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjF,   // forward adjListVec
    const std::vector<std::vector<std::pair<int, double>>>& adjR,   // reverse adjListVec (bei ungerichtet = adjF)
    const std::vector<std::vector<double>>& tables_forward,
    const std::vector<std::vector<double>>& tables_backward,
    const std::vector<int>& landmarkIds,
    const std::unordered_map<int, Node>& nodeMap,
    std::vector<double>& distances,
    std::vector<int>& previous,
    QueryBucket bucket = QueryBucket::Unknown,
    bool log_decision = true
) {
    // Luftlinien-Schalter (gleiche Logik wie zuvor)
    const double dd = haversine_source_target_m(source, target, nodeMap);
#ifdef OR_LOG_DECISION
    if (log_decision) {
        std::cout << "[ALT-OR] great-circle=" << std::llround(dd)
            << " m | thr=" << std::llround(OR_SWITCH_METERS)
            << " → mode=" << (dd > OR_SWITCH_METERS ? "BI" : "UNI") << "\n";
    }
#else
    (void)log_decision; // falls Macro nicht gesetzt ist
#endif

    if (dd > OR_SWITCH_METERS) {
        // Weit: bidirektionales ALT auf Vektor-Adjazenz
        return alt_bidirectional(
            source, target,
            adjF, adjR,
            tables_forward, tables_backward,
            landmarkIds, nodeMap,
            distances, previous
        );
    }

    // Nah: unidirektionales ALT auf Vektor-Adjazenz
    return alt(
        source, target,
        adjF,
        tables_forward, tables_backward,
        landmarkIds, nodeMap,
        distances, previous
    );
}


// --- Params, die du ggf. oben zentralisieren willst ---
static constexpr double TIGHT_SYMMETRY_THR = 0.20; // <= 20% Unterschied -> "balanced" -> BI
static constexpr double MIN_TIGHT_FRAC_H = 0.35; // BI nur, wenn beide Tightness ~>= 35% der Luftlinie
// ------------------------------------------------------



