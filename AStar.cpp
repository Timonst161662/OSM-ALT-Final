#include <vector>
#include <queue>
#include <unordered_map>
#include <functional>   // std::greater
#include <cmath>
#include <chrono>
#include <limits>
#include <iostream> 

#include "AStar.h"
#include "GraphUtils.h"
#include "Geo.h"

// kein #define PI, kein routing.hpp

using namespace std;

bool astar(int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous,
    const std::unordered_map<int, Node>& nodeMap)
{
    if (source < 0 || target < 0 ||
        source >= (int)adjList.size() || target >= (int)adjList.size())
        return false;

    const int n = (int)adjList.size();
    distances.assign(n, std::numeric_limits<double>::infinity());
    previous.assign(n, -1);

    auto h = [&](int u)->double {
        const auto& A = nodeMap.at(u);
        const auto& B = nodeMap.at(target);
        return geo::haversine(A.lat, A.lon, B.lat, B.lon);
        };

    using Item = std::pair<double, int>; // (f = g+h, node)
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    distances[source] = 0.0;
    pq.emplace(h(source), source);

    while (!pq.empty()) {
        auto [f_u, u] = pq.top(); pq.pop();
        // stale-Check: echter f* = g[u] + h(u)
        if (f_u > distances[u] + h(u)) continue;

        if (u == target) {
            return true; // früh raus — KEIN Print
        }

        for (const auto& [v, w] : adjList[u]) {
            if (!std::isfinite(w) || w < 0) continue;
            double g_v = distances[u] + w;
            if (g_v < distances[v]) {
                distances[v] = g_v;
                previous[v] = u;
                pq.emplace(g_v + h(v), v);
            }
        }
    }
    return false; // ebenfalls stumm
}

