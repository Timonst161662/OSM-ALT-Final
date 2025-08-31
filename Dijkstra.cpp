#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <limits>
#include <utility>
#include <iterator>
#include <algorithm>
#include <chrono>
#include "GraphUtils.h"


using namespace std;

int num_nodes = 0;
int num_edges = 0;

//Dijkstra logic with priority_queue
// Dijkstra (early exit) — still, kein Print, keine interne Zeitmessung
bool dijkstra_early_exit(int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous)
{
    if (source < 0 || target < 0 ||
        source >= (int)adjList.size() || target >= (int)adjList.size())
        return false;

    const int n = (int)adjList.size();
    distances.assign(n, std::numeric_limits<double>::infinity());
    previous.assign(n, -1);

    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    distances[source] = 0.0;
    pq.emplace(0.0, source);

    while (!pq.empty()) {
        auto [du, u] = pq.top(); pq.pop();
        if (du > distances[u]) continue;
        if (u == target) return true;  // früh raus, KEIN Print

        for (const auto& [v, w] : adjList[u]) {
            if (w < 0 || !std::isfinite(w)) continue;
            double nd = du + w;
            if (nd < distances[v]) {
                distances[v] = nd;
                previous[v] = u;
                pq.emplace(nd, v);
            }
        }
    }
    return false;
}



vector<double> dijkstraAll(int source, const vector<vector<pair<int, double>>>& adjList) {
    int n = adjList.size();
    vector<double> distances(n, numeric_limits<double>::infinity());
    vector<int> previous(n, -1); // optional, falls du Pfade willst

    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq;

    distances[source] = 0.0;
    pq.emplace(0.0, source);

    while (!pq.empty()) {
        auto [dist_u, u] = pq.top();
        pq.pop();

        if (dist_u > distances[u]) continue;

        for (const auto& [v, weight] : adjList[u]) {
            double newDist = distances[u] + weight;
            if (newDist < distances[v]) {
                distances[v] = newDist;
                previous[v] = u; // nur wenn Pfade gebraucht
                pq.emplace(newDist, v);
            }
        }
    }

    return distances;
}


//Stores all nodes which are along the shortest path
vector<int> reconstruct_path(int source, int target, const vector<int>& previous) {
    vector<int> path;
    for (int at = target; at != -1; at = previous[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());
    if (path[0] == source) {
        return path;
    }
    return {};
}

/*int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Wrong Usage! Necessary: " << argv[0] << " <input_file> <source_node_id> <target_node_id>" << endl;
        return 1;
    }
    string filePath = argv[1];
    int source = stoi(argv[2]);
    int target = stoi(argv[3]);

    auto [nodes, edges] = parseFMIFile(filePath);
    if (nodes.empty() || edges.empty()) {
        cerr << "Error reading graph data from file: No Nodes or Edges found." << endl;
    }

    vector<double> distances;
    vector<int> previous;

    vector<vector<pair<int, double>>> adjList(num_nodes + 1);
    for (const auto& e : edges) {
        adjList[e.source_id].emplace_back(e.target_id, e.distance);
        adjList[e.target_id].emplace_back(e.source_id, e.distance);
    }

    bool found = dijkstra_early_exit(source, target, adjList, distances, previous);

    if (found) {
        cout << "Shortest distance from " << source << " to " << target << ": " << distances[target] << endl;
        auto path = reconstruct_path(source, target, previous);
        cout << "Path: ";
        for (int node : path) {
            cout << node << " ";
        }
        cout << endl;
    }
    else {
        cout << "No path from " << source << " to " << target << " exists." << endl;
    }

    return 0;
}*/