#pragma once

#include <vector>
#include <utility>

// Dijkstra-Algorithmus mit Early-Exit
bool dijkstra_early_exit(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous);

// Pfadrekonstruktion aus dem Ergebnis von Dijkstra
std::vector<int> reconstruct_path(
    int source, int target,
    const std::vector<int>& previous);

// Dijkstra von einem Punkt zu allen anderen (für Landmark-Tabellen)
std::vector<double> dijkstraAll(
    int source,
    const std::vector<std::vector<std::pair<int, double>>>& adjList);
