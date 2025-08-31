#pragma once
#include <vector>
#include <utility>
#include <unordered_map>

// Vorwärtsdeklaration reicht für den Funktions-Prototyp
struct Node;

// Deklaration der A*-Funktion
bool astar(
    int source, int target,
    const std::vector<std::vector<std::pair<int, double>>>& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous,
    const std::unordered_map<int, Node>& nodeMap);
