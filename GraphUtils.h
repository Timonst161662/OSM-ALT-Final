#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Geo.h" // geo::Point, geo::haversine, geo::PI, ...

// ===== Externe (global) ALT-Strukturen =====
extern std::vector<double> g_alt_distances;
extern std::vector<int>    g_alt_previous;
extern std::vector<int>    g_alt_touched;
extern thread_local std::vector<int> g_alt_lastUsedLandmarks;

// ===== Graph-Grundtypen =====
struct Node { int id; double lat; double lon; };
struct Edge { int a; int b; double length; };

// Adjazenz-Varianten
using AdjacencyList = std::unordered_map<int, std::vector<int>>;
using WeightedAdjList = std::unordered_map<int, std::vector<std::pair<int, double>>>;

// ===== Graph I/O =====
bool loadGraph(const std::string& filename,
    std::vector<Node>& nodes,
    AdjacencyList& adj);

bool loadGraphWeighted(const std::string& filename,
    std::vector<Node>& nodes,
    WeightedAdjList& adj);

std::vector<std::pair<int, int>> loadQueries(const std::string& filename);


void generateRandomQueries(const std::string& filename,
    const std::vector<Node>& nodes,
    int numQueries);

// FMI-Kram
std::pair<std::vector<Node>, std::vector<Edge>>
convertGraphToFmiFormat(const std::vector<geo::Point>& points,
    const std::unordered_map<int, std::unordered_set<int>>& graph);

void saveToFmi(const std::string& filename,
    const std::vector<Node>& nodes,
    const std::vector<Edge>& edges);

std::pair<
    std::unordered_map<int, Node>,
    std::vector<std::vector<std::pair<int, double>>>>
    parseFMIFile(const std::string& filePath);

// ===== Landmark-Tabellen (I/O) =====
void saveDistancesToBinary(const std::string& filename,
    const std::vector<double>& distances);

std::vector<double> loadDistancesFromBinary(const std::string& filename);
bool saveDistancesToBinaryF32(const std::string& path,
    const std::vector<double>& dist);

// ===== ALT / Search-Strukturen =====
void initAltGlobals(std::size_t n);
void resetAltTouched();
void resetAltNode(int v);

// ===== Such-Initialisierung =====
void initializeSearchStructures(int source,
    int target,
    const WeightedAdjList& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous);

// ===== Pfad-/CSV-Utilities (praktisch & klein) =====
std::vector<int>  rebuildPath(int source, int target, const std::vector<int>& prev);
std::string       pathIdsCsv(const std::vector<int>& path);
std::string       pathCoordsCsv(const std::vector<int>& path,
    const std::unordered_map<int, Node>& nodeMap);

void append_csv_row(const std::string& filepath,
    int qidx,
    const char* algo,
    int source, int target,
    bool found,
    long long time_ms,
    double dist_target,
    size_t path_nodes,
    double source_lat, double source_lon,
    double target_lat, double target_lon);
