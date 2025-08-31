#include "GraphUtils.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

using namespace std;

// ===== Globale ALT-Strukturen =====
std::vector<double> g_alt_distances;
std::vector<int>    g_alt_previous;
std::vector<int>    g_alt_touched;
thread_local std::vector<int> g_alt_lastUsedLandmarks;


// ===== Graph I/O (FMI-ähnlich) =====
bool loadGraph(const std::string& filename,
    std::vector<Node>& nodes,
    AdjacencyList& adj)
{
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "[ERROR] Kann Datei nicht öffnen: " << filename << "\n";
        return false;
    }

    int numNodes, numEdges;
    if (!(in >> numNodes >> numEdges)) {
        std::cerr << "[ERROR] Fehler beim Lesen von Kopfzeile\n";
        return false;
    }
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    nodes.resize(numNodes);

    std::string line;
    for (int i = 0; i < numNodes; ++i) {
        if (!std::getline(in, line)) {
            std::cerr << "[ERROR] Datei zu kurz: Erwartet " << numNodes << " Knotenzeilen\n";
            return false;
        }
        std::istringstream iss(line);
        int id; double lat, lon;
        if (!(iss >> id >> lat >> lon)) {
            std::cerr << "[ERROR] Ungültige Knotenzeile bei i=" << i << ": " << line << "\n";
            return false;
        }
        if (id < 0 || id >= numNodes) {
            std::cerr << "[ERROR] Node-ID außerhalb des erlaubten Bereichs: " << id << "\n";
            return false;
        }
        nodes[id] = Node{ id, lat, lon };
    }

    for (int i = 0; i < numEdges; ++i) {
        if (!std::getline(in, line)) {
            std::cerr << "[ERROR] Datei zu kurz: Erwartet " << numEdges << " Kantenzeilen\n";
            return false;
        }
        std::istringstream iss(line);
        int src, tgt; double dist;
        if (!(iss >> src >> tgt >> dist)) {
            std::cerr << "[ERROR] Ungültige Kantenzeile bei i=" << i << ": " << line << "\n";
            return false;
        }
        adj[src].push_back(tgt);
    }
    return true;
}

bool loadGraphWeighted(const std::string& filename,
    std::vector<Node>& nodes,
    WeightedAdjList& adj)
{
    std::ifstream in(filename);
    if (!in) return false;

    int numNodes, numEdges;
    in >> numNodes >> numEdges;

    nodes.resize(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        int id; double lat, lon;
        in >> id >> lat >> lon;
        nodes[i] = { id, lat, lon };
    }

    for (int i = 0; i < numEdges; ++i) {
        int src, tgt; double dist;
        in >> src >> tgt >> dist;
        adj[src].emplace_back(tgt, dist);
        // falls ungerichtet: adj[tgt].emplace_back(src, dist);
    }
    return true;
}

// ===== Queries / CSV =====
std::vector<std::pair<int, int>> loadQueries(const std::string& filename) {
    std::vector<std::pair<int, int>> queries;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "[ERROR] Kann Query-Datei nicht öffnen: " << filename << "\n";
        return queries;
    }
    int s, t;
    while (in >> s >> t) queries.emplace_back(s, t);
    return queries;
}

void generateRandomQueries(const std::string& filename,
    const std::vector<Node>& nodes,
    int numQueries)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "[ERROR] Kann Datei nicht schreiben: " << filename << "\n";
        return;
    }
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, (int)nodes.size() - 1);

    for (int i = 0; i < numQueries; ++i) {
        int a = dist(rng), b = dist(rng);
        while (a == b) b = dist(rng);
        out << nodes[a].id << " " << nodes[b].id << "\n";
    }
    std::cout << "[INFO] " << numQueries << " Zufalls-Queries gespeichert in " << filename << "\n";
}

void append_csv_row(const std::string& filepath,
    int qidx,
    const char* algo,
    int source, int target,
    bool found,
    long long time_ms,
    double dist_target,
    size_t path_nodes,
    double source_lat, double source_lon,
    double target_lat, double target_lon)
{
    bool need_header = true;
    {
        std::ifstream in(filepath);
        if (in.good() && in.peek() != std::ifstream::traits_type::eof())
            need_header = false;
    }

    std::ofstream out(filepath, std::ios::app);
    if (!out) return;

    if (need_header) {
        out << "qidx,algo,source,target,source_lat,source_lon,target_lat,target_lon,"
            "found,time_ms,dist_target,path_nodes\n";
    }

    out << qidx << ","
        << algo << ","
        << source << ","
        << target << ","
        << std::setprecision(10) << source_lat << ","
        << std::setprecision(10) << source_lon << ","
        << std::setprecision(10) << target_lat << ","
        << std::setprecision(10) << target_lon << ","
        << (found ? 1 : 0) << ","
        << time_ms << ","
        << std::setprecision(10) << dist_target << ","
        << path_nodes << "\n";
}

// Pfad aus predecessor-Array rekonstruieren
std::vector<int> rebuildPath(int source, int target, const std::vector<int>& prev) {
    std::vector<int> path;
    if (target < 0 || target >= (int)prev.size()) return path;
    int cur = target;
    while (cur != -1) {
        path.push_back(cur);
        if (cur == source) break;
        cur = prev[cur];
        if (cur < 0) break;
    }
    if (path.empty() || path.back() != source) { path.clear(); return path; }
    std::reverse(path.begin(), path.end());
    return path;
}

std::string pathIdsCsv(const std::vector<int>& path) {
    std::ostringstream oss;
    for (size_t i = 0; i < path.size(); ++i) {
        if (i) oss << ' ';
        oss << path[i];
    }
    return oss.str();
}

std::string pathCoordsCsv(const std::vector<int>& path,
    const std::unordered_map<int, Node>& nodeMap)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(7);
    bool first = true;
    for (int id : path) {
        auto it = nodeMap.find(id);
        if (it == nodeMap.end()) continue;
        if (!first) oss << ' ';
        first = false;
        oss << it->second.lat << ":" << it->second.lon;
    }
    return oss.str();
}

// ===== FMI-Konvertierung / Parsing =====
std::pair<std::vector<Node>, std::vector<Edge>>
convertGraphToFmiFormat(const std::vector<geo::Point>& points,
    const std::unordered_map<int, std::unordered_set<int>>& graph)
{
    std::vector<Node> nodes;
    std::vector<Edge> edges;

    // Punkte (in RAD) -> Nodes (in DEG)
    nodes.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        nodes.push_back(Node{
            static_cast<int>(i),
            points[i].y * 180.0 / geo::PI,  // lat
            points[i].x * 180.0 / geo::PI   // lon
            });
    }

    // Kanten
    for (const auto& entry : graph) {
        int a = entry.first;
        const auto& neighbors = entry.second;
        for (int b : neighbors) {
            if (a < b) {
                double dist = geo::haversine(
                    points[a].y * 180.0 / geo::PI, points[a].x * 180.0 / geo::PI,
                    points[b].y * 180.0 / geo::PI, points[b].x * 180.0 / geo::PI
                );
                edges.push_back(Edge{ a, b, dist });
                edges.push_back(Edge{ b, a, dist });
            }
        }
    }
    return { nodes, edges };
}

std::pair<
    std::unordered_map<int, Node>,
    std::vector<std::vector<std::pair<int, double>>>>
    parseFMIFile(const std::string& filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open File: " << filePath << std::endl;
        return { {}, {} };
    }

    std::unordered_map<int, Node> nodeMap;
    std::vector<std::vector<std::pair<int, double>>> adjList;
    std::string line;
    int num_nodes = 0, num_edges = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        num_nodes = std::stoi(line);
        break;
    }
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        num_edges = std::stoi(line);
        break;
    }

    adjList.resize(num_nodes + 1);

    int nodes_read = 0;
    while (nodes_read < num_nodes && std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        Node n;
        iss >> n.id >> n.lat >> n.lon;
        nodeMap[n.id] = n;
        nodes_read++;
    }

    int edges_read = 0;
    while (edges_read < num_edges && std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int source_id, target_id; double distance;
        iss >> source_id >> target_id >> distance;
        adjList[source_id].emplace_back(target_id, distance);
        edges_read++;
    }

    file.close();
    std::cout << "[DEBUG] FMI geladen: nodes=" << num_nodes << ", edges=" << num_edges << std::endl;
    std::cout << "[DEBUG] Tatsächlich gelesen: nodes=" << nodes_read << ", edges=" << edges_read << std::endl;
    return { nodeMap, adjList };
}

// ===== FMI speichern =====
void saveToFmi(const std::string& filename,
    const std::vector<Node>& nodes,
    const std::vector<Edge>& edges)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "[FATAL] Konnte Datei nicht öffnen: " << filename << std::endl;
        return;
    }

    out << nodes.size() << "\n";
    out << edges.size() << "\n";

    for (size_t i = 0; i < nodes.size(); ++i) {
        out << i << " " << nodes[i].lat << " " << nodes[i].lon << "\n";
    }
    for (const auto& e : edges) {
        out << e.a << " " << e.b << " " << e.length << "\n";
    }

    std::cout << "[ERGEBNIS] Knoten: " << nodes.size() << std::endl;
    std::cout << "[ERGEBNIS] Kanten: " << edges.size() << std::endl;
    std::cout << "[ERGEBNIS] Ø Kanten pro Knoten (bidirektional): "
        << (edges.size() / static_cast<double>(nodes.size())) << std::endl;
}

// ===== Landmark-Tabellen I/O =====
void saveDistancesToBinary(const std::string& filename,
    const std::vector<double>& distances)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "[ERROR] Konnte Datei nicht schreiben: " << filename << "\n";
        return;
    }
    size_t n = distances.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(distances.data()), n * sizeof(double));
    std::cout << "[INFO] Distanztabelle gespeichert: " << filename
        << " (" << n << " Werte, " << (n * sizeof(double) / 1024 / 1024.0) << " MB)\n";
}

// Speichert Distanzvektor als float32, binär, ohne Header.
// Werte: finite double -> float; INF/NAN bleiben als +inf/NaN in float.
bool saveDistancesToBinaryF32(const std::string& path,
    const std::vector<double>& dist)
{
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "[ERROR] Kann Datei nicht schreiben: " << path << "\n";
        return false;
    }
    // Cast-on-the-fly in Blöcken, um Peak-RAM klein zu halten
    constexpr size_t BLK = 1 << 16; // 65k
    std::vector<float> buf; buf.resize(std::min(BLK, dist.size()));
    size_t i = 0;
    while (i < dist.size()) {
        size_t n = std::min(BLK, dist.size() - i);
        for (size_t j = 0; j < n; ++j) {
            double d = dist[i + j];
            // bewusst keine Clamp — float::inf/NaN sind ok & eindeutig
            buf[j] = static_cast<float>(d);
        }
        out.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(n * sizeof(float)));
        i += n;
    }
    return true;
}


// Liest sowohl alte (Header + double) als auch neue (roh float32) Tabellen.
std::vector<double> loadDistancesFromBinary(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary | std::ios::ate);
    if (!in) {
        // keine Datei → still return {}
        return {};
    }

    const std::streamsize bytes = in.tellg();
    if (bytes <= 0) return {};
    in.seekg(0, std::ios::beg);

    // Variante A: altes Format (Header + double)
    {
        if (bytes > static_cast<std::streamsize>(sizeof(size_t))) {
            size_t n_header = 0;
            std::streampos p0 = in.tellg();
            in.read(reinterpret_cast<char*>(&n_header), sizeof(size_t));
            const std::streamsize remain = bytes - static_cast<std::streamsize>(sizeof(size_t));
            if (n_header > 0 && remain == static_cast<std::streamsize>(n_header * sizeof(double))) {
                std::vector<double> out(n_header);
                in.read(reinterpret_cast<char*>(out.data()), remain);
                return out;
            }
            // Header passt nicht → zurückspulen für die Roh-Fälle
            in.clear();
            in.seekg(p0, std::ios::beg);
        }
    }

    // Variante B: raw Float32 (neues Format, ohne Header)
    if (bytes % 4 == 0) {
        const size_t n = static_cast<size_t>(bytes / 4);
        std::vector<float> buf(n);
        if (in.read(reinterpret_cast<char*>(buf.data()), bytes)) {
            std::vector<double> out(n);
            for (size_t i = 0; i < n; ++i) out[i] = static_cast<double>(buf[i]);
            return out;
        }
        in.clear(); in.seekg(0, std::ios::beg);
    }

    // Variante C: raw Double (ohne Header)
    if (bytes % 8 == 0) {
        const size_t n = static_cast<size_t>(bytes / 8);
        std::vector<double> out(n);
        if (in.read(reinterpret_cast<char*>(out.data()), bytes)) {
            return out;
        }
    }

    // wenn keine Variante gepasst hat → leeres Ergebnis
    return {};
}


// ===== ALT / Search =====
void initAltGlobals(std::size_t n) {
    g_alt_distances.assign(n, std::numeric_limits<double>::infinity());
    g_alt_previous.assign(n, -1);
    g_alt_touched.clear();
    g_alt_touched.reserve(n / 10);
}

void resetAltTouched() {
    for (int v : g_alt_touched) {
        g_alt_distances[v] = std::numeric_limits<double>::infinity();
        g_alt_previous[v] = -1;
    }
    g_alt_touched.clear();
}

void resetAltNode(int v) {
    g_alt_distances[v] = std::numeric_limits<double>::infinity();
    g_alt_previous[v] = -1;
    g_alt_touched.push_back(v);
}

// ===== Such-Initialisierung =====
void initializeSearchStructures(int source,
    int target,
    const WeightedAdjList& adjList,
    std::vector<double>& distances,
    std::vector<int>& previous)
{
    int maxId = std::max(source, target);
    for (const auto& kv : adjList) {
        int u = kv.first;
        const auto& nbrs = kv.second;
        maxId = std::max(maxId, u);
        for (const auto& [v, _] : nbrs) {
            (void)_; // Gewicht nicht genutzt
            maxId = std::max(maxId, v);
        }
    }
    maxId += 1;

    distances.assign(maxId, std::numeric_limits<double>::infinity());
    previous.assign(maxId, -1);

    if (source >= maxId || target >= maxId) {
        std::cerr << "[FATAL] Source/Target außerhalb des Arrays! "
            << "source=" << source << ", target=" << target
            << ", maxIndex=" << maxId << "\n";
        std::exit(1);
    }
}
