#include "graph.h"
#include <queue>
#include <algorithm>
#include <set>

void CampusGraph::addSpot(const Spot& spot) {
    spots_.push_back(spot);
    codeToId_[spot.code] = spot.id;
    if (adj_.size() <= static_cast<size_t>(spot.id)) {
        adj_.resize(spot.id + 1);
    }
}

void CampusGraph::addEdge(int from, int to, double distance) {
    size_t maxId = std::max(from, to);
    if (adj_.size() <= maxId) {
        adj_.resize(maxId + 1);
    }
    adj_[from].push_back({to, distance});
    adj_[to].push_back({from, distance});
}

const std::vector<Spot>& CampusGraph::getSpots() const {
    return spots_;
}

const Spot* CampusGraph::getSpotById(int id) const {
    for (const auto& spot : spots_) {
        if (spot.id == id) return &spot;
    }
    return nullptr;
}

const Spot* CampusGraph::getSpotByCode(const std::string& code) const {
    auto it = codeToId_.find(code);
    if (it != codeToId_.end()) {
        return getSpotById(it->second);
    }
    return nullptr;
}

PathResult CampusGraph::dijkstra(int startId, int endId) const {
    PathResult result;
    result.found = false;
    result.totalDistance = 0;
    result.algorithmName = "";
    result.elapsedMs = 0;
    
    if (startId < 0 || startId >= static_cast<int>(adj_.size()) ||
        endId < 0 || endId >= static_cast<int>(adj_.size())) {
        return result;
    }
    
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> dist(adj_.size(), INF);
    std::vector<int> prev(adj_.size(), -1);
    
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> pq;
    
    dist[startId] = 0;
    pq.push({0, startId});
    
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (d > dist[u]) continue;
        if (u == endId) break;
        
        for (const auto& [v, w] : adj_[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
    
    if (dist[endId] == INF) return result;
    
    result.found = true;
    result.totalDistance = dist[endId];
    
    int current = endId;
    while (current != -1) {
        result.path.push_back(current);
        current = prev[current];
    }
    std::reverse(result.path.begin(), result.path.end());
    
    return result;
}

PathResult CampusGraph::astar(int startId, int endId) const {
    PathResult result;
    result.found = false;
    result.totalDistance = 0;
    
    if (startId < 0 || startId >= static_cast<int>(adj_.size()) ||
        endId < 0 || endId >= static_cast<int>(adj_.size())) {
        return result;
    }
    
    const Spot* startSpot = getSpotById(startId);
    const Spot* endSpot = getSpotById(endId);
    if (!startSpot || !endSpot) {
        return result;
    }
    
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> gScore(adj_.size(), INF);  // 从起点到当前节点的实际距离
    std::vector<int> prev(adj_.size(), -1);
    
    // 优先队列，存储 (f_score, 节点ID)，f = g + h
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> pq;
    
    gScore[startId] = 0;
    double hScore = std::hypot(startSpot->x - endSpot->x, startSpot->y - endSpot->y);
    pq.push({hScore, startId});  // f = g + h = 0 + h
    
    while (!pq.empty()) {
        auto [currentF, u] = pq.top();
        pq.pop();
        
        if (u == endId) break;
        
        // 如果当前 f_score 比记录的大，说明已经有更优路径了
        // 因为 f = g + h，而 g[u] 是最小的已知实际距离
        // 所以如果 currentF > g[u] + h(u)，说明这是一个过期的节点
        const Spot* currentSpot = getSpotById(u);
        if (!currentSpot) continue;
        
        double currentH = std::hypot(currentSpot->x - endSpot->x, currentSpot->y - endSpot->y);
        if (currentF > gScore[u] + currentH) continue;
        
        for (const auto& [v, w] : adj_[u]) {
            double tentativeG = gScore[u] + w;
            
            if (tentativeG < gScore[v]) {
                gScore[v] = tentativeG;
                prev[v] = u;
                
                const Spot* currentSpot = getSpotById(v);
                if (currentSpot) {
                    hScore = std::hypot(currentSpot->x - endSpot->x, currentSpot->y - endSpot->y);
                    double fScore = tentativeG + hScore;
                    pq.push({fScore, v});
                }
            }
        }
    }
    
    if (gScore[endId] == INF) return result;
    
    result.found = true;
    result.totalDistance = gScore[endId];
    
    int current = endId;
    while (current != -1) {
        result.path.push_back(current);
        current = prev[current];
    }
    std::reverse(result.path.begin(), result.path.end());
    
    return result;
}

PathResult CampusGraph::multiSpotPath(const std::vector<int>& spotIds) const {
    PathResult result;
    result.found = false;
    result.totalDistance = 0;
    
    if (spotIds.empty()) return result;
    if (spotIds.size() == 1) {
        result.found = true;
        result.path = spotIds;
        return result;
    }
    
    std::set<int> unvisited(spotIds.begin() + 1, spotIds.end());
    int current = spotIds[0];
    result.path.push_back(current);
    
    while (!unvisited.empty()) {
        int nearest = -1;
        double minDist = std::numeric_limits<double>::infinity();
        PathResult bestPath;
        
        for (int next : unvisited) {
            PathResult pr = dijkstra(current, next);
            if (pr.found && pr.totalDistance < minDist) {
                minDist = pr.totalDistance;
                nearest = next;
                bestPath = pr;
            }
        }
        
        if (nearest == -1) return result;
        
        for (size_t i = 1; i < bestPath.path.size(); ++i) {
            result.path.push_back(bestPath.path[i]);
        }
        result.totalDistance += minDist;
        
        unvisited.erase(nearest);
        current = nearest;
    }
    
    result.found = true;
    return result;
}

std::vector<Edge> CampusGraph::getEdges() const {
    std::vector<Edge> edges;
    std::set<std::pair<int, int>> added;
    
    for (size_t i = 0; i < adj_.size(); ++i) {
        for (const auto& [j, dist] : adj_[i]) {
            int minId = std::min(static_cast<int>(i), j);
            int maxId = std::max(static_cast<int>(i), j);
            if (added.find({minId, maxId}) == added.end()) {
                edges.push_back({static_cast<int>(i), j, dist});
                added.insert({minId, maxId});
            }
        }
    }
    
    return edges;
}
