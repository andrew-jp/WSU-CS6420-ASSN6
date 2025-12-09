#ifndef DYNAMIC_DIJKSTRA_H
#define DYNAMIC_DIJKSTRA_H

#include "alg_weighted_digraph.h"
#include <vector>

struct DijkstraState {
    std::vector<double> dist;
    std::vector<int> parent; // parent[v] = predecessor on SPT; -1 for root/unreached
};

class DynamicDijkstra {
private:
    const WeightedDigraph& G;
    int source;
    static constexpr double INF = std::numeric_limits<double>::infinity();

public:
    DynamicDijkstra(const WeightedDigraph& graph, int s);

    // Full recompute baseline
    DijkstraState recompute() const;

    // Incremental update for DECREASE or insertion on edge (u,v)
    void incrementalDecrease(DijkstraState& st, int u, int v) const;

    // For INCREASE on edge (u,v):
    // recompute required iff edge is on current SPT (parent[v] == u)
    bool requiresRecomputeOnIncrease(const DijkstraState& st, int u, int v) const;
};

#endif // DYNAMIC_DIJKSTRA_H
