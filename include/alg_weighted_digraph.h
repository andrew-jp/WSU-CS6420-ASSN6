#ifndef WEIGHTED_DIGRAPH_H
#define WEIGHTED_DIGRAPH_H

#include <vector>
#include <unordered_map>
#include <limits>

struct Edge {
    int to;
    double w;
    Edge(int t, double wt) : to(t), w(wt) {}
};

class WeightedDigraph {
private:
    int n;
    std::vector<std::vector<Edge>> adj;

    // Fast updates: key = u*n + v  => index in adj[u]
    std::unordered_map<long long, int> edgeIndex;

    long long key(int u, int v) const;

public:
    explicit WeightedDigraph(int nVertices);

    int size() const;

    const std::vector<Edge>& neighbors(int u) const;

    bool hasEdge(int u, int v) const;

    double getWeight(int u, int v) const;

    void addEdge(int u, int v, double w);

    // Returns old weight if existed, else INF
    double updateEdgeWeight(int u, int v, double newW);

    void removeEdge(int u, int v);
};

#endif // WEIGHTED_DIGRAPH_H
