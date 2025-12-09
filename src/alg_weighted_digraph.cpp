#include "alg_weighted_digraph.h"

using namespace std;

long long WeightedDigraph::key(int u, int v) const {
    return 1LL * u * n + v;
}

WeightedDigraph::WeightedDigraph(int nVertices)
    : n(nVertices), adj(nVertices) {}

int WeightedDigraph::size() const {
    return n;
}

const vector<Edge>& WeightedDigraph::neighbors(int u) const {
    return adj[u];
}

bool WeightedDigraph::hasEdge(int u, int v) const {
    return edgeIndex.find(key(u, v)) != edgeIndex.end();
}

double WeightedDigraph::getWeight(int u, int v) const {
    auto it = edgeIndex.find(key(u, v));
    if (it == edgeIndex.end()) return numeric_limits<double>::infinity();
    return adj[u][it->second].w;
}

void WeightedDigraph::addEdge(int u, int v, double w) {
    long long k = key(u, v);
    if (edgeIndex.find(k) != edgeIndex.end()) {
        updateEdgeWeight(u, v, w);
        return;
    }
    edgeIndex[k] = (int)adj[u].size();
    adj[u].emplace_back(v, w);
}

double WeightedDigraph::updateEdgeWeight(int u, int v, double newW) {
    long long k = key(u, v);
    auto it = edgeIndex.find(k);
    if (it == edgeIndex.end()) {
        addEdge(u, v, newW);
        return numeric_limits<double>::infinity();
    }
    int idx = it->second;
    double oldW = adj[u][idx].w;
    adj[u][idx].w = newW;
    return oldW;
}

void WeightedDigraph::removeEdge(int u, int v) {
    long long k = key(u, v);
    auto it = edgeIndex.find(k);
    if (it == edgeIndex.end()) return;

    int idx = it->second;
    int lastIdx = (int)adj[u].size() - 1;

    if (idx != lastIdx) {
        adj[u][idx] = adj[u][lastIdx];
        edgeIndex[key(u, adj[u][idx].to)] = idx;
    }

    adj[u].pop_back();
    edgeIndex.erase(it);
}
