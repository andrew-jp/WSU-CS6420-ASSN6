#include "alg_dynamic_dijkstra.h"
#include <queue>
#include <cmath>

using namespace std;

DynamicDijkstra::DynamicDijkstra(const WeightedDigraph& graph, int s)
    : G(graph), source(s) {}

DijkstraState DynamicDijkstra::recompute() const {
    int n = G.size();
    DijkstraState st;
    st.dist.assign(n, INF);
    st.parent.assign(n, -1);

    using P = pair<double,int>;
    priority_queue<P, vector<P>, greater<P>> pq;

    st.dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        auto [d,u] = pq.top(); pq.pop();
        if (d != st.dist[u]) continue;

        for (const auto& e : G.neighbors(u)) {
            int v = e.to;
            double nd = d + e.w;
            if (nd < st.dist[v]) {
                st.dist[v] = nd;
                st.parent[v] = u;
                pq.push({nd, v});
            }
        }
    }

    return st;
}

void DynamicDijkstra::incrementalDecrease(DijkstraState& st, int u, int v) const {
    using P = pair<double,int>;
    priority_queue<P, vector<P>, greater<P>> pq;

    double w = G.getWeight(u, v);
    if (isinf(w) || isinf(st.dist[u])) return;

    double candidate = st.dist[u] + w;
    if (candidate >= st.dist[v]) return;

    st.dist[v] = candidate;
    st.parent[v] = u;
    pq.push({candidate, v});

    while (!pq.empty()) {
        auto [d,x] = pq.top(); pq.pop();
        if (d != st.dist[x]) continue;

        for (const auto& e : G.neighbors(x)) {
            int y = e.to;
            double nd = d + e.w;
            if (nd < st.dist[y]) {
                st.dist[y] = nd;
                st.parent[y] = x;
                pq.push({nd, y});
            }
        }
    }
}

bool DynamicDijkstra::requiresRecomputeOnIncrease(const DijkstraState& st, int u, int v) const {
    return st.parent[v] == u;
}
