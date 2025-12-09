#include <vector>
#include <unordered_map>
#include <functional>
#include <queue>
#include <limits>
#include <cmath>
#include <cstdint>
#include <random>
#include <chrono>
#include <iostream>

#include <alg_dynamic_dijkstra.h>
#include <alg_weighted_digraph.h>

using namespace std;

// ---------- Synthetic workload generation + experiment harness ----------

struct Update {
    int u, v;
    double newW;
    double oldW;
    bool existed;
};

static uint64_t rngSeed = 1234567;
static mt19937_64 rng(rngSeed);

WeightedDigraph generateRandomGraph(int n, int m, double maxW) {
    WeightedDigraph G(n);
    uniform_int_distribution<int> nodeDist(0, n-1);
    uniform_real_distribution<double> wDist(1.0, maxW);

    // Connectivity spine (simple chain)
    for (int i = 0; i < n-1; i++) {
        double w = wDist(rng);
        G.addEdge(i, i+1, w);
    }

    int added = n-1;
    while (added < m) {
        int u = nodeDist(rng);
        int v = nodeDist(rng);
        if (u == v) continue;
        if (G.hasEdge(u, v)) continue;

        double w = wDist(rng);
        G.addEdge(u, v, w);
        added++;
    }
    return G;
}

Update randomEdgeUpdate(WeightedDigraph& G, double maxW, double maxDelta) {
    int n = G.size();
    uniform_int_distribution<int> nodeDist(0, n-1);
    uniform_real_distribution<double> deltaDist(-maxDelta, maxDelta);

    // Pick an existing edge by random trial
    int u=-1, v=-1;
    for (;;) {
        u = nodeDist(rng);
        const auto& nbrs = G.neighbors(u);
        if (nbrs.empty()) continue;

        uniform_int_distribution<int> idxDist(0, (int)nbrs.size()-1);
        int idx = idxDist(rng);
        v = nbrs[idx].to;
        break;
    }

    double oldW = G.getWeight(u, v);
    double delta = deltaDist(rng);
    double newW = oldW + delta;

    if (newW < 0.1) newW = 0.1;
    if (newW > maxW) newW = maxW;

    G.updateEdgeWeight(u, v, newW);

    Update up;
    up.u=u; up.v=v; up.oldW=oldW; up.newW=newW;
    up.existed=true;
    return up;
}

bool approxEqual(const vector<double>& a, const vector<double>& b, double eps=1e-9) {
    if (a.size() != b.size()) return false;
    for (size_t i=0; i<a.size(); i++) {
        double x=a[i], y=b[i];
        if (isinf(x) && isinf(y)) continue;
        if (fabs(x-y) > eps) return false;
    }
    return true;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Experiment parameters
    int n = 2000;
    int m = 30000;          // density lever
    int source = 0;
    int numUpdates = 1000; // streaming updates
    double maxW = 100.0;
    double maxDelta = 30.0;

    cout << "Generating graph n=" << n << " m=" << m << "\n";
    WeightedDigraph G = generateRandomGraph(n, m, maxW);

    DynamicDijkstra engine(G, source);

    cout << "Initial full recompute...\n";
    auto t0 = chrono::high_resolution_clock::now();
    DijkstraState fullState = engine.recompute();
    auto t1 = chrono::high_resolution_clock::now();
    DijkstraState incState = fullState;

    double initialMs = chrono::duration<double, milli>(t1-t0).count();
    cout << "Initial Dijkstra runtime: " << initialMs << " ms\n\n";

    long long fullMsTotal = 0;
    long long incMsTotal = 0;
    int incApplied = 0;
    int incFallbacks = 0;

    for (int i=1; i<=numUpdates; i++) {
        Update up = randomEdgeUpdate(G, maxW, maxDelta);

        bool isDecrease = (up.newW < up.oldW);

        // ----- Baseline: full recompute -----
        auto f0 = chrono::high_resolution_clock::now();
        fullState = engine.recompute();
        auto f1 = chrono::high_resolution_clock::now();
        fullMsTotal += (long long)chrono::duration<double, milli>(f1-f0).count();

        // ----- Incremental path: decrease-localized or increase-fallback -----
        auto i0 = chrono::high_resolution_clock::now();

        if (isDecrease) {
            engine.incrementalDecrease(incState, up.u, up.v);
            incApplied++;
        } else {
            bool needsFallback = engine.requiresRecomputeOnIncrease(incState, up.u, up.v);
            if (needsFallback) {
                incState = engine.recompute();
                incFallbacks++;
            }
        }

        auto i1 = chrono::high_resolution_clock::now();
        incMsTotal += (long long)chrono::duration<double, milli>(i1-i0).count();

        // Correctness gate
        if (!approxEqual(fullState.dist, incState.dist)) {
            cerr << "Mismatch detected at update " << i
                 << " edge(" << up.u << "," << up.v << ")"
                 << " oldW=" << up.oldW << " newW=" << up.newW << "\n";
            return 1;
        }

        if (i % 200 == 0) {
            cout << "[Checkpoint " << i << "] cumulative full=" << fullMsTotal
                 << " ms, incremental=" << incMsTotal << " ms\n";
        }
    }

    cout << "\n=== Experiment Summary ===\n";
    cout << "Updates processed: " << numUpdates << "\n";
    cout << "Incremental decreases applied: " << incApplied << "\n";
    cout << "Incremental increase fallbacks: " << incFallbacks << "\n";
    cout << "Total full recompute time: " << fullMsTotal << " ms\n";
    cout << "Total incremental time:    " << incMsTotal << " ms\n";

    double speedup = (incMsTotal > 0) ? (double)fullMsTotal / (double)incMsTotal : 0.0;
    cout << "Observed speedup (full/inc): " << speedup << "x\n";

    return 0;
}
