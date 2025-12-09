#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "alg_weighted_digraph.h"
#include "alg_dynamic_dijkstra.h"

#include <cmath>
#include <limits>

using namespace std;

static bool is_inf(double x) {
    return std::isinf(x);
}

static bool approx_equal(double a, double b, double eps = 1e-9) {
    if (is_inf(a) && is_inf(b)) return true;
    return std::fabs(a - b) < eps;
}

TEST_CASE("WeightedDigraph basic construction and edge operations", "[graph]") {
    int n = 4;
    WeightedDigraph G(n);

    SECTION("Graph size and initial state") {
        REQUIRE(G.size() == n);
        for (int u = 0; u < n; ++u) {
            REQUIRE(G.neighbors(u).empty());
        }
    }

    SECTION("Add edge and query weight") {
        G.addEdge(0, 1, 3.5);
        REQUIRE(G.hasEdge(0, 1));
        REQUIRE_FALSE(G.hasEdge(1, 0));

        double w01 = G.getWeight(0, 1);
        CHECK(approx_equal(w01, 3.5));

        // Nonexistent edge should be INF
        double w10 = G.getWeight(1, 0);
        CHECK(is_inf(w10));
    }

    SECTION("Update edge weight and reuse") {
        G.addEdge(0, 1, 5.0);
        double oldW = G.updateEdgeWeight(0, 1, 2.0);
        CHECK(approx_equal(oldW, 5.0));
        CHECK(approx_equal(G.getWeight(0, 1), 2.0));

        // Updating a non-existing edge should behave like add
        double oldW2 = G.updateEdgeWeight(1, 2, 4.0);
        CHECK(is_inf(oldW2));
        CHECK(G.hasEdge(1, 2));
        CHECK(approx_equal(G.getWeight(1, 2), 4.0));
    }

    SECTION("Remove edge") {
        G.addEdge(0, 1, 1.0);
        G.addEdge(0, 2, 2.0);
        REQUIRE(G.hasEdge(0, 1));
        REQUIRE(G.hasEdge(0, 2));

        G.removeEdge(0, 1);
        CHECK_FALSE(G.hasEdge(0, 1));
        CHECK(G.hasEdge(0, 2)); // other edge still intact

        CHECK(is_inf(G.getWeight(0, 1)));
        CHECK(approx_equal(G.getWeight(0, 2), 2.0));
    }
}

TEST_CASE("Dijkstra full recompute on a simple graph", "[dijkstra][baseline]") {
    // Graph:
    // 0 -> 1 (1)
    // 0 -> 2 (10)
    // 1 -> 2 (5)
    //
    // Shortest distances from 0:
    // dist[0] = 0
    // dist[1] = 1
    // dist[2] = 6  via 0 -> 1 -> 2
    WeightedDigraph G(3);
    G.addEdge(0, 1, 1.0);
    G.addEdge(0, 2, 10.0);
    G.addEdge(1, 2, 5.0);

    DynamicDijkstra engine(G, 0);
    DijkstraState st = engine.recompute();

    REQUIRE(st.dist.size() == 3);
    REQUIRE(st.parent.size() == 3);

    CHECK(approx_equal(st.dist[0], 0.0));
    CHECK(approx_equal(st.dist[1], 1.0));
    CHECK(approx_equal(st.dist[2], 6.0));

    // Parent relationships (one valid SPT option):
    CHECK(st.parent[0] == -1);
    CHECK(st.parent[1] == 0);
    CHECK(st.parent[2] == 1);
}

TEST_CASE("Incremental decrease updates downstream shortest paths correctly", "[dijkstra][incremental]") {
    // Initial graph:
    // 0 -> 1 (1)
    // 0 -> 2 (10)
    // 1 -> 2 (5)
    //
    // Initial shortest paths from 0:
    // dist[0] = 0
    // dist[1] = 1
    // dist[2] = 6 (via 0->1->2)
    //
    // Then decrease edge (1,2) from 5 to 1, new best:
    // dist[2] = 2 (via 0->1->2)
    WeightedDigraph G(3);
    G.addEdge(0, 1, 1.0);
    G.addEdge(0, 2, 10.0);
    G.addEdge(1, 2, 5.0);

    DynamicDijkstra engine(G, 0);
    DijkstraState st = engine.recompute();

    REQUIRE(approx_equal(st.dist[2], 6.0));

    // Apply a weight decrease on (1,2)
    double oldW = G.updateEdgeWeight(1, 2, 1.0);
    REQUIRE(approx_equal(oldW, 5.0));

    engine.incrementalDecrease(st, 1, 2);

    // Validate updated distances
    CHECK(approx_equal(st.dist[0], 0.0));
    CHECK(approx_equal(st.dist[1], 1.0));
    CHECK(approx_equal(st.dist[2], 2.0));

    CHECK(st.parent[0] == -1);
}