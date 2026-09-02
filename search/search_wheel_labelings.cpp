// C++ port of search_wheel_labelings.py — the same annealer, the same two-phase
// method (Proposition 3.2 frame-completion), the same CLI and output format, built
// for the long open cases W_14..W_40 where interpreter overhead dominates.
//
//   g++ -O3 -march=native -std=c++17 -o wheel_search_cpp search_wheel_labelings.cpp
//
// Parity with the Python script, and the three deliberate divergences:
//   * The algorithm, cost function, move set, temperature schedule (step-based,
//     cycle 20000), quotas/assembly/backtracking and rim completion mirror the
//     Python code function for function.
//   * RNG streams differ (std::mt19937_64 here, Python's Mersenne Twister
//     there), so trajectories are not step-identical across the two programs;
//     each program is still fully deterministic given --seed.
//   * Where Python iterates dicts in insertion order (quota tie order among
//     equal-sized collision groups), this port uses ascending value order.
//   * The "(a labeling not in Table 5)" annotation is not printed here; use the
//     Python script for that comparison.
// Every labeling this program prints has been re-verified internally by
// recomputing all four geodesic weight classes from scratch (is_geodesic_leech
// below, a direct port of the Python check). Re-verify externally with the
// Python is_geodesic_leech as well; a search engine should not be trusted, only
// its certificates.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <vector>

using std::vector;
using i64 = long long;

static int geodesic_count(int m) { return m * (m + 3) / 2; }  // t_gp(W_{m+1})
static bool g_dump_sets = false;  // --dump-sets: print every cost-zero set

// ---------------------------------------------------------------------------
// RNG: deterministic given seed; stream differs from Python's random.Random.
// ---------------------------------------------------------------------------
struct Rng {
    std::mt19937_64 g;
    explicit Rng(uint64_t seed) : g(seed) {}
    int randrange(int n) { return (int)(g() % (uint64_t)n); }          // [0, n)
    int randint(int a, int b) { return a + (int)(g() % (uint64_t)(b - a + 1)); }
    double random() { return (g() >> 11) * (1.0 / 9007199254740992.0); }
    template <class T> void shuffle(vector<T>& v) {
        for (int i = (int)v.size() - 1; i > 0; --i) {
            int j = randrange(i + 1);
            std::swap(v[i], v[j]);
        }
    }
    vector<int> sample_1_to(int N, int m) {  // m distinct values from [1, N]
        vector<int> vals(N);
        for (int i = 0; i < N; ++i) vals[i] = i + 1;
        for (int i = 0; i < m; ++i) std::swap(vals[i], vals[i + randrange(N - i)]);
        vals.resize(m);
        return vals;
    }
};

// ---------------------------------------------------------------------------
// Phase 1a: annealing over spoke sets (port of class SpokeSet)
// ---------------------------------------------------------------------------
struct SpokeSet {
    int m, N;
    vector<int> s;               // the m labels
    vector<uint8_t> in_set;      // [0..N]
    vector<int> total;           // [0..2N]: pairs with that sum
    vector<int> forced_count;    // [0..2N]: forced pairs with that sum
    vector<vector<int>> pairs_of;  // [0..2N] -> packed pairs (i * 64 + j, i < j)
    vector<uint8_t> forced;      // m*m, index i*m+j for i < j
    vector<int> forced_degree;   // [0..m)
    i64 required = 0;
    i64 label_sum = 0;
    int residue_target;
    int cost = 0;

    SpokeSet(int m_, int N_, const vector<int>& values)
        : m(m_), N(N_), s(values), in_set(N_ + 1, 0), total(2 * N_ + 1, 0),
          forced_count(2 * N_ + 1, 0), pairs_of(2 * N_ + 1), forced(m_ * m_, 0),
          forced_degree(m_, 0) {
        for (int v : s) in_set[v] = 1;
        label_sum = 0;
        for (int v : s) label_sum += v;
        residue_target = (int)(((i64)N * (N + 1) / 2) % 3);
        for (int i = 0; i < m; ++i)
            for (int j = i + 1; j < m; ++j) add_pair(i, j);
        cost = recompute_cost();
    }

    i64 required_of(int value) const {
        int quota = total[value] > 1 ? total[value] - 1 : 0;
        return forced_count[value] > quota ? forced_count[value] : quota;
    }

    void add_pair(int i, int j) {  // requires i < j
        int value = s[i] + s[j];
        i64 before = required_of(value);
        total[value] += 1;
        pairs_of[value].push_back(i * 64 + j);
        bool f = value > N || in_set[value];
        forced[i * m + j] = f;
        if (f) {
            forced_count[value] += 1;
            forced_degree[i] += 1;
            forced_degree[j] += 1;
        }
        required += required_of(value) - before;
    }

    void remove_pair(int i, int j) {  // requires i < j
        int value = s[i] + s[j];
        i64 before = required_of(value);
        total[value] -= 1;
        auto& bucket = pairs_of[value];
        int packed = i * 64 + j;
        for (size_t k = 0; k < bucket.size(); ++k)
            if (bucket[k] == packed) {
                bucket[k] = bucket.back();
                bucket.pop_back();
                break;
            }
        if (forced[i * m + j]) {
            forced[i * m + j] = 0;
            forced_count[value] -= 1;
            forced_degree[i] -= 1;
            forced_degree[j] -= 1;
        }
        required += required_of(value) - before;
    }

    void refresh_forced(int value) {  // membership of `value` changed
        if (value > N) return;
        bool should_be = in_set[value];
        for (int packed : pairs_of[value]) {
            int i = packed >> 6, j = packed & 63;
            if ((bool)forced[i * m + j] == should_be) continue;
            i64 before = required_of(value);
            forced[i * m + j] = should_be;
            int step = should_be ? 1 : -1;
            forced_count[value] += step;
            forced_degree[i] += step;
            forced_degree[j] += step;
            required += required_of(value) - before;
        }
    }

    int recompute_cost() const {
        i64 over_budget = required > m ? required - m : 0;
        int over_degree = 0;
        for (int d : forced_degree)
            if (d > 2) over_degree += d - 2;
        int residue = ((i64)(m - 2) * label_sum) % 3 == residue_target ? 0 : 1;
        return (int)over_budget + over_degree + residue;
    }

    int replace(int index, int value) {  // value must not already be present
        int previous = s[index];
        for (int j = 0; j < m; ++j)
            if (j != index) remove_pair(std::min(index, j), std::max(index, j));
        in_set[previous] = 0;
        s[index] = value;
        in_set[value] = 1;
        for (int j = 0; j < m; ++j)
            if (j != index) add_pair(std::min(index, j), std::max(index, j));
        refresh_forced(previous);
        refresh_forced(value);
        label_sum += value - previous;
        cost = recompute_cost();
        return previous;
    }
};

// ---------------------------------------------------------------------------
// Phase 1b: assembling a feasible set into cyclic orders (ports of PathSet,
// edge_selections, cyclic_orders)
// ---------------------------------------------------------------------------
struct PathSet {
    int m;
    vector<int> adj;  // m*2 neighbor slots, -1 = empty
    vector<int> deg;
    explicit PathSet(int m_) : m(m_), adj(m_ * 2, -1), deg(m_, 0) {}

    bool reaches(int i, int j) const {
        int previous = -1, current = i;
        for (;;) {
            int forward = -1;
            for (int k = 0; k < deg[current]; ++k)
                if (adj[current * 2 + k] != previous) {
                    forward = adj[current * 2 + k];
                    break;
                }
            if (forward < 0) return false;
            previous = current;
            current = forward;
            if (current == j || current == i) return true;
        }
    }
    bool can_add(int i, int j) const {
        return deg[i] < 2 && deg[j] < 2 && !reaches(i, j);
    }
    void add(int i, int j) {
        adj[i * 2 + deg[i]++] = j;
        adj[j * 2 + deg[j]++] = i;
    }
    void remove(int i, int j) {
        auto drop = [&](int a, int b) {
            if (adj[a * 2] == b) adj[a * 2] = adj[a * 2 + 1];
            adj[a * 2 + 1] = -1;
            --deg[a];
        };
        drop(i, j);
        drop(j, i);
    }
    vector<vector<int>> components() const {
        vector<uint8_t> seen(m, 0);
        vector<vector<int>> out;
        for (int start = 0; start < m; ++start) {
            if (seen[start] || deg[start] == 2) continue;  // start only at ends
            vector<int> path{start};
            seen[start] = 1;
            int previous = -1, current = start;
            for (;;) {
                int forward = -1;
                for (int k = 0; k < deg[current]; ++k)
                    if (adj[current * 2 + k] != previous) {
                        forward = adj[current * 2 + k];
                        break;
                    }
                if (forward < 0) break;
                previous = current;
                current = forward;
                seen[current] = 1;
                path.push_back(current);
            }
            out.push_back(std::move(path));
        }
        return out;
    }
};

// Yield up to `limit` edge selections; cb returns true to stop early.
template <class Cb>
static void edge_selections(const SpokeSet& state, int limit, Cb&& cb) {
    int m = state.m;
    PathSet paths(m);
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j) {
            if (!state.forced[i * m + j]) continue;
            if (!paths.can_add(i, j)) return;
            paths.add(i, j);
        }

    // (deficit, free pairs) per collision group, smallest group first.
    vector<std::pair<int, vector<std::pair<int, int>>>> quotas;
    for (int value = 0; value <= 2 * state.N; ++value) {
        if (state.total[value] == 0) continue;
        int deficit = (state.total[value] - 1) - state.forced_count[value];
        if (deficit <= 0) continue;
        vector<std::pair<int, int>> free;
        for (int packed : state.pairs_of[value]) {
            int i = packed >> 6, j = packed & 63;
            if (!state.forced[i * m + j]) free.emplace_back(i, j);
        }
        std::sort(free.begin(), free.end());
        if ((int)free.size() < deficit) return;  // cost should have caught this
        quotas.emplace_back(deficit, std::move(free));
    }
    std::stable_sort(quotas.begin(), quotas.end(),
                     [](const auto& a, const auto& b) { return a.second.size() < b.second.size(); });

    int produced = 0;
    bool stop = false;

    // Backtrack over combinations(free, deficit) per group.
    std::function<void(size_t)> descend = [&](size_t k) {
        if (stop || produced >= limit) return;
        if (k == quotas.size()) {
            ++produced;
            if (cb(paths)) stop = true;
            return;
        }
        int deficit = quotas[k].first;
        const auto& free = quotas[k].second;
        int n = (int)free.size();
        vector<int> pick(deficit);
        for (int t = 0; t < deficit; ++t) pick[t] = t;
        for (;;) {
            vector<std::pair<int, int>> added;
            bool feasible = true;
            for (int t = 0; t < deficit; ++t) {
                auto [i, j] = free[pick[t]];
                if (!paths.can_add(i, j)) {
                    feasible = false;
                    break;
                }
                paths.add(i, j);
                added.emplace_back(i, j);
            }
            if (feasible) descend(k + 1);
            for (auto it = added.rbegin(); it != added.rend(); ++it)
                paths.remove(it->first, it->second);
            if (stop || produced >= limit) return;
            // next combination in lexicographic order
            int t = deficit - 1;
            while (t >= 0 && pick[t] == n - deficit + t) --t;
            if (t < 0) break;
            ++pick[t];
            for (int u = t + 1; u < deficit; ++u) pick[u] = pick[u - 1] + 1;
        }
    };
    descend(0);
}

// Yield admissible-candidate cyclic orders of state.s; cb returns true to stop.
template <class Cb>
static void cyclic_orders(const SpokeSet& state, Rng& rng, int selections, Cb&& cb) {
    const int linkings = 4;
    edge_selections(state, selections, [&](const PathSet& paths) {
        auto components = paths.components();
        for (int rep = 0; rep < linkings; ++rep) {
            vector<int> order;
            vector<int> idx(components.size());
            for (size_t i = 0; i < idx.size(); ++i) idx[i] = (int)i;
            rng.shuffle(idx);
            for (int ci : idx) {
                const auto& path = components[ci];
                if (rng.random() < 0.5)
                    order.insert(order.end(), path.rbegin(), path.rend());
                else
                    order.insert(order.end(), path.begin(), path.end());
            }
            if ((int)order.size() != state.m) break;  // should not happen
            vector<int> labels(state.m);
            for (int i = 0; i < state.m; ++i) labels[i] = state.s[order[i]];
            if (cb(labels)) return true;
        }
        return false;
    });
}

// ---------------------------------------------------------------------------
// Phase 2: rim completion (ports of residual_set, complete_rim)
// ---------------------------------------------------------------------------
static vector<int> residual_set(const vector<int>& spokes) {
    int m = (int)spokes.size();
    int N = geodesic_count(m);
    vector<uint8_t> phi(N + 1, 0);
    for (int v : spokes)
        if (v <= N) phi[v] = 1;
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j) {
            bool adjacent = (j == i + 1) || (i == 0 && j == m - 1);
            if (adjacent) continue;
            int v = spokes[i] + spokes[j];
            if (v <= N) phi[v] = 1;
        }
    vector<int> out;
    for (int v = 1; v <= N; ++v)
        if (!phi[v]) out.push_back(v);
    return out;
}

static std::optional<vector<int>> complete_rim(const vector<int>& residual, int m) {
    if ((int)residual.size() != 2 * m) return std::nullopt;
    i64 total = 0;
    for (int v : residual) total += v;
    if (total % 3) return std::nullopt;  // sum(T) = 3 sum(b), Appendix B
    i64 rim_sum = total / 3;

    int maxv = residual.back();
    vector<uint8_t> available(maxv + 2, 0);
    for (int v : residual) available[v] = 1;
    int avail_count = (int)residual.size();
    int smallest = residual.front();  // min(T) cannot be a two-edge sum
    vector<int> rim{smallest};
    available[smallest] = 0;
    --avail_count;

    std::function<std::optional<vector<int>>(int, i64)> descend =
        [&](int k, i64 running) -> std::optional<vector<int>> {
        if (running > rim_sum) return std::nullopt;
        if (k == m) {
            int closing = rim[m - 1] + rim[0];
            if (avail_count == 1 && closing <= maxv && available[closing]) {
                if (m < 3 || rim[1] < rim[m - 1]) return rim;
            }
            return std::nullopt;
        }
        int previous = rim.back();
        for (int value = 1; value <= maxv; ++value) {
            if (!available[value]) continue;
            int link = previous + value;
            if (link > maxv || !available[link]) continue;
            available[value] = 0;
            available[link] = 0;
            avail_count -= 2;
            rim.push_back(value);
            auto found = descend(k + 1, running + value);
            rim.pop_back();
            available[value] = 1;
            available[link] = 1;
            avail_count += 2;
            if (found) return found;
        }
        return std::nullopt;
    };
    return descend(1, smallest);
}

// ---------------------------------------------------------------------------
// Independent check of a finished labeling (port of is_geodesic_leech)
// ---------------------------------------------------------------------------
static bool is_geodesic_leech(const vector<int>& spokes, const vector<int>& rims) {
    int m = (int)spokes.size();
    if ((int)rims.size() != m) return false;
    int N = geodesic_count(m);
    vector<int> weights;
    weights.reserve(N);
    for (int v : spokes) weights.push_back(v);
    for (int v : rims) weights.push_back(v);
    for (int i = 0; i < m; ++i) weights.push_back(rims[i] + rims[(i + 1) % m]);
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j) {
            bool adjacent = (j == i + 1) || (i == 0 && j == m - 1);
            if (!adjacent) weights.push_back(spokes[i] + spokes[j]);
        }
    if ((int)weights.size() != N) return false;
    std::sort(weights.begin(), weights.end());
    for (int v = 1; v <= N; ++v)
        if (weights[v - 1] != v) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Driver (ports of search / report / main)
// ---------------------------------------------------------------------------
struct Stats {
    i64 sets = 0, orders = 0, steps = 0;
    double seconds = 0.0;
};
using Labeling = std::pair<vector<int>, vector<int>>;

static std::string with_commas(i64 v) {
    std::string raw = std::to_string(v), out;
    int count = 0;
    for (int i = (int)raw.size() - 1; i >= 0; --i) {
        out += raw[i];
        if (++count == 3 && i > 0) {
            out += ',';
            count = 0;
        }
    }
    std::reverse(out.begin(), out.end());
    return out;
}

static std::optional<Labeling> search_wheel(int n, uint64_t seed, double time_limit,
                                            int selections, bool verbose, Stats& stats) {
    int m = n - 1;
    int N = geodesic_count(m);
    Rng rng(seed);
    SpokeSet state(m, N, rng.sample_1_to(N, m));

    auto started = std::chrono::steady_clock::now();
    auto elapsed = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    };
    double last_report = 0.0;

    std::optional<Labeling> solution;
    const double t_hot = 2.0, t_cold = 0.05;
    const int cycle = 20000;
    const double ratio = t_cold / t_hot;

    auto on_set = [&](SpokeSet& st) -> bool {
        stats.sets += 1;
        if (g_dump_sets) {
            vector<int> sorted_labels(st.s);
            std::sort(sorted_labels.begin(), sorted_labels.end());
            std::printf("SET");
            for (int v : sorted_labels) std::printf(" %d", v);
            std::printf("\n");
        }
        std::set<vector<int>> seen;
        bool done = false;
        cyclic_orders(st, rng, selections, [&](const vector<int>& order) {
            if (!seen.insert(order).second) return false;
            stats.orders += 1;
            auto residual = residual_set(order);
            auto rim = complete_rim(residual, m);
            if (!rim || !is_geodesic_leech(order, *rim)) return false;
            solution = Labeling(order, *rim);
            done = true;
            return true;
        });
        return done;
    };

    for (;;) {
        double now = elapsed();
        if (now >= time_limit) break;
        if (verbose && now - last_report >= 15.0) {
            last_report = now;
            std::fprintf(stderr,
                         "    W_%d: %5.0fs  %10s steps  %7lld feasible sets  %8lld orders tried"
                         "  current cost %d\n",
                         n, now, with_commas(stats.steps).c_str(), stats.sets, stats.orders,
                         state.cost);
            std::fflush(stderr);
        }
        for (int chunk = 0; chunk < 200; ++chunk) {
            stats.steps += 1;
            double temperature =
                t_hot * std::pow(ratio, (double)(stats.steps % cycle) / cycle);
            int index = rng.randrange(m);
            int value = rng.randint(1, N);
            if (state.in_set[value]) continue;
            int before = state.cost;
            int previous = state.replace(index, value);
            int rise = state.cost - before;
            if (rise > 0 && rng.random() > std::exp(-rise / temperature))
                state.replace(index, previous);
            if (state.cost == 0) {
                if (on_set(state)) {
                    stats.seconds = elapsed();
                    return solution;
                }
                index = rng.randrange(m);
                value = rng.randint(1, N);
                if (!state.in_set[value]) state.replace(index, value);
            }
        }
    }
    stats.seconds = elapsed();
    return solution;
}

static void print_list(const char* name, const vector<int>& v) {
    std::printf("      %s = [", name);
    for (size_t i = 0; i < v.size(); ++i)
        std::printf(i ? ", %d" : "%d", v[i]);
    std::printf("]\n");
}

int main(int argc, char** argv) {
    vector<int> wheels;
    uint64_t seed = 20260831;
    double time_limit = 120.0;
    int selections = 8;
    bool verbose = false;
    bool assemble_stdin = false;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "-n") {
            while (a + 1 < argc && argv[a + 1][0] != '-') wheels.push_back(std::atoi(argv[++a]));
        } else if (arg == "--seed" && a + 1 < argc) {
            seed = std::strtoull(argv[++a], nullptr, 10);
        } else if (arg == "--time-limit" && a + 1 < argc) {
            time_limit = std::atof(argv[++a]);
        } else if (arg == "--selections" && a + 1 < argc) {
            selections = std::atoi(argv[++a]);
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--assemble-stdin") {
            assemble_stdin = true;
        } else if (arg == "--dump-sets") {
            g_dump_sets = true;
        } else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 2;
        }
    }
    if (assemble_stdin) {
        // Self-test hook: each stdin line is a set of m labels (one per line);
        // prints "cost selections" so the cost function and the assembly stage
        // can be compared line by line against the Python implementation.
        std::string line;
        while (std::getline(std::cin, line)) {
            std::istringstream in(line);
            vector<int> vals;
            int x;
            while (in >> x) vals.push_back(x);
            if (vals.empty()) continue;
            int m = (int)vals.size(), N = geodesic_count(m);
            SpokeSet st(m, N, vals);
            int count = 0;
            if (st.cost == 0)
                edge_selections(st, selections, [&](const PathSet&) { ++count; return false; });
            std::printf("%d %d\n", st.cost, count);
        }
        return 0;
    }
    if (wheels.empty())
        for (int n = 7; n <= 13; ++n) wheels.push_back(n);

    std::printf("seed=%llu, time limit %gs per wheel\n\n", (unsigned long long)seed, time_limit);

    int found = 0;
    for (int n : wheels) {
        if (n < 5) {
            std::fprintf(stderr, "the wheel W_n is defined here for n >= 5\n");
            return 1;
        }
        Stats stats;
        auto result = search_wheel(n, seed + (uint64_t)n, time_limit, selections, verbose, stats);
        if (!result) {
            std::printf("W_%-3d not found   %.0fs, %lld feasible sets, %lld orders tried\n", n,
                        stats.seconds, stats.sets, stats.orders);
        } else {
            bool ok = is_geodesic_leech(result->first, result->second);
            std::printf("W_%-3d %s      %.1fs, %lld feasible sets\n", n,
                        ok ? "PASS" : "INVALID", stats.seconds, stats.sets);
            print_list("A", result->first);
            print_list("B", result->second);
            if (ok) ++found;
        }
        std::fflush(stdout);
    }
    std::printf("\n%d/%zu wheels solved and independently verified.\n", found, wheels.size());
    if (found < (int)wheels.size())
        std::printf(
            "A wheel not solved here is undecided, not shown to be impossible: "
            "this search is incomplete by design.\n");
    return 0;
}
