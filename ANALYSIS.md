# Assignment 3 – Graph Algorithms: Analysis

Language: Python 3.12+ (developed and measured on 3.13). Source: `src/graph_algos.py` (algorithms), `src/main.py` (JSON Lines adapter).

## Build / run / test commands

```bash
# no build step (pure Python, standard library only)
python src/main.py < requests.jsonl > results.jsonl     # JSON Lines adapter
python -m unittest discover -s tests -v                 # own tests (32 tests, ~5 s)
python bench/bench_mst.py                               # Kruskal vs Prim benchmark
python tools/run.py --setup && python tools/run.py --check   # template public checks (36 cases)
```

## State reset and deep DFS

* **State reset.** Every algorithm is a plain function whose working arrays (`parent`, `size`, `adj`, `visited`, `state`, ...) are created inside the call. There are no module-level or global mutable variables, so request *k* cannot influence request *k+1*. `main.py` handles each line independently and `ping` is answered without touching any state. Test `TestAdapter.test_mixed_repeated_and_reset` sends mixed and repeated requests (including a cyclic DAG request followed by an acyclic one) and checks identical requests give identical results.
* **Deep DFS.** No function is recursive. DFS in `topo` and the first pass of `scc` keeps an explicit stack of vertices plus a per-vertex edge pointer `ptr[u]` (the iterative equivalent of "resume the loop over neighbours"). The second Kosaraju pass is a plain stack traversal. A chain of 100 000 vertices therefore uses O(n) heap memory and no interpreter stack (tests `test_deep_chain*`).

## Methods

Notation: n vertices, m edges (or union operations). Memory counts the input as given.

### dsu
* **Steps.** `parent[i]=i`, `size[i]=1`, `components=n`, `largest=1` (0 if n=0). For each pair `[u,v]` find both roots (iterative, two-pass path compression); if different, attach the smaller tree under the larger, add sizes, decrement `components`, update `largest`; append `[components, largest]` (also for no-op unions).
* **Invariants.** `parent` forms a forest; `size[r]` is the exact size of the tree rooted at `r`; `components` = number of roots; `largest` = max over root sizes (sizes only grow, so a running max is exact).
* **Correctness.** Union of two different sets merges exactly two components; same-set unions (including loops `[u,u]`) change nothing.
* **Boundaries.** n=0, empty operations, loops, repeated pairs, n=1.
* **Complexity.** O((n+m)·α(n)) time, O(n) extra memory (plus O(m) output).

### kruskal
* **Steps.** Sort edge IDs by weight (stable); scan with a DSU, taking an edge when its endpoints are in different trees; `components = n − |taken|`.
* **Invariants.** Taken edges form a forest; after processing a prefix of the sorted order the DSU equals the connectivity of all edges seen so far.
* **Correctness.** Cut/cycle property: the lightest edge crossing any cut that respects the current forest is in some MST, so the greedy choice is safe; connectivity is preserved because an edge is skipped only if its endpoints are already connected. Negative weights need no special care.
* **Boundaries.** Loops (same root → skipped), parallel edges, isolates (counted in `components`), n=0, weights ±10⁹ (Python ints, no overflow; test with sum −2·10⁹).
* **Complexity.** O(m log m) time, O(n + m) memory.

### prim
* **Steps.** Adjacency lists `(w, id, v)` (loops dropped). For each unvisited start vertex: new component, mark visited, push its edges to a binary heap; pop the lightest, skip if the far end is visited, else take it and push the new vertex's edges to unvisited neighbours. Restart in every component, isolates become one-vertex components.
* **Invariants.** The visited set of the current component is the vertex set of a tree; the heap contains (at least) every edge leaving it, stale entries are discarded on pop.
* **Correctness.** The popped valid edge is the lightest crossing the cut (visited, rest) – safe by the cut property. Restarting handles disconnected graphs, giving a spanning *forest*.
* **Boundaries.** As Kruskal; parallel edges are both pushed and the lighter pops first.
* **Complexity.** O(m log m) time (lazy deletion, heap holds ≤ 2m entries), O(n + m) memory.

### scc (Kosaraju)
* **Steps.** (1) Build `adj` and transpose `radj`. (2) Iterative DFS over all vertices on `adj`, record vertices in finish order. (3) Process vertices in reverse finish order; each unlabeled vertex starts a traversal over `radj` that labels one component. (4) `condensation` = sorted set of distinct `(comp[u], comp[v])` with different labels.
* **Invariants.** In pass 2, a traversal from `s` over `radj` can only reach vertices that are in no earlier component, because the vertex with the latest finish time lies in a source component of the condensation.
* **Correctness.** Standard Kosaraju argument: the component containing the latest-finishing vertex is a source in the condensation, so in the transpose it is a sink and the traversal stays inside it; remove it and repeat.
* **Boundaries.** n=0, loops (same component, no condensation pair), parallel edges (merged into one pair), isolates (singleton components), 100 000-vertex chain and ring.
* **Complexity.** O(n + m) time (plus O(c log c) to sort pairs), O(n + m) memory.

### topo
* **Steps.** Three-state DFS (white 0 / gray 1 / black 2), iterative. An edge to a gray vertex is a back edge → return `acyclic=false, order=[]`. Otherwise vertices are appended on turning black; the reversed list is the order.
* **Invariants.** Gray vertices are exactly the current DFS stack; black vertices have all descendants black and already in `finish`.
* **Correctness.** A directed graph has a cycle iff DFS finds a back edge; in a DAG every edge `u→v` finishes `v` before `u`, so reverse finish order is topological. Self loops are back edges (`u` is gray when its own edge is examined).
* **Boundaries.** n=0, loops, parallel edges, cycles in any component, 100 000-vertex chain.
* **Complexity.** O(n + m) time and memory.

### dag (shortest paths)
* **Steps.** Run `topo` on **all** edges first (cycle anywhere → `acyclic=false`, `distance=[]`, `parent_edge=[]`). Otherwise `dist[source]=0`, others `None`; walk vertices in topological order, skip those with `dist=None`, relax outgoing edges with strict `<`, storing the original edge ID in `parent_edge`.
* **Invariants.** When `u` is processed, `dist[u]` is final, because all its predecessors precede it in the order. `parent_edge[v]` is an edge `(u,v)` with `dist[u]+w = dist[v]`.
* **Correctness.** Induction along the topological order; negative weights are fine because no vertex is relaxed after being processed. Unreachable vertices stay `None`/`-1`; the source has `0`/`-1`.
* **Boundaries.** n=1, loop (cycle), a cycle unreachable from the source (still reported as cyclic – test `test_cycle_unreachable_from_source`), parallel/tied edges (first one kept), weights ±10⁹ over 100 000 edges (−9.9999·10¹³, exceeds 32 bits).
* **Complexity.** O(n + m) time and memory.

### Route helper `reconstruct_route(source, target, parent_edge, edges)`
Walks `parent_edge` from `target` back to `source`, collecting original edge IDs and vertices, then reverses. Returns `{"vertices": [...], "edge_ids": [...]}`, or `None` when the target is unreachable. O(path length).

Example calls (from `tests/test_algos.py::TestRoute`), graph `edges=[[0,1,2],[0,1,2],[1,2,3]]`, n=4, source 0:

| Call | Result | Why |
|---|---|---|
| `reconstruct_route(0, 2, pe, edges)` (reachable) | `vertices [0,1,2]`, `edge_ids [0,2]`, total weight 5 = `distance[2]` | path via parent edges 2 then 0 |
| `reconstruct_route(0, 0, pe, edges)` (source) | `vertices [0]`, `edge_ids []` | zero-length route |
| `reconstruct_route(0, 3, pe, edges)` (unreachable) | `None` | `parent_edge[3] = -1` |
| `reconstruct_route(0, 1, pe, edges)` (tied/parallel) | `edge_ids [0]` | edges 0 and 1 tie; strict `<` keeps the first. `[1]` would be equally valid |

## Testing summary (`tests/test_algos.py`, 32 tests, all pass)

Inputs, hand-derived expected answers (with the reasoning in comments) and oracles are in the file; run with the command above. Non-unique answers are checked **by validity**, not equality:

* MST: acyclic, ids distinct and in range, `cost = Σ w`, component count equals that of the whole graph, cost equals exhaustive-search optimum on 300 random tiny graphs; the invalid witness (a cycle) is rejected (`test_invalid_witness_detected`); alternative valid answer: tied triangle.
* SCC: labels compared with the mutual-reachability relation; condensation pairs compared with the edge set (300 random graphs).
* Topo: positions respect every edge; alternative valid order `[0,2,1]` accepted, invalid `[1,0,2]` rejected; cyclic flag compared with reachability.
* DAG: distances against Bellman–Ford; `parent_edge[v]` must be a tight edge into `v`.
* Limits: n=100 000 with m=250 000 for DSU, Kruskal/Prim (also cross-checked: same cost and components), DAG; deep chains/rings of 100 000 vertices for SCC, topo, DAG.
* Adapter: mixed and repeated JSON Lines requests via a subprocess, one result per line, empty stderr, no state leak.

## Benchmark: Kruskal vs Prim

* **Command:** `python bench/bench_mst.py` (raw output saved in `bench/results.txt`).
* **Inputs:** random multigraphs (loops and parallel edges possible), weights uniform in [−10⁹, 10⁹], generated with `random.Random(12345 + n)`; the *same* edge list is given to both methods at each size. Sizes (n, m): (1 000, 5 000), (10 000, 50 000), (100 000, 250 000).
* **Environment:** CPython 3.13.17, Linux 6.18 x86_64, Intel Xeon @ 2.10 GHz, 2 vCPU (shared cloud sandbox, so absolute numbers are noisy).
* **Timer / warm-up / boundary:** `time.perf_counter()`; 1 untimed warm-up call per method and size; 5 timed runs; median reported. The timed region is one call `kruskal(n, edges)` or `prim(n, edges)` (adjacency construction/sorting included); graph generation and JSON parsing are excluded. No compilation exists in Python.

| n | m | Kruskal runs (s) | Kruskal median | Prim runs (s) | Prim median |
|---|---|---|---|---|---|
| 1 000 | 5 000 | 0.0032, 0.0028, 0.0027, 0.0027, 0.0026 | **0.0027** | 0.0041, 0.0043, 0.0046, 0.0045, 0.0050 | **0.0045** |
| 10 000 | 50 000 | 0.0566, 0.0564, 0.0608, 0.0725, 0.0496 | **0.0566** | 0.1323, 0.1692, 0.1106, 0.1121, 0.1196 | **0.1196** |
| 100 000 | 250 000 | 0.5487, 0.5263, 0.5000, 0.4745, 0.4856 | **0.5000** | 1.4482, 1.4011, 1.2784, 1.3429, 1.3077 | **1.3429** |

**Observation.** Both are O(m log m) and returned the same cost at every size (asserted in the script), but Kruskal was about 1.7–2.7× faster here. In CPython the sort uses a C implementation and the DSU loop is simple, whereas Prim builds 2m adjacency tuples and performs up to 2m Python-level `heappush`/`heappop` calls. On sparse graphs Kruskal's single sort is cheaper; Prim is expected to catch up on dense graphs, which were not measured. This is a single-machine measurement, so only the ranking, not the exact ratios, should be trusted.

## Sources and AI assistance

* Algorithms follow standard textbook descriptions (Cormen et al., *Introduction to Algorithms*, chapters on disjoint sets, MST, DFS/SCC, shortest paths in DAGs) and the course lectures.
* **AI assistance:** the first version of the code, the tests, the benchmark script and this document were drafted with Claude (Anthropic). I reviewed and ran them myself. **<-- Edit this paragraph so it matches the syllabus rules exactly (what was generated, what you changed, what you verified).**
