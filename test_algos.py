"""Own tests. Run from the project root:  python -m unittest discover -s tests -v

Non-unique answers (MST edge set, SCC labels, topological order, shortest-path
witness) are checked by VALIDITY against brute-force oracles, not by equality.
"""
import io
import itertools
import json
import os
import random
import subprocess
import sys
import time
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from graph_algos import dsu, kruskal, prim, scc, topo, dag, reconstruct_route  # noqa: E402


# ------------------------------------------------------------ oracles
def brute_mst_cost(n, edges):
    """Cost of a minimum spanning forest by exhaustive subset search (tiny m)."""
    comps = _components(n, edges, range(len(edges)))
    best = None
    for k in range(len(edges) + 1):
        for sub in itertools.combinations(range(len(edges)), k):
            if _acyclic(n, edges, sub) and _components(n, edges, sub) == comps:
                c = sum(edges[i][2] for i in sub)
                if best is None or c < best:
                    best = c
    return best, comps


def _components(n, edges, ids):
    p = list(range(n))

    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    c = n
    for i in ids:
        a, b = f(edges[i][0]), f(edges[i][1])
        if a != b:
            p[a] = b
            c -= 1
    return c


def _acyclic(n, edges, ids):
    p = list(range(n))

    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for i in ids:
        a, b = f(edges[i][0]), f(edges[i][1])
        if a == b:
            return False
        p[a] = b
    return True


def reach(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
    R = []
    for s in range(n):
        seen = {s}
        st = [s]
        while st:
            u = st.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    st.append(v)
        R.append(seen)
    return R


def bellman_ford(n, edges, s):
    INF = None
    d = [INF] * n
    d[s] = 0
    for _ in range(n):
        for u, v, w in edges:
            if d[u] is not None and (d[v] is None or d[u] + w < d[v]):
                d[v] = d[u] + w
    return d


def check_forest(tc, n, edges, res):
    ids = res["edge_ids"]
    tc.assertEqual(len(set(ids)), len(ids))
    tc.assertTrue(all(0 <= i < len(edges) for i in ids))
    tc.assertTrue(_acyclic(n, edges, ids))
    tc.assertEqual(res["components"], _components(n, edges, range(len(edges))))
    tc.assertEqual(n - len(ids), res["components"])
    tc.assertEqual(res["cost"], sum(edges[i][2] for i in ids))


def rand_wgraph(rng, n, m, lo=-5, hi=5):
    return [[rng.randrange(n), rng.randrange(n), rng.randint(lo, hi)] for _ in range(m)]


# ------------------------------------------------------------ DSU
class TestDSU(unittest.TestCase):
    def test_basic(self):
        # expected by hand: n=5; (0,1)->4 comps largest 2; (1,0) no-op; (2,3)->3,2; (0,2)->2,4; (4,4) loop no-op
        self.assertEqual(dsu(5, [[0, 1], [1, 0], [2, 3], [0, 2], [4, 4]]),
                         [[4, 2], [4, 2], [3, 2], [2, 4], [2, 4]])

    def test_empty_and_isolates(self):
        self.assertEqual(dsu(0, []), [])
        self.assertEqual(dsu(3, []), [])
        self.assertEqual(dsu(1, [[0, 0]]), [[1, 1]])

    def test_full_limit_chain(self):
        n = 100000
        ops = [[i, i + 1] for i in range(n - 1)] + [[0, n - 1]] * (250000 - n + 1)
        out = dsu(n, ops)
        self.assertEqual(out[n - 2], [1, n])
        self.assertEqual(out[-1], [1, n])
        self.assertEqual(len(out), 250000)


# ------------------------------------------------------------ MST
class TestMST(unittest.TestCase):
    def test_hand_example(self):
        # triangle 0-1(1) 1-2(2) 0-2(3) + isolate 3 + loop + parallel edge
        edges = [[0, 1, 1], [1, 2, 2], [0, 2, 3], [2, 2, -7], [0, 1, 5]]
        for f in (kruskal, prim):
            r = f(4, edges)
            self.assertEqual(r["cost"], 3)       # 1 + 2, the loop must not be used
            self.assertEqual(r["components"], 2)  # {0,1,2} and {3}
            check_forest(self, 4, edges, r)

    def test_negative_weights(self):
        edges = [[0, 1, -10 ** 9], [1, 2, -10 ** 9], [0, 2, 10 ** 9]]
        for f in (kruskal, prim):
            self.assertEqual(f(3, edges)["cost"], -2 * 10 ** 9)   # beyond 32-bit

    def test_empty(self):
        for f in (kruskal, prim):
            self.assertEqual(f(0, []), {"cost": 0, "edge_ids": [], "components": 0})
            self.assertEqual(f(3, [])["components"], 3)

    def test_random_vs_bruteforce(self):
        rng = random.Random(1)
        for _ in range(300):
            n = rng.randint(1, 6)
            edges = rand_wgraph(rng, n, rng.randint(0, 8))
            best, comps = brute_mst_cost(n, edges)
            for f in (kruskal, prim):
                r = f(n, edges)
                check_forest(self, n, edges, r)
                self.assertEqual(r["cost"], best, (f.__name__, n, edges))
                self.assertEqual(r["components"], comps)

    def test_invalid_witness_detected(self):
        edges = [[0, 1, 1], [1, 2, 1], [0, 2, 1]]
        bad = {"cost": 3, "edge_ids": [0, 1, 2], "components": 1}   # contains a cycle
        with self.assertRaises(AssertionError):
            check_forest(self, 3, edges, bad)

    def test_ties_alternative_valid(self):
        edges = [[0, 1, 1], [1, 2, 1], [0, 2, 1]]   # any 2 of 3 edges are valid
        for f in (kruskal, prim):
            r = f(3, edges)
            check_forest(self, 3, edges, r)
            self.assertEqual(r["cost"], 2)

    def test_full_limit_agree(self):
        rng = random.Random(7)
        n, m = 100000, 250000
        edges = [[rng.randrange(n), rng.randrange(n), rng.randint(-10 ** 9, 10 ** 9)]
                 for _ in range(m)]
        a, b = kruskal(n, edges), prim(n, edges)
        self.assertEqual(a["cost"], b["cost"])
        self.assertEqual(a["components"], b["components"])
        check_forest(self, n, edges, a)
        check_forest(self, n, edges, b)


# ------------------------------------------------------------ SCC
def check_scc(tc, n, edges, res):
    comp, cond = res["component"], res["condensation"]
    R = reach(n, edges)
    for u in range(n):
        for v in range(n):
            same = v in R[u] and u in R[v]
            tc.assertEqual(comp[u] == comp[v], same)
    exp = {(comp[u], comp[v]) for u, v in edges if comp[u] != comp[v]}
    tc.assertEqual({tuple(p) for p in cond}, exp)
    tc.assertEqual(len(cond), len({tuple(p) for p in cond}))     # distinct pairs


class TestSCC(unittest.TestCase):
    def test_hand(self):
        # cycle 0->1->2->0, tail 2->3, 3->4->3, parallel 2->3, loop 5->5, isolate 5
        edges = [[0, 1], [1, 2], [2, 0], [2, 3], [2, 3], [3, 4], [4, 3], [5, 5]]
        r = scc(6, edges)
        self.assertEqual(r["component"][0], r["component"][2])
        self.assertEqual(r["component"][3], r["component"][4])
        self.assertEqual(len(set(r["component"])), 3)
        self.assertEqual(len(r["condensation"]), 1)              # parallel edges merge
        check_scc(self, 6, edges, r)

    def test_empty(self):
        self.assertEqual(scc(0, []), {"component": [], "condensation": []})

    def test_random(self):
        rng = random.Random(2)
        for _ in range(300):
            n = rng.randint(1, 8)
            edges = [[rng.randrange(n), rng.randrange(n)] for _ in range(rng.randint(0, 14))]
            check_scc(self, n, edges, scc(n, edges))

    def test_deep_chain_and_cycle(self):
        n = 100000
        chain = [[i, i + 1] for i in range(n - 1)]
        r = scc(n, chain)
        self.assertEqual(len(set(r["component"])), n)
        self.assertEqual(len(r["condensation"]), n - 1)
        r2 = scc(n, chain + [[n - 1, 0]])
        self.assertEqual(len(set(r2["component"])), 1)
        self.assertEqual(r2["condensation"], [])


# ------------------------------------------------------------ topo
def valid_order(n, edges, order):
    if sorted(order) != list(range(n)):
        return False
    pos = {v: i for i, v in enumerate(order)}
    return all(pos[u] < pos[v] for u, v in edges)


class TestTopo(unittest.TestCase):
    def test_acyclic(self):
        edges = [[0, 1], [0, 2], [1, 3], [2, 3], [0, 3]]
        r = topo(5, edges)
        self.assertTrue(r["acyclic"])
        self.assertTrue(valid_order(5, edges, r["order"]))

    def test_alternative_valid_and_invalid(self):
        edges = [[0, 1], [0, 2]]
        self.assertTrue(valid_order(3, edges, [0, 1, 2]))
        self.assertTrue(valid_order(3, edges, [0, 2, 1]))      # alternative valid answer
        self.assertFalse(valid_order(3, edges, [1, 0, 2]))     # invalid witness

    def test_cycles(self):
        for n, edges in [(1, [[0, 0]]), (2, [[0, 1], [1, 0]]), (4, [[0, 1], [2, 3], [3, 2]])]:
            self.assertEqual(topo(n, edges), {"acyclic": False, "order": []})

    def test_empty(self):
        self.assertEqual(topo(0, []), {"acyclic": True, "order": []})

    def test_random_vs_reachability(self):
        rng = random.Random(3)
        for _ in range(300):
            n = rng.randint(1, 7)
            edges = [[rng.randrange(n), rng.randrange(n)] for _ in range(rng.randint(0, 9))]
            R = reach(n, edges)
            # a cycle exists iff some edge (u,v) has u reachable from v (loops: u in R[u])
            has_cycle = any(u in R[v] for u, v in edges)
            r = topo(n, edges)
            self.assertEqual(r["acyclic"], not has_cycle)
            if r["acyclic"]:
                self.assertTrue(valid_order(n, edges, r["order"]))
            else:
                self.assertEqual(r["order"], [])

    def test_deep_chain(self):
        n = 100000
        edges = [[i, i + 1] for i in range(n - 1)]
        r = topo(n, edges)
        self.assertTrue(r["acyclic"])
        self.assertEqual(r["order"], list(range(n)))
        self.assertFalse(topo(n, edges + [[n - 1, 0]])["acyclic"])


# ------------------------------------------------------------ DAG
def check_dag(tc, n, edges, s, res):
    exp = bellman_ford(n, edges, s)
    tc.assertEqual(res["distance"], exp)
    pe = res["parent_edge"]
    tc.assertEqual(pe[s], -1)
    for v in range(n):
        if v == s or exp[v] is None:
            tc.assertEqual(pe[v], -1)
        else:                                   # witness must be a tight edge
            u, x, w = edges[pe[v]]
            tc.assertEqual(x, v)
            tc.assertEqual(exp[u] + w, exp[v])


class TestDAG(unittest.TestCase):
    def test_hand(self):
        # 0->1 (4), 0->2 (1), 2->1 (2), 1->3 (-5); vertex 4 unreachable
        edges = [[0, 1, 4], [0, 2, 1], [2, 1, 2], [1, 3, -5]]
        r = dag(5, edges, 0)
        self.assertEqual(r["distance"], [0, 3, 1, -2, None])
        self.assertEqual(r["parent_edge"], [-1, 2, 1, 3, -1])
        check_dag(self, 5, edges, 0, r)

    def test_cycle_unreachable_from_source(self):
        edges = [[0, 1, 1], [2, 3, 1], [3, 2, 1]]    # cycle 2<->3 not reachable from 0
        self.assertEqual(dag(4, edges, 0), {"acyclic": False, "distance": [], "parent_edge": []})

    def test_loop(self):
        self.assertFalse(dag(1, [[0, 0, 5]], 0)["acyclic"])

    def test_single_vertex(self):
        self.assertEqual(dag(1, [], 0), {"acyclic": True, "distance": [0], "parent_edge": [-1]})

    def test_random(self):
        rng = random.Random(4)
        for _ in range(300):
            n = rng.randint(1, 8)
            perm = list(range(n))
            rng.shuffle(perm)
            edges = []
            for _ in range(rng.randint(0, 14)):
                a, b = sorted(rng.sample(range(n), 2)) if n > 1 else (0, 0)
                if a != b:
                    edges.append([perm[a], perm[b], rng.randint(-5, 5)])
            s = rng.randrange(n)
            check_dag(self, n, edges, s, dag(n, edges, s))

    def test_deep_chain_negative(self):
        n = 100000
        edges = [[i, i + 1, -10 ** 9] for i in range(n - 1)]
        r = dag(n, edges, 0)
        self.assertEqual(r["distance"][-1], -(n - 1) * 10 ** 9)   # needs > 32 bits

    def test_full_limit(self):
        rng = random.Random(5)
        n, m = 100000, 250000
        edges = []
        for _ in range(m):
            a, b = rng.randrange(n), rng.randrange(n)
            if a == b:
                b = (a + 1) % n if a + 1 < n else a - 1
            a, b = min(a, b), max(a, b)
            edges.append([a, b, rng.randint(-10 ** 9, 10 ** 9)])
        r = dag(n, edges, 0)
        self.assertTrue(r["acyclic"])
        self.assertEqual(len(r["distance"]), n)


# ------------------------------------------------------------ route helper
class TestRoute(unittest.TestCase):
    def setUp(self):
        # parallel tied edges 0->1 (ids 0 and 1, both w=2); 1->2 (id 2); 3 unreachable
        self.edges = [[0, 1, 2], [0, 1, 2], [1, 2, 3]]
        self.res = dag(4, self.edges, 0)

    def test_reachable(self):
        r = reconstruct_route(0, 2, self.res["parent_edge"], self.edges)
        self.assertEqual(r["vertices"], [0, 1, 2])
        self.assertEqual(r["edge_ids"], [0, 2])           # first tied edge kept (strict <)
        self.assertEqual(sum(self.edges[i][2] for i in r["edge_ids"]), self.res["distance"][2])

    def test_source(self):
        self.assertEqual(reconstruct_route(0, 0, self.res["parent_edge"], self.edges),
                         {"vertices": [0], "edge_ids": []})

    def test_unreachable(self):
        self.assertIsNone(reconstruct_route(0, 3, self.res["parent_edge"], self.edges))

    def test_tied_parallel(self):
        r = reconstruct_route(0, 1, self.res["parent_edge"], self.edges)
        self.assertIn(r["edge_ids"], ([0], [1]))          # both are valid shortest witnesses


# ------------------------------------------------------------ JSONL adapter
class TestAdapter(unittest.TestCase):
    def run_main(self, reqs):
        data = "\n".join(json.dumps(r) for r in reqs) + "\n"
        p = subprocess.run([sys.executable, os.path.join(ROOT, "src", "main.py")],
                           input=data, capture_output=True, text=True, timeout=120)
        self.assertEqual(p.stderr, "")
        lines = p.stdout.splitlines()
        self.assertEqual(len(lines), len(reqs))
        return [json.loads(x) for x in lines]

    def test_mixed_repeated_and_reset(self):
        reqs = [
            {"task": "ping"},
            {"task": "dsu", "n": 3, "operations": [[0, 1]]},
            {"task": "kruskal", "n": 3, "edges": [[0, 1, 1]]},
            {"task": "dsu", "n": 3, "operations": [[0, 1]]},      # identical -> identical result
            {"task": "dag", "n": 2, "source": 0, "edges": [[0, 1, 1], [1, 0, 1]]},
            {"task": "dag", "n": 2, "source": 0, "edges": [[0, 1, 1]]},   # no leak from cyclic one
            {"task": "prim", "n": 3, "edges": [[0, 1, 1]]},
            {"task": "scc", "n": 2, "edges": [[0, 1]]},
            {"task": "topo", "n": 2, "edges": [[0, 1]]},
            {"task": "ping"},
        ]
        out = self.run_main(reqs)
        self.assertEqual(out[1], out[3])
        self.assertEqual(out[0], out[9])
        self.assertFalse(out[4]["acyclic"])
        self.assertEqual(out[5]["distance"], [0, 1])
        self.assertEqual(out[2]["cost"], 1)
        self.assertEqual(out[6]["components"], 2)


if __name__ == "__main__":
    unittest.main()
