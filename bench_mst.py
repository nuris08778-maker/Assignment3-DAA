"""Kruskal vs Prim on identical graphs. Run from the project root:
    python bench/bench_mst.py

Seed 12345; graphs are built once per size OUTSIDE the timed region.
Timer: time.perf_counter(); 1 warm-up call per method per size (not counted);
5 timed runs per method per size; median reported. The timed region is exactly
one call kruskal(n, edges) / prim(n, edges) (includes adjacency building and
sorting/heap work, excludes graph generation and JSON parsing).
"""
import os
import platform
import random
import statistics
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from graph_algos import kruskal, prim  # noqa: E402

SEED = 12345
SIZES = [(1000, 5000), (10000, 50000), (100000, 250000)]
RUNS = 5


def make(n, m, rng):
    edges = [[rng.randrange(n), rng.randrange(n), rng.randint(-10 ** 9, 10 ** 9)] for _ in range(m)]
    return edges


def timed(f, n, edges):
    t = time.perf_counter()
    r = f(n, edges)
    return time.perf_counter() - t, r


def main():
    print(f"python {platform.python_version()} on {platform.platform()} ({platform.machine()})")
    print(f"seed={SEED} runs={RUNS} warmup=1 timer=time.perf_counter")
    print("| n | m | method | run1 s | run2 s | run3 s | run4 s | run5 s | median s |")
    print("|---|---|---|---|---|---|---|---|---|")
    for n, m in SIZES:
        rng = random.Random(SEED + n)
        edges = make(n, m, rng)
        costs = {}
        for name, f in (("Kruskal", kruskal), ("Prim", prim)):
            f(n, edges)                                  # warm-up
            ts = []
            for _ in range(RUNS):
                dt, r = timed(f, n, edges)
                ts.append(dt)
            costs[name] = r["cost"]
            cells = " | ".join(f"{x:.4f}" for x in ts)
            print(f"| {n} | {m} | {name} | {cells} | {statistics.median(ts):.4f} |")
        assert costs["Kruskal"] == costs["Prim"], "MST costs differ!"


if __name__ == "__main__":
    main()
