"""Graph algorithms for Assignment 3 (DAA). Pure Python, no recursion.

All traversals are iterative so that chains of 100000 vertices do not
overflow the interpreter stack. Every function is stateless: all working
arrays are created inside the call, so nothing leaks between requests.
"""
import heapq


# ---------------------------------------------------------------- DSU
def _find(parent, x):
    root = x
    while parent[root] != root:
        root = parent[root]
    while parent[x] != root:          # path compression
        parent[x], x = root, parent[x]
    return root


def dsu(n, operations):
    """Return [components, largest] after each union operation."""
    parent = list(range(n))
    size = [1] * n
    components = n
    largest = 1 if n else 0
    states = []
    for u, v in operations:
        ru, rv = _find(parent, u), _find(parent, v)
        if ru != rv:
            if size[ru] < size[rv]:   # union by size
                ru, rv = rv, ru
            parent[rv] = ru
            size[ru] += size[rv]
            components -= 1
            if size[ru] > largest:
                largest = size[ru]
        states.append([components, largest])
    return states


# ---------------------------------------------------------------- MST
def kruskal(n, edges):
    """Minimum spanning forest: sort edges, add if it joins two trees."""
    order = sorted(range(len(edges)), key=lambda i: edges[i][2])
    parent = list(range(n))
    size = [1] * n
    cost = 0
    ids = []
    for i in order:
        u, v, w = edges[i]
        ru, rv = _find(parent, u), _find(parent, v)
        if ru == rv:                  # loops and cycle-closing edges
            continue
        if size[ru] < size[rv]:
            ru, rv = rv, ru
        parent[rv] = ru
        size[ru] += size[rv]
        cost += w
        ids.append(i)
    return {"cost": cost, "edge_ids": ids, "components": n - len(ids)}


def prim(n, edges):
    """Minimum spanning forest: lazy Prim with a heap, restarted per component."""
    adj = [[] for _ in range(n)]
    for i, (u, v, w) in enumerate(edges):
        if u != v:
            adj[u].append((w, i, v))
            adj[v].append((w, i, u))
    visited = [False] * n
    cost = 0
    ids = []
    components = 0
    for s in range(n):
        if visited[s]:
            continue
        components += 1
        visited[s] = True
        heap = list(adj[s])
        heapq.heapify(heap)
        while heap:
            w, i, v = heapq.heappop(heap)
            if visited[v]:
                continue
            visited[v] = True
            cost += w
            ids.append(i)
            for item in adj[v]:
                if not visited[item[2]]:
                    heapq.heappush(heap, item)
    return {"cost": cost, "edge_ids": ids, "components": components}


# ---------------------------------------------------------------- SCC
def _dfs_finish_order(n, adj):
    """Iterative DFS over all vertices; returns vertices by finish time."""
    seen = [False] * n
    ptr = [0] * n
    order = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [s]
        while stack:
            u = stack[-1]
            nbrs = adj[u]
            if ptr[u] < len(nbrs):
                v = nbrs[ptr[u]]
                ptr[u] += 1
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
            else:
                stack.pop()
                order.append(u)
    return order


def scc(n, edges):
    """Kosaraju: DFS for finish order, then DFS on the transpose."""
    adj = [[] for _ in range(n)]
    radj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)
    order = _dfs_finish_order(n, adj)
    comp = [-1] * n
    c = 0
    for s in reversed(order):
        if comp[s] != -1:
            continue
        comp[s] = c
        stack = [s]
        while stack:
            u = stack.pop()
            for v in radj[u]:
                if comp[v] == -1:
                    comp[v] = c
                    stack.append(v)
        c += 1
    pairs = {(comp[u], comp[v]) for u, v in edges if comp[u] != comp[v]}
    return {"component": comp, "condensation": [list(p) for p in sorted(pairs)]}


# ---------------------------------------------------------------- topo
def topo(n, edges):
    """Three-state DFS (0 white, 1 gray, 2 black). Cycle -> acyclic False."""
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
    state = [0] * n
    ptr = [0] * n
    finish = []
    for s in range(n):
        if state[s]:
            continue
        state[s] = 1
        stack = [s]
        while stack:
            u = stack[-1]
            if ptr[u] < len(adj[u]):
                v = adj[u][ptr[u]]
                ptr[u] += 1
                if state[v] == 1:     # back edge (also a self loop)
                    return {"acyclic": False, "order": []}
                if state[v] == 0:
                    state[v] = 1
                    stack.append(v)
            else:
                state[u] = 2
                stack.pop()
                finish.append(u)
    finish.reverse()
    return {"acyclic": True, "order": finish}


# ---------------------------------------------------------------- DAG
def dag(n, edges, source):
    """Shortest paths from source in a DAG (negative weights allowed)."""
    t = topo(n, [(e[0], e[1]) for e in edges])   # whole-graph cycle check
    if not t["acyclic"]:
        return {"acyclic": False, "distance": [], "parent_edge": []}
    adj = [[] for _ in range(n)]
    for i, (u, v, w) in enumerate(edges):
        adj[u].append((v, w, i))
    dist = [None] * n
    parent_edge = [-1] * n
    dist[source] = 0
    for u in t["order"]:
        du = dist[u]
        if du is None:
            continue
        for v, w, i in adj[u]:
            nd = du + w
            if dist[v] is None or nd < dist[v]:
                dist[v] = nd
                parent_edge[v] = i
    return {"acyclic": True, "distance": dist, "parent_edge": parent_edge}


def reconstruct_route(source, target, parent_edge, edges):
    """Return {"vertices": [...], "edge_ids": [...]} or None if unreachable.

    source -> {"vertices": [source], "edge_ids": []}.
    Walks parent edges back from target; uses the stored original edge IDs,
    so parallel/tied edges are reported exactly as the DAG result chose them.
    """
    if target == source:
        return {"vertices": [source], "edge_ids": []}
    if parent_edge[target] == -1:
        return None
    vertices, ids = [target], []
    cur = target
    while cur != source:
        e = parent_edge[cur]
        if e == -1:
            return None
        ids.append(e)
        cur = edges[e][0]
        vertices.append(cur)
    vertices.reverse()
    ids.reverse()
    return {"vertices": vertices, "edge_ids": ids}
