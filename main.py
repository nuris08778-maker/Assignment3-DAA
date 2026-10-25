"""JSON Lines adapter: one request per input line, one result per output line.

Only results are written to stdout. State is never kept between requests.
NOTE: field names follow the assignment text; adjust to the template's
contract if it differs (see README in the template repo).
"""
import json
import sys

from graph_algos import dsu, kruskal, prim, scc, topo, dag


def handle(req):
    task = req.get("task")
    if task == "ping":
        return {"task": "ping", "ok": True}
    n = req.get("n", 0)
    if task == "dsu":
        return {"task": task, "states": dsu(n, req.get("operations", []))}
    if task == "kruskal":
        return {"task": task, **kruskal(n, req.get("edges", []))}
    if task == "prim":
        return {"task": task, **prim(n, req.get("edges", []))}
    if task == "scc":
        return {"task": task, **scc(n, req.get("edges", []))}
    if task == "topo":
        return {"task": task, **topo(n, req.get("edges", []))}
    if task == "dag":
        return {"task": task, **dag(n, req.get("edges", []), req["source"])}
    return {"task": task, "error": "unknown task"}


def main():
    out = sys.stdout
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            res = handle(json.loads(line))
        except Exception as exc:  # keep one result per request line
            res = {"error": str(exc)}
        out.write(json.dumps(res, separators=(",", ":")) + "\n")
        out.flush()


if __name__ == "__main__":
    main()
