"""Standalone Branch & Bound - proves global optimum for the QAP."""
import time

N = 14
points_raw = [
    (1.3, 14.3), (1.3, 12.3), (1.3, 9.7), (4.7, 14.7),
    (4.7, 12.3), (4.7, 9.7), (10.5, 13.0), (10.5, 10.5),
    (3.0, 1.5), (1.5, 4.0), (4.5, 4.0), (10.5, 6.7),
    (13.5, 6.3), (13.5, 13.0),
]
ROLES = ["a1","a2","a3","a4","a5","a6","b1","b2","b3","c1","c2","c3","d1","d2"]
R = {n: i for i, n in enumerate(ROLES)}
PX = [int(round(p[0]*10)) for p in points_raw]
PY = [int(round(p[1]*10)) for p in points_raw]
DIST = [[abs(PX[i]-PX[j]) + abs(PY[i]-PY[j]) for j in range(N)] for i in range(N)]

OBJ_TERMS = [
    (1, R["a5"], R["c3"]), (1, R["a3"], R["b2"]), (1, R["a6"], R["b2"]),
    (1, R["a3"], R["c1"]), (3, R["a1"], R["c1"]), (2, R["a1"], R["b1"]),
    (1, R["a5"], R["d1"]), (1, R["d1"], R["b1"]), (1, R["c3"], R["b2"]),
    (2, R["a5"], R["c1"]), (1, R["a5"], R["d2"]), (2, R["a5"], R["a6"]),
    (1, R["a1"], R["b3"]), (1, R["a4"], R["a2"]), (2, R["a4"], R["c2"]),
    (2, R["c3"], R["c2"]), (1, R["d2"], R["a6"]), (2, R["a2"], R["c2"]),
    (1, R["d1"], R["b3"]), (1, R["b3"], R["b1"]), (1, R["d2"], R["b2"]),
    (1, R["a5"], R["b2"]),
]

def d1(p, q): return abs(p[0]-q[0]) + abs(p[1]-q[1])
def evaluate(a): return sum(w*d1(a[ROLES[i]], a[ROLES[j]]) for w,i,j in OBJ_TERMS)

# ── Branch order: greedy connectivity ────────────────────────────────────────
branch_order = [R["a5"], R["c1"], R["a1"], R["a6"], R["b2"],
                R["d2"], R["c3"], R["c2"], R["b1"], R["a4"],
                R["a2"], R["a3"], R["d1"], R["b3"]]
branch_pos = [0]*N
for pos, role in enumerate(branch_order):
    branch_pos[role] = pos

pair_weight = [[0]*N for _ in range(N)]
for w, ri, rj in OBJ_TERMS:
    pair_weight[ri][rj] = w
    pair_weight[rj][ri] = w

new_cost_terms = [[] for _ in range(N)]
for k in range(N):
    rk = branch_order[k]
    for j in range(k):
        rj2 = branch_order[j]
        w = pair_weight[rk][rj2]
        if w:
            new_cost_terms[k].append((w, j))

half_terms = [[] for _ in range(N+1)]
for k in range(N+1):
    for w, ri, rj in OBJ_TERMS:
        pi, pj = branch_pos[ri], branch_pos[rj]
        if pi < k <= pj:
            half_terms[k].append((w, pi))
        elif pj < k <= pi:
            half_terms[k].append((w, pj))

sorted_dists = [sorted((DIST[p][q], q) for q in range(N) if q != p) for p in range(N)]

# Tuples for speed
new_cost_terms = [tuple(t) for t in new_cost_terms]
half_terms = [tuple(t) for t in half_terms]
sorted_dists_t = [tuple(sd) for sd in sorted_dists]

# ── Upper bound from CP-SAT ─────────────────────────────────────────────────
# CP-SAT found f = 141.8 (scaled: 1418). Use UB = 1419 so B&B can independently
# find & confirm the optimal solution at 1418.
best = [1419]
best_perm = [None]
nodes = [0]

def search(depth, perm, used, cost):
    nodes[0] += 1
    if depth == N:
        if cost < best[0]:
            best[0] = cost
            best_perm[0] = perm[:]
            print(f"  *** New best: {cost} ({cost/10.0})")
        return

    lb = 0
    ht = half_terms[depth]
    for i in range(len(ht)):
        w, ad = ht[i]
        p = perm[ad]
        for dq, q in sorted_dists_t[p]:
            if not (used & (1 << q)):
                lb += w * dq
                break
    if cost + lb >= best[0]:
        return

    nct = new_cost_terms[depth]
    for pt in range(N):
        if used & (1 << pt):
            continue
        delta = 0
        for i in range(len(nct)):
            w, j = nct[i]
            delta += w * DIST[pt][perm[j]]
        new_cost = cost + delta
        if new_cost >= best[0]:
            continue

        lb2 = 0
        ht2 = half_terms[depth + 1]
        new_used = used | (1 << pt)
        for i in range(len(ht2)):
            w, ad = ht2[i]
            p2 = pt if ad == depth else perm[ad]
            for dq, q in sorted_dists_t[p2]:
                if not (new_used & (1 << q)):
                    lb2 += w * dq
                    break
        if new_cost + lb2 >= best[0]:
            continue

        perm[depth] = pt
        search(depth + 1, perm, new_used, new_cost)
    perm[depth] = -1

# ── Run ──────────────────────────────────────────────────────────────────────
perm = [-1]*N
print(f"Starting Branch & Bound  (UB = {best[0]} = {best[0]/10.})")
print(f"This will PROVE the global optimum by exhaustive search with pruning.\n")
t0 = time.time()
search(0, perm, 0, 0)
elapsed = time.time() - t0

print(f"\n{'='*60}")
print(f"B&B COMPLETE")
print(f"  Nodes explored : {nodes[0]:,}")
print(f"  Time           : {elapsed:.1f} s")
print(f"  Optimal cost   : {best[0]} (×10)  =>  f* = {best[0]/10.0}")
print(f"{'='*60}")

if best_perm[0]:
    print(f"\nPROVEN optimal assignment (f* = {best[0]/10.0}):")
    result = {}
    for depth in range(N):
        role = branch_order[depth]
        pt = best_perm[0][depth]
        result[ROLES[role]] = points_raw[pt]
        print(f"  {ROLES[role]:>2s} = {points_raw[pt]}")
    check = evaluate(result)
    print(f"\nVerification: f = {check}")
    print(f"\n*** f* = {best[0]/10.0} IS THE PROVEN GLOBAL OPTIMUM ***")
    print("(Branch & Bound explored all feasible branches; no better permutation exists.)")
else:
    print("No feasible solution found (should not happen).")
