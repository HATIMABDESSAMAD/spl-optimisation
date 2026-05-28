"""
Exact global optimum solver for QAP on 14 points.
Combines:
  1. CP-SAT with boolean encoding (tighter propagation)
  2. Pure-Python Branch & Bound with smart ordering
"""

from ortools.sat.python import cp_model
import time

# ── Data ─────────────────────────────────────────────────────────────────────
points_raw = [
    (1.3, 14.3),   # 0
    (1.3, 12.3),   # 1
    (1.3, 9.7),    # 2
    (4.7, 14.7),   # 3
    (4.7, 12.3),   # 4
    (4.7, 9.7),    # 5
    (10.5, 13.0),  # 6
    (10.5, 10.5),  # 7
    (3.0, 1.5),    # 8
    (1.5, 4.0),    # 9
    (4.5, 4.0),    # 10
    (10.5, 6.7),   # 11
    (13.5, 6.3),   # 12
    (13.5, 13.0),  # 13
]
N = 14
ROLES = ["a1","a2","a3","a4","a5","a6","b1","b2","b3","c1","c2","c3","d1","d2"]
R = {name: i for i, name in enumerate(ROLES)}

PX = [int(round(p[0]*10)) for p in points_raw]
PY = [int(round(p[1]*10)) for p in points_raw]
DIST = [[abs(PX[i]-PX[j]) + abs(PY[i]-PY[j]) for j in range(N)] for i in range(N)]
MAX_D = max(DIST[i][j] for i in range(N) for j in range(N))

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

def d1(p, q):
    return abs(p[0]-q[0]) + abs(p[1]-q[1])

def evaluate(asgn):
    return sum(w * d1(asgn[ROLES[ri]], asgn[ROLES[rj]]) for w, ri, rj in OBJ_TERMS)

# ── 1. CP-SAT with Boolean encoding ─────────────────────────────────────────
def solve_cpsat():
    print("=" * 60)
    print("CP-SAT solver (boolean encoding)")
    print("=" * 60)
    model = cp_model.CpModel()

    # Boolean assignment variables x[r][p]
    x = [[model.new_bool_var(f'x_{r}_{p}') for p in range(N)] for r in range(N)]
    for r in range(N):
        model.add_exactly_one(x[r])
    for p in range(N):
        model.add_exactly_one([x[r][p] for r in range(N)])

    # For each objective term, create distance via conditional constraints
    term_vars = []
    for t_idx, (w, ri, rj) in enumerate(OBJ_TERMS):
        dist_var = model.new_int_var(0, MAX_D, f'd_{t_idx}')
        for p in range(N):
            expr = sum(DIST[p][q] * x[rj][q] for q in range(N))
            model.add(dist_var == expr).only_enforce_if(x[ri][p])
        term_vars.append((w, dist_var))

    model.minimize(sum(w * dv for w, dv in term_vars))

    # Hint from best known solution (141.8 from previous CP-SAT run)
    cpsat_map = {0:7, 1:8, 2:13, 3:10, 4:4, 5:3, 6:11, 7:1, 8:12, 9:6, 10:9, 11:2, 12:5, 13:0}
    for r in range(N):
        for p in range(N):
            model.add_hint(x[r][p], 1 if cpsat_map[r] == p else 0)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 600
    solver.parameters.num_workers = 8
    solver.parameters.log_search_progress = True

    t0 = time.time()
    status = solver.solve(model)
    elapsed = time.time() - t0

    names = {cp_model.OPTIMAL: "OPTIMAL", cp_model.FEASIBLE: "FEASIBLE",
             cp_model.INFEASIBLE: "INFEASIBLE", cp_model.UNKNOWN: "UNKNOWN"}
    print(f"\nStatus: {names.get(status, status)}  Time: {elapsed:.1f}s")

    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        obj = solver.objective_value / 10.0
        print(f"Cost: {obj}")
        if status == cp_model.FEASIBLE:
            print(f"Bound: {solver.best_objective_bound / 10.0}")
        result = {}
        for r in range(N):
            for p in range(N):
                if solver.value(x[r][p]):
                    result[ROLES[r]] = points_raw[p]
                    print(f"  {ROLES[r]:>2s} -> pt {p:>2d}  {points_raw[p]}")
        check = evaluate(result)
        print(f"Verification: {check}")
        if status == cp_model.OPTIMAL:
            print(f"\n*** {obj} IS THE PROVEN GLOBAL OPTIMUM ***")
        return obj, result, (status == cp_model.OPTIMAL)
    return None, None, False


# ── 2. Branch & Bound ────────────────────────────────────────────────────────
def solve_bb(upper_bound_scaled):
    print("\n" + "=" * 60)
    print(f"Branch & Bound (upper bound = {upper_bound_scaled/10.0})")
    print("=" * 60)

    # Branching order: greedy connectivity schedule
    branch_order = [R["a5"], R["c1"], R["a1"], R["a6"], R["b2"],
                    R["d2"], R["c3"], R["c2"], R["b1"], R["a4"],
                    R["a2"], R["a3"], R["d1"], R["b3"]]
    branch_pos = [0]*N
    for pos, role in enumerate(branch_order):
        branch_pos[role] = pos

    # pair_weight[ri][rj] = weight of term connecting roles ri, rj (0 if none)
    pair_weight = [[0]*N for _ in range(N)]
    for w, ri, rj in OBJ_TERMS:
        pair_weight[ri][rj] = w
        pair_weight[rj][ri] = w

    # For each depth k: fully-constrained terms formed at that depth
    # new_cost_terms[k] = [(weight, earlier_depth), ...]
    new_cost_terms = [[] for _ in range(N)]
    for k in range(N):
        rk = branch_order[k]
        for j in range(k):
            rj = branch_order[j]
            w = pair_weight[rk][rj]
            if w:
                new_cost_terms[k].append((w, j))

    # Half-constrained terms at each depth
    # half_terms[k] = [(weight, assigned_depth_index), ...]
    half_terms = [[] for _ in range(N+1)]
    for k in range(N+1):
        for w, ri, rj in OBJ_TERMS:
            pi, pj = branch_pos[ri], branch_pos[rj]
            if pi < k <= pj:
                half_terms[k].append((w, pi))
            elif pj < k <= pi:
                half_terms[k].append((w, pj))

    # Sorted distances for fast min-over-remaining lookup
    sorted_dists = [sorted((DIST[p][q], q) for q in range(N) if q != p) for p in range(N)]

    # Convert to tuples for speed
    new_cost_terms = [tuple(t) for t in new_cost_terms]
    half_terms = [tuple(t) for t in half_terms]
    sorted_dists_t = [tuple(sd) for sd in sorted_dists]
    DIST_flat = DIST  # use direct indexing

    best = [int(upper_bound_scaled)]
    best_perm = [None]
    nodes = [0]
    last_print = [time.time()]

    def search(depth, perm, used, cost):
        nodes[0] += 1

        if nodes[0] % 5_000_000 == 0:
            now = time.time()
            rate = nodes[0] / (now - last_print[0] + 1e-9) if last_print[0] else 0
            print(f"  ... {nodes[0]:>12,} nodes  best={best[0]/10.0}  "
                  f"depth={depth}  rate={rate/1e6:.1f}M/s")

        if depth == N:
            if cost < best[0]:
                best[0] = cost
                best_perm[0] = perm[:]
                elapsed = time.time() - t0
                print(f"  *** New best: {cost} ({cost/10.0}) at node {nodes[0]:,} "
                      f"[{elapsed:.1f}s]")
            return

        # Lower bound: min-over-remaining for each half-constrained term
        lb = 0
        ht = half_terms[depth]
        for i in range(len(ht)):
            w, ad = ht[i]
            p = perm[ad]
            sd = sorted_dists_t[p]
            for dq, q in sd:
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
                delta += w * DIST_flat[pt][perm[j]]

            new_cost = cost + delta
            if new_cost >= best[0]:
                continue

            # Tighter lb after this assignment
            lb2 = 0
            ht2 = half_terms[depth + 1]
            new_used = used | (1 << pt)
            for i in range(len(ht2)):
                w, ad = ht2[i]
                p2 = pt if ad == depth else perm[ad]
                sd = sorted_dists_t[p2]
                for dq, q in sd:
                    if not (new_used & (1 << q)):
                        lb2 += w * dq
                        break

            if new_cost + lb2 >= best[0]:
                continue

            perm[depth] = pt
            search(depth + 1, perm, new_used, new_cost)

        perm[depth] = -1

    perm = [-1]*N
    t0 = time.time()
    search(0, perm, 0, 0)
    elapsed = time.time() - t0

    print(f"\nB&B complete: {elapsed:.1f}s, {nodes[0]:,} nodes")
    print(f"Optimal cost (×10): {best[0]}  =>  f* = {best[0]/10.0}")

    if best_perm[0]:
        print("Optimal assignment:")
        result = {}
        for depth in range(N):
            role = branch_order[depth]
            pt = best_perm[0][depth]
            result[ROLES[role]] = points_raw[pt]
            print(f"  {ROLES[role]:>2s} -> pt {pt:>2d}  {points_raw[pt]}")
        check = evaluate(result)
        print(f"Verification: f = {check}")
        return best[0]/10.0, result
    return None, None


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # First verify heuristic
    heuristic = {
        "a1":(1.3,9.7), "a2":(10.5,6.7), "a3":(13.5,13), "a4":(13.5,6.3),
        "a5":(1.3,14.3), "a6":(4.7,14.7), "b1":(1.5,4), "b2":(4.7,12.3),
        "b3":(4.5,4), "c1":(1.3,12.3), "c2":(10.5,10.5), "c3":(10.5,13),
        "d1":(3,1.5), "d2":(4.7,9.7),
    }
    print(f"Heuristic cost: {evaluate(heuristic)}\n")

    # Run CP-SAT first
    cpsat_cost, cpsat_result, cpsat_optimal = solve_cpsat()

    if cpsat_optimal:
        print(f"\nDone. Global optimum = {cpsat_cost}")
    else:
        # Use CP-SAT result as upper bound for B&B
        ub = int(round(cpsat_cost * 10)) if cpsat_cost else 1588
        print(f"\nCP-SAT did not prove optimality. Running B&B with UB={ub/10.0}...")
        bb_cost, bb_result = solve_bb(ub)
        if bb_cost is not None:
            print(f"\n*** PROVEN GLOBAL OPTIMUM: f* = {bb_cost} ***")
