"""
Exact global optimum solver for the QAP (Quadratic Assignment Problem)
on 14 points using OR-Tools CP-SAT.

Problem: Assign 14 fixed 2D points to 14 roles (a1..a6, b1..b3, c1..c3, d1..d2)
to minimize a weighted sum of L1 (Manhattan) distances between role pairs.
"""

from ortools.sat.python import cp_model
import time

# ── 14 fixed points ──────────────────────────────────────────────────────────
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
ROLES = ["a1", "a2", "a3", "a4", "a5", "a6",
         "b1", "b2", "b3",
         "c1", "c2", "c3",
         "d1", "d2"]
R = {name: i for i, name in enumerate(ROLES)}

# Scale coordinates ×10 to work with integers (all coords have 1 decimal)
PX = [int(round(p[0] * 10)) for p in points_raw]
PY = [int(round(p[1] * 10)) for p in points_raw]

# Precompute all-pairs L1 distances (scaled ×10)
DIST = [[abs(PX[i] - PX[j]) + abs(PY[i] - PY[j]) for j in range(N)] for i in range(N)]
FLAT_DIST = [DIST[i][j] for i in range(N) for j in range(N)]
MAX_D = max(FLAT_DIST)

# ── Objective terms: (weight, role_i_name, role_j_name) ──────────────────────
OBJ_TERMS = [
    (1, "a5", "c3"),
    (1, "a3", "b2"),
    (1, "a6", "b2"),
    (1, "a3", "c1"),
    (3, "a1", "c1"),
    (2, "a1", "b1"),
    (1, "a5", "d1"),
    (1, "d1", "b1"),
    (1, "c3", "b2"),
    (2, "a5", "c1"),
    (1, "a5", "d2"),
    (2, "a5", "a6"),
    (1, "a1", "b3"),
    (1, "a4", "a2"),
    (2, "a4", "c2"),
    (2, "c3", "c2"),
    (1, "d2", "a6"),
    (2, "a2", "c2"),
    (1, "d1", "b3"),
    (1, "b3", "b1"),
    (1, "d2", "b2"),
    (1, "a5", "b2"),
]

# ── Helper: evaluate f for a given assignment dict {role_name: point_tuple} ──
def d1(p, q):
    return abs(p[0] - q[0]) + abs(p[1] - q[1])

def evaluate(assign):
    total = 0.0
    for w, rn1, rn2 in OBJ_TERMS:
        total += w * d1(assign[rn1], assign[rn2])
    return total

# ── Verify the heuristic solution (f = 158.8) ───────────────────────────────
heuristic = {
    "a1": (1.3, 9.7),   "a2": (10.5, 6.7),  "a3": (13.5, 13.0),
    "a4": (13.5, 6.3),  "a5": (1.3, 14.3),  "a6": (4.7, 14.7),
    "b1": (1.5, 4.0),   "b2": (4.7, 12.3),  "b3": (4.5, 4.0),
    "c1": (1.3, 12.3),  "c2": (10.5, 10.5), "c3": (10.5, 13.0),
    "d1": (3.0, 1.5),   "d2": (4.7, 9.7),
}
heur_cost = evaluate(heuristic)
print(f"Heuristic solution cost: {heur_cost}")
print()

# ── Build CP-SAT model ──────────────────────────────────────────────────────
print("Building CP-SAT model...")
model = cp_model.CpModel()

# Decision variables: assign[r] ∈ {0..13} = point index for role r
assign = [model.new_int_var(0, N - 1, ROLES[r]) for r in range(N)]
model.add_all_different(assign)

# Provide the heuristic solution as an initial hint
heur_perm = []
for r in range(N):
    role_name = ROLES[r]
    pt = heuristic[role_name]
    idx = points_raw.index(pt)
    heur_perm.append(idx)
    model.add_hint(assign[r], idx)

# For each objective term, create distance variable via Element constraint
term_vars = []
for t_idx, (w, rn1, rn2) in enumerate(OBJ_TERMS):
    ri, rj = R[rn1], R[rn2]

    # combined_idx = assign[ri] * N + assign[rj]
    combined = model.new_int_var(0, N * N - 1, f"comb_{t_idx}")
    model.add(combined == assign[ri] * N + assign[rj])

    # dist_var = FLAT_DIST[combined]
    dist_var = model.new_int_var(0, MAX_D, f"dist_{t_idx}")
    model.add_element(combined, FLAT_DIST, dist_var)

    term_vars.append((w, dist_var))

# Objective: minimise ∑ w_k · dist_k  (all in ×10 scale)
model.minimize(sum(w * dv for w, dv in term_vars))

# ── Solve ────────────────────────────────────────────────────────────────────
print("Solving (proving global optimum)...\n")
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 600   # up to 10 min
solver.parameters.num_workers = 8
solver.parameters.log_search_progress = True

start = time.time()
status = solver.solve(model)
elapsed = time.time() - start

STATUS_NAME = {
    cp_model.OPTIMAL: "OPTIMAL",
    cp_model.FEASIBLE: "FEASIBLE (not proven optimal)",
    cp_model.INFEASIBLE: "INFEASIBLE",
    cp_model.MODEL_INVALID: "MODEL_INVALID",
    cp_model.UNKNOWN: "UNKNOWN",
}

print(f"\n{'='*60}")
print(f"Status : {STATUS_NAME.get(status, status)}")
print(f"Time   : {elapsed:.2f} s")

if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    obj_scaled = solver.objective_value
    obj_real = obj_scaled / 10.0
    print(f"Cost   : {obj_real}  (scaled ×10 = {obj_scaled})")

    if status == cp_model.FEASIBLE:
        bound_real = solver.best_objective_bound / 10.0
        print(f"Bound  : {bound_real}")

    print(f"\nOptimal assignment:")
    result = {}
    for r in range(N):
        p_idx = solver.value(assign[r])
        result[ROLES[r]] = points_raw[p_idx]
        print(f"  {ROLES[r]:>2s} -> point {p_idx:>2d}  {points_raw[p_idx]}")

    # Double-check with our Python evaluator
    check = evaluate(result)
    print(f"\nVerification: f = {check}")

    if status == cp_model.OPTIMAL:
        print(f"\n*** {obj_real} IS THE PROVEN GLOBAL OPTIMUM ***")
else:
    print("No solution found.")

print(f"{'='*60}")
