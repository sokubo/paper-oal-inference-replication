#!/usr/bin/env python3
"""make_paper_tables4.py -- render the four battery tables of paper.qmd (Table 1 = tbl-battery,
Web Appendix C tables tbl-appC1, tbl-appC1b, tbl-appC2) from analysis/output/battery4_table.csv and
battery4_selection.csv, in exactly the markdown layout used in paper.qmd.

    python3 make_paper_tables4.py [output_dir] [--splice paper.qmd]

Without --splice the four tables are printed; with --splice the table bodies in paper.qmd are
replaced in place (a backup is not made -- commit or copy first). Nothing is typed by hand: every
number comes from the CSVs, so a re-run of sim_battery4.R followed by this script refreshes the paper.
"""
import sys, re, math
import pandas as pd

out_dir = "output"
splice = None
args = sys.argv[1:]
if "--splice" in args:
    i = args.index("--splice"); splice = args[i + 1]; args = args[:i] + args[i + 2:]
if args: out_dir = args[0]

tab = pd.read_csv(f"{out_dir}/battery4_table.csv", comment="#")
sel = pd.read_csv(f"{out_dir}/battery4_selection.csv", comment="#")
W = tab.pivot_table(index=["cell", "arm"], columns="measure", values="value", aggfunc="first")
S = sel.set_index(["cell", "arm"])

MAIN = ["I", "II", "III", "M", "P1", "P2", "C", "Bd", "H", "IV6"]
CELL_LAB = {"I": "I ($n{=}500, p{=}20$)", "II": "II ($n{=}200, p{=}50$)", "III": "III (local boundary)",
            "M": "M (misspecified OM)", "P1": "P1 ($n{=}500, p{=}100$)", "P2": "P2 ($n{=}500, p{=}200$)",
            "C": "C (AR(1) $\\rho{=}.5$)", "Bd": "Bd (bounded)", "H": "H (heteroscedastic)", "IV6": "IV6 (six instruments)"}
ARM_LAB = {"oracle": "oracle set", "out_thr": "outcome-only thr.", "oads_thr": "ODS-thr", "oads_bic": "ODS-bic",
           "bch": "BCH-DS", "full": "full set", "oal_ipw": "naive OAL-IPW", "goal_ipw": "GOAL-IPW",
           "ols_oadsthr": "naive OLS on ODS-thr set", "ctmle": "C-TMLE"}

def v(cell, arm, m):
    try: return W.loc[(cell, arm), m]
    except KeyError: return float("nan")
def num(x, d, dash="–"):
    if x is None or (isinstance(x, float) and math.isnan(x)): return dash
    s = f"{x:.{d}f}"
    return s.replace("-", "−")
def pct(x, d=1):
    return num(x * 100, d)
def cov_native(cell, arm):  # practitioner pipelines: native sandwich = "emp" slot
    return v(cell, arm, "cov_stk") if arm not in ("oal_ipw", "goal_ipw", "ols_oadsthr") else v(cell, arm, "cov_emp")

# ---- Table 1 (tbl-battery)
t1 = ["| Cell | naive OAL-IPW | GOAL-IPW | C-TMLE | outcome-only thr. | ODS-thr | ODS-bic | BCH-DS |", "|---|---|---|---|---|---|---|---|"]
for c in MAIN:
    row = [CELL_LAB[c], pct(v(c, "oal_ipw", "cov_emp")), pct(v(c, "goal_ipw", "cov_emp")), pct(v(c, "ctmle", "cov_emp")), pct(v(c, "out_thr", "cov_stk")),
           "**" + pct(v(c, "oads_thr", "cov_stk")) + "**", pct(v(c, "oads_bic", "cov_stk")), pct(v(c, "bch", "cov_stk"))]
    t1.append("| " + " | ".join(row) + " |")

# ---- Web Appendix C.1 (tbl-appC1)
c1 = ["| Cell | Arm | bias | SD | SE/SD emp | SE/SD stk | cov emp | cov stk | width stk |", "|---|---|---|---|---|---|---|---|---|"]
ARMS_C1 = ["oracle", "out_thr", "oads_thr", "oads_bic", "bch", "full", "oal_ipw", "goal_ipw", "ols_oadsthr", "ctmle"]
for c in MAIN:
    for a in ARMS_C1:
        if a == "full" and c == "P2": continue
        if a == "ctmle" and math.isnan(v(c, a, "bias")): continue
        stk = a not in ("oal_ipw", "goal_ipw", "ols_oadsthr", "ctmle")
        row = [c, ARM_LAB[a], num(v(c, a, "bias"), 3), num(v(c, a, "emp_sd"), 3), num(v(c, a, "se_ratio_emp"), 2),
               num(v(c, a, "se_ratio_stk"), 2) if stk else "–", pct(v(c, a, "cov_emp")),
               pct(v(c, a, "cov_stk")) if stk else "–", num(v(c, a, "width_stk"), 2) if stk else "–"]
        c1.append("| " + " | ".join(row) + " |")

# ---- Web Appendix C.1b (tbl-appC1b)
c1b = ["| Cell | Arm | mean size | P(superset of oracle) | P(boundary retained) | P(instruments excluded) | P(OAL arm adds) |", "|---|---|---|---|---|---|---|"]
for c in MAIN:
    for a in ["out_thr", "oads_thr", "oads_bic", "bch"]:
        r = S.loc[(c, a)]
        oal = pct(r["p_oal_extra"]) if a in ("oads_thr", "oads_bic") else "–"
        row = [c, ARM_LAB[a], num(r["mean_size"], 1), pct(r["p_superset_oracle"]), pct(r["p_bdry_retained"]),
               pct(r["p_iv_excluded"]), oal]
        c1b.append("| " + " | ".join(row) + " |")

# ---- Web Appendix C.2 (tbl-appC2): boundary sequence
c2 = ["| $n$ | $\\kappa$ | $b_n$ | $b_n/\\kappa_n$ | ret. OT | ret. ODS | bias OT | bias/SD | cov OR | cov OT | cov ODS | cov BIC | cov BCH |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for k in ["1.50", "1.10", "0.75"]:
    for n in [500, 2000, 8000]:
        c = f"B{n}_k{k}"
        d = tab[tab.cell == c].design.iloc[0]
        bn = float(re.search(r"b_n=([0-9.]+)", d).group(1))
        n_tr = n * 4 / 5                      # kappa_n evaluated at the training size (K = 5)
        kap = math.log(n_tr) / math.sqrt(n_tr)
        bias = v(c, "out_thr", "bias"); sd = v(c, "out_thr", "emp_sd")
        row = [f"{n:,}", k.rstrip("0").rstrip(".") if k != "1.50" else "1.5", f"{bn:.3f}", f"{bn / kap:.2f}",
               pct(S.loc[(c, "out_thr"), "p_bdry_retained"]), pct(S.loc[(c, "oads_thr"), "p_bdry_retained"]),
               num(bias, 3), num(bias / sd, 2), pct(v(c, "oracle", "cov_stk")), pct(v(c, "out_thr", "cov_stk")),
               pct(v(c, "oads_thr", "cov_stk")), pct(v(c, "oads_bic", "cov_stk")), pct(v(c, "bch", "cov_stk"))]
        c2.append("| " + " | ".join(row) + " |")

blocks = {"tbl-battery": t1, "tbl-appC1": c1, "tbl-appC1b": c1b, "tbl-appC2": c2}
if splice is None:
    for k, b in blocks.items():
        print(f"\n### {k}\n" + "\n".join(b))
    sys.exit(0)

src = open(splice, encoding="utf-8").read()
for k, b in blocks.items():
    # the table body is the run of '|' lines immediately preceding the caption line that carries {#k}
    m = re.search(r"((?:^\|.*\n)+)\n(: [^\n]*\{#" + re.escape(k) + r"[ }])", src, flags=re.M)
    if not m: sys.exit(f"table {k} not found in {splice}")
    src = src[:m.start(1)] + "\n".join(b) + "\n" + src[m.end(1):]
open(splice, "w", encoding="utf-8").write(src)
print("spliced:", ", ".join(blocks))
