#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXTENSÕES para os dois papers desmembrados.
  Paper A (fiscal/Laffer): sensibilidade do efeito Laffer ao parâmetro-chave phi
     (elasticidade do capital humano ao gasto) — % erodido, break-even e custo
     fiscal líquido em função de phi.
  Paper B (desigualdade regional): TRAJETÓRIA temporal do Gini salarial entre UF
     e DECOMPOSIÇÃO de Theil (entre vs. dentro de macrorregiões).

Reimplementa um solver compacto do DSGE (parametrizado por phi/tau) para as
sensibilidades, e reutiliza os CSVs do bloco empírico para a desigualdade.
Requer que 01_ e 02_ tenham sido executados antes (para os CSVs em out/).
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import root
import matplotlib.pyplot as plt

import _common as C

C.set_style()

# ---------------------------------------------------------------------------
# Solver compacto do DSGE, parametrizado por phi e tau
# ---------------------------------------------------------------------------
BASE = dict(alpha=0.40, beta=0.96, dk=0.10, dh=0.04, tau=0.33, theta=1.75,
            phi=0.30, gsh=0.05)
TT = 220


def build(phi, tau):
    p = dict(BASE); p["phi"] = phi; p["tau"] = tau
    a, b, dk, dh, th, gsh = (p["alpha"], p["beta"], p["dk"], p["dh"],
                             p["theta"], p["gsh"])
    Rtil = (1 / b - 1 + dk) / (1 - tau); KY = a / Rtil
    B = (a / Rtil) ** (a / (1 - a)); cshare = 1 - dk * KY - gsh
    M = (1 - tau) * (1 - a) / cshare; l = M / (th + M)
    Y = B * l; K = KY * Y; Ge = gsh * Y
    A = dh / (Ge ** phi)
    ss = dict(H=1, K=K, l=l, Y=Y, Ge=Ge, C=cshare * Y, T=tau * Y, S=(1 - a) * Y / l)
    return p, a, b, dk, dh, th, KY, B, A, ss


def solve_cut(phi, tau, cut):
    p, a, b, dk, dh, th, KY, B, A, ss = build(phi, tau)
    Ge = np.full(TT + 1, ss["Ge"] * (1 - cut))
    H = np.empty(TT + 2); H[0] = 1
    for t in range(TT + 1):
        H[t + 1] = (1 - dh) * H[t] + A * (Ge[min(t, TT)] ** phi)
    # EE terminal (nível fixo de gasto)
    Gt = Ge[TT]; Ht = A * Gt ** phi / dh
    lt = ((1 - tau) * (1 - a) * B * Ht + th * Gt) / (
        B * Ht * (th * (1 - dk * KY) + (1 - tau) * (1 - a)))
    Ct = B * Ht * lt * (1 - dk * KY) - Gt; Kt = KY * B * Ht * lt
    Kss = ss["K"]

    def Yf(K, Hh, l):
        return K ** a * (Hh * l) ** (1 - a)

    def res(x):
        Cc = x[:TT]; l = x[TT:2 * TT]; Kn = x[2 * TT:]
        K = np.empty(TT + 1); K[0] = Kss; K[1:] = Kn
        Y = Yf(K[:TT], H[:TT], l)
        r1 = th / (1 - l) - (1 - tau) * (1 - a) * Y / (Cc * l)
        r2 = Cc + Kn - (1 - dk) * K[:TT] + Ge[:TT] - Y
        Yn = Yf(Kn, H[1:TT + 1], np.append(l[1:], lt))
        r3 = 1.0 / Cc - b / np.append(Cc[1:], Ct) * ((1 - tau) * a * Yn / Kn + 1 - dk)
        return np.concatenate([r1, r2, r3])

    x0 = np.concatenate([np.full(TT, Ct), np.full(TT, lt), np.full(TT, Kt)])
    sol = root(res, x0, method="lm", options={"maxiter": 8000, "xtol": 1e-11})
    Cc = sol.x[:TT]; l = sol.x[TT:2 * TT]; Kn = sol.x[2 * TT:]
    K = np.append(Kss, Kn); Y = Yf(K[:TT], H[:TT], l)
    S = (1 - a) * Y / l; Tax = tau * Y
    return dict(Y=Y, S=S, Tax=Tax, Ge=Ge[:TT], H=H[:TT], ss=ss,
                S_dev=100 * (S / ss["S"] - 1))


# ===========================================================================
# PAPER A — Sensibilidade do efeito Laffer a phi
# ===========================================================================
PIB = C.PIB_BR_BI
phis = np.array([0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
cutA = 0.10
rowsA = []
for phi in phis:
    r = solve_cut(phi, 0.33, cutA)
    kappa = PIB / r["ss"]["Y"]
    dGe = (r["Ge"] - r["ss"]["Ge"]) * kappa
    dTax = (r["Tax"] - r["ss"]["T"]) * kappa
    sav = -dGe; los = -dTax
    w30 = C.DISC ** np.arange(30)
    pv_s = np.sum(w30 * sav[:30]); pv_l = np.sum(w30 * los[:30])
    cum_s = np.cumsum(sav); cum_l = np.cumsum(los)
    be = np.where(cum_l > cum_s)[0]
    rowsA.append(dict(phi=phi, Y_LP=100 * (r["Y"][-1] / r["ss"]["Y"] - 1),
                      eroded30=100 * pv_l / pv_s, net30=pv_s - pv_l,
                      breakeven=int(be[0]) if len(be) else np.nan))
dfA = pd.DataFrame(rowsA)
print("Paper A — sensibilidade a phi (corte 10%):")
print(dfA.to_string(index=False))

fig, ax1 = plt.subplots(figsize=(8, 4.6))
ax1.plot(dfA["phi"], dfA["eroded30"], "o-", color="#b2182b", lw=2.2,
         label="% da economia erodido em 30 anos")
ax1.set_xlabel("$\\phi$ — elasticidade do capital humano ao gasto")
ax1.set_ylabel("% erodido (30 anos)", color="#b2182b")
ax1.tick_params(axis="y", labelcolor="#b2182b")
ax1.axhline(100, color="gray", ls=":", lw=1)
ax2 = ax1.twinx()
ax2.plot(dfA["phi"], dfA["breakeven"], "s--", color="#2166ac", lw=2,
         label="Ano de break-even fiscal")
ax2.set_ylabel("Ano de break-even fiscal", color="#2166ac")
ax2.tick_params(axis="y", labelcolor="#2166ac")
ax2.grid(False)
ax1.set_title("Efeito Laffer da educação: sensibilidade ao parâmetro $\\phi$ (corte de 10%)")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "laffer_phi.png")); plt.close()

with open(os.path.join(C.TAB, "laffer_phi.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
    f.write("$\\phi$ & $\\Delta Y_{LP}$ (\\%) & \\% erodido 30a & "
            "Break-even (ano) & Líquido 30a (R\\$ bi) \\\\\n\\midrule\n")
    for _, r in dfA.iterrows():
        be = "--" if np.isnan(r["breakeven"]) else f"{int(r['breakeven'])}"
        f.write(f"{C.num(r['phi'],2)} & {C.num(r['Y_LP'])} & "
                f"{C.num(r['eroded30'],1)}\\% & {be} & {C.brl(r['net30'])} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ===========================================================================
# PAPER B — Trajetória do Gini e decomposição de Theil
# ===========================================================================
el = pd.read_csv(os.path.join(C.OUT, "elasticidades_uf.csv")).set_index("uf")
painel = pd.read_csv(os.path.join(C.OUT, "painel_calibrado.csv"))
base_wage = painel[painel["ano"] == painel["ano"].max()].set_index("uf")["salario_medio"]
el["base_wage"] = base_wage
el["pop"] = pd.Series(C.POP_UF)
el["regiao"] = pd.Series(C.REGIAO)
norm = el["elast_eb"].mean()

# resposta salarial agregada por corte (via solver compacto)
paths = {cut: solve_cut(0.30, 0.33, cut)["S_dev"] for cut in [0.05, 0.10, 0.25]}


def gini_at(pct):
    w_after = el["base_wage"].values * (1 + pct)
    return C.gini(w_after, el["pop"].values)


fig, ax = plt.subplots(figsize=(8, 4.6))
horiz = 30
g0 = C.gini(el["base_wage"].values, el["pop"].values)
for cut, col in zip([0.05, 0.10, 0.25], ["#ffb14e", "#ea5f94", "#9d02d7"]):
    S = paths[cut]
    gtraj = [g0] + [gini_at(S[t] / 100 * (el["elast_eb"].values / norm))
                    for t in range(horiz)]
    ax.plot(range(-1, horiz), gtraj, color=col, lw=2.2, label=f"Corte {int(cut*100)}%")
ax.axhline(g0, color="black", ls=":", lw=1)
ax.text(horiz - 8, g0, "Gini inicial", fontsize=8, va="bottom")
ax.set_xlabel("Anos após o corte"); ax.set_ylabel("Gini salarial entre UF")
ax.set_title("Trajetória temporal da desigualdade salarial regional após o corte")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "gini_trajetoria.png")); plt.close()


def theil(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    p = w / w.sum(); mu = np.sum(p * x)
    return float(np.sum(p * (x / mu) * np.log(x / mu)))


def theil_decomp(wage):
    w = el["pop"].values; reg = el["regiao"].values
    Ttot = theil(wage, w)
    p = w / w.sum(); mu = np.sum(p * wage)
    between = 0.0; within = 0.0
    for R in np.unique(reg):
        m = reg == R; wr = w[m]; pr = wr.sum() / w.sum()
        mur = np.sum((wr / wr.sum()) * wage[m])
        between += pr * (mur / mu) * np.log(mur / mu)
        within += pr * (mur / mu) * theil(wage[m], wr)
    return Ttot, between, within


w_before = el["base_wage"].values
S25 = paths[0.25]
w_after = w_before * (1 + S25[-1] / 100 * (el["elast_eb"].values / norm))
Tb = theil_decomp(w_before); Ta = theil_decomp(w_after)
print(f"\nPaper B — Theil antes: total={Tb[0]:.4f} (entre {Tb[1]:.4f}, dentro {Tb[2]:.4f})")
print(f"Paper B — Theil depois (corte 25%): total={Ta[0]:.4f} "
      f"(entre {Ta[1]:.4f}, dentro {Ta[2]:.4f})")

with open(os.path.join(C.TAB, "theil_decomp.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrr}\n\\toprule\n")
    f.write("Cenário & Theil total & Entre regiões & Dentro das regiões \\\\\n\\midrule\n")
    f.write(f"Antes do corte & {C.num(Tb[0],4)} & {C.num(Tb[1],4)} & {C.num(Tb[2],4)} \\\\\n")
    f.write(f"Após corte de 25\\% & {C.num(Ta[0],4)} & {C.num(Ta[1],4)} & {C.num(Ta[2],4)} \\\\\n")
    f.write(f"Variação (\\%) & {C.num(100*(Ta[0]/Tb[0]-1),1)} & "
            f"{C.num(100*(Ta[1]/Tb[1]-1),1)} & {C.num(100*(Ta[2]/Tb[2]-1),1)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

print("\nExtensões concluídas: laffer_phi.(png|tex), gini_trajetoria.png, theil_decomp.tex")
