#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLOCO DSGE — Modelo neoclássico com capital físico e humano, tributação
distorciva e gasto público em educação como insumo produtivo da acumulação
de capital humano (Uzawa-Lucas / Glomm-Ravikumar com governo).
Solução de transição por PREVISÃO PERFEITA NÃO-LINEAR (stacked-time Newton).

IMPACTOS QUANTIFICADOS (versão ampliada):
  (1) Cenários de corte 1/5/10/25% + teto de gasto (EC 95/2016).
  (2) Arrecadação em R$ bilhões e CUSTO FISCAL LÍQUIDO (Laffer da educação).
  (3) DESIGUALDADE REGIONAL: Gini/Lorenz e razão salarial Sudeste-Nordeste.
  (4) BEM-ESTAR (equivalente de consumo), MULTIPLICADORES e trajetórias 30 anos.
"""
import os
import json

import numpy as np
import pandas as pd
from scipy.optimize import root
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import _common as C

C.set_style()

# ===========================================================================
# PARÂMETROS (calibração anual, Brasil)
# ===========================================================================
P = dict(
    alpha=0.40,   # participação do capital físico
    beta=0.96,    # desconto intertemporal (anual)
    dk=0.10,      # depreciação do capital físico
    dh=0.04,      # depreciação do capital humano
    tau=0.33,     # carga tributária efetiva (~33% do PIB)
    theta=1.75,   # peso do lazer na utilidade
    phi=0.30,     # elasticidade do capital humano ao gasto público
    gsh=0.05,     # gasto público em educação / PIB no EE
    rhoG=0.90,    # persistência do choque fiscal (cenário transitório)
)
a, b, dk, dh = P["alpha"], P["beta"], P["dk"], P["dh"]
tau, th, phi = P["tau"], P["theta"], P["phi"]


def steady_state(P):
    a, b, dk, dh, tau, th, phi, gsh = (P["alpha"], P["beta"], P["dk"], P["dh"],
                                       P["tau"], P["theta"], P["phi"], P["gsh"])
    H = 1.0
    Rtil = (1 / b - 1 + dk) / (1 - tau)
    KY = a / Rtil
    coef = (a / Rtil) ** (a / (1 - a))
    cshare = 1 - dk * KY - gsh
    M = (1 - tau) * (1 - a) / cshare
    l = M / (th + M)
    Y = coef * l
    K = KY * Y
    Ge = gsh * Y
    Cc = cshare * Y
    A = dh / (Ge ** phi)
    I = dk * K
    S = (1 - a) * Y / l
    T = tau * Y
    r = a * Y / K
    return dict(H=H, K=K, l=l, Y=Y, Ge=Ge, C=Cc, A=A, I=I, S=S, T=T, r=r, KY=KY)


ss = steady_state(P)
A = ss["A"]
Kss, Css, lss = ss["K"], ss["C"], ss["l"]
KAPPA = C.PIB_BR_BI / ss["Y"]           # fator de escala modelo -> R$ bilhões
Rtil = (1 / b - 1 + dk) / (1 - tau)
KY = a / Rtil
Bcoef = (a / Rtil) ** (a / (1 - a))     # Y = Bcoef * H * l no estado estacionário
print("Estado estacionário:", {k: round(ss[k], 4) for k in
      ["Y", "K", "H", "l", "C", "I", "Ge", "S", "T", "r"]})


def steady_state_for_Ge(Ge_level):
    """Estado estacionário para um NÍVEL FIXO de gasto educacional Ge (solução
    fechada). Usado como condição terminal na transição, essencial para cortes
    PERMANENTES (o EE muda). Reproduz o EE base quando Ge = Ge_ss."""
    H = A * Ge_level ** phi / dh
    l = ((1 - tau) * (1 - a) * Bcoef * H + th * Ge_level) / (
        Bcoef * H * (th * (1 - dk * KY) + (1 - tau) * (1 - a)))
    Y = Bcoef * H * l
    K = KY * Y
    Cc = Y * (1 - dk * KY) - Ge_level
    I = dk * K
    S = (1 - a) * Y / l
    T = tau * Y
    return dict(H=H, K=K, l=l, Y=Y, Ge=Ge_level, C=Cc, I=I, S=S, T=T)

# ===========================================================================
# SOLUÇÃO POR PREVISÃO PERFEITA NÃO-LINEAR
# ===========================================================================
TT = 220


def Yfun(K, H, l):
    return K ** a * (H * l) ** (1 - a)


def human_capital_path(Ge_path):
    H = np.empty(TT + 2)
    H[0] = ss["H"]
    for t in range(TT + 1):
        H[t + 1] = (1 - dh) * H[t] + A * (Ge_path[min(t, TT)] ** phi)
    return H


def solve_path(Ge_path):
    """Resolve a transição; retorna (dev%, níveis, sucesso, resíduo).
    Condição terminal = estado estacionário do NÍVEL FINAL de gasto (crucial
    para cortes permanentes, em que o EE se desloca)."""
    H_path = human_capital_path(Ge_path)
    tss = steady_state_for_Ge(Ge_path[TT])       # EE terminal (novo)
    Cterm, lterm = tss["C"], tss["l"]

    def residuals(x):
        Cc = x[0:TT]; l = x[TT:2 * TT]; Kn = x[2 * TT:3 * TT]
        K = np.empty(TT + 1); K[0] = Kss; K[1:] = Kn
        Y = Yfun(K[:TT], H_path[:TT], l)
        res = np.empty(3 * TT)
        res[0:TT] = th / (1 - l) - (1 - tau) * (1 - a) * Y / (Cc * l)
        res[TT:2 * TT] = Cc + Kn - (1 - dk) * K[:TT] + Ge_path[:TT] - Y
        l_next = np.append(l[1:], lterm); C_next = np.append(Cc[1:], Cterm)
        Ynext = Yfun(Kn, H_path[1:TT + 1], l_next)
        Rnext = (1 - tau) * a * Ynext / Kn + (1 - dk)
        res[2 * TT:3 * TT] = 1.0 / Cc - b / C_next * Rnext
        return res

    x0 = np.concatenate([np.full(TT, Cterm), np.full(TT, lterm),
                         np.full(TT, tss["K"])])
    sol = root(residuals, x0, method="lm", options={"maxiter": 8000, "xtol": 1e-11})
    Cc = sol.x[0:TT]; l = sol.x[TT:2 * TT]; Kn = sol.x[2 * TT:3 * TT]
    K = np.append(Kss, Kn)
    Y = Yfun(K[:TT], H_path[:TT], l)
    I = Kn - (1 - dk) * K[:TT]
    S = (1 - a) * Y / l
    Tax = tau * Y
    Ge = Ge_path[:TT]
    Hh = H_path[:TT]

    def dev(x, x0):
        return 100 * (x / x0 - 1)

    devs = pd.DataFrame({
        "t": np.arange(TT),
        "Produto": dev(Y, ss["Y"]), "Investimento": dev(I, ss["I"]),
        "Salario": dev(S, ss["S"]), "Consumo": dev(Cc, ss["C"]),
        "Arrecadacao": dev(Tax, ss["T"]), "CapitalHumano": dev(Hh, ss["H"]),
        "GastoEduc": dev(Ge, ss["Ge"]), "Emprego": dev(l, ss["l"]),
    })
    levels = pd.DataFrame({
        "t": np.arange(TT), "Y": Y, "I": I, "S": S, "C": Cc, "Tax": Tax,
        "H": Hh, "Ge": Ge, "l": l, "K": K[:TT],
    })
    return devs, levels, sol.success, float(np.max(np.abs(sol.fun)))


# --- Construção dos caminhos de gasto -------------------------------------
def path_permanent(cut):
    return np.full(TT + 1, ss["Ge"] * (1 - cut))


def path_transitory(cut, rho):
    return ss["Ge"] * np.exp(-cut * rho ** np.arange(TT + 1))


def path_ceiling(total=0.10, ramp=8):
    """Teto de gasto (EC 95): queda gradual do gasto educ. pc até -total% em
    `ramp` anos, permanecendo comprimido depois."""
    frac = np.minimum(np.arange(TT + 1) / ramp, 1.0)
    return ss["Ge"] * (1 - total * frac)


# ===========================================================================
# (1) CENÁRIOS DE CORTE
# ===========================================================================
CUTS = [0.01, 0.05, 0.10, 0.25]
results = {}   # cut -> (devs, levels)
for cut in CUTS:
    d, lv, ok, r = solve_path(path_permanent(cut))
    results[cut] = (d, lv)
    print(f"Corte permanente {int(cut*100):>2d}%: conv={ok}, max|res|={r:.1e}, "
          f"H_LP={d['CapitalHumano'].iloc[-1]:.3f}%, Y_LP={d['Produto'].iloc[-1]:.3f}%")

# transitório (destaque, corte 1%) e teto de gasto
devs_tr, lev_tr, _, _ = solve_path(path_transitory(0.01, P["rhoG"]))
devs_ec, lev_ec, _, _ = solve_path(path_ceiling(0.10, 8))

# cenário-base de destaque = corte permanente de 1%
irf, lev1 = results[0.01]
irf.to_csv(os.path.join(C.OUT, "irf_dsge.csv"), index=False)
devs_tr.to_csv(os.path.join(C.OUT, "irf_dsge_transitorio.csv"), index=False)

VARS = ["Produto", "Investimento", "Salario", "Consumo", "Arrecadacao"]


def cum_loss(series, n):
    return C.cum_disc(series.values, n)


# ===========================================================================
# (4) BEM-ESTAR (equivalente de consumo) E MULTIPLICADORES
# ===========================================================================
def welfare_cev(levels):
    """Perda de bem-estar em equivalente de consumo (% ), utilidade
    U = ln C + theta ln(1-l), com cauda no EE terminal."""
    Cc = levels["C"].values; l = levels["l"].values
    u = np.log(Cc) + th * np.log(1 - l)
    disc = b ** np.arange(TT)
    W = float(np.sum(disc * u)) + b ** TT * u[-1] / (1 - b)
    u_ss = np.log(Css) + th * np.log(1 - lss)
    W_ss = u_ss / (1 - b)
    lam = np.exp((W - W_ss) * (1 - b)) - 1.0    # <0 => perda
    return -lam * 100.0                          # perda em % do consumo permanente


def multipliers(levels, n_list=(1, 5, 10, 20, 30)):
    dY = levels["Y"].values - ss["Y"]
    dGe = levels["Ge"].values - ss["Ge"]
    out = {"impacto": float(dY[0] / dGe[0])}
    for n in n_list:
        w = C.DISC ** np.arange(n)
        out[f"cum{n}"] = float(np.sum(w * dY[:n]) / np.sum(w * dGe[:n]))
    return out


# ===========================================================================
# (2) CUSTO FISCAL EM R$ E CUSTO LÍQUIDO (LAFFER DA EDUCAÇÃO)
# ===========================================================================
def fiscal_accounts(levels, horizons=(10, 20, 30)):
    dGe = (levels["Ge"].values - ss["Ge"]) * KAPPA      # R$ bi (negativo)
    dTax = (levels["Tax"].values - ss["T"]) * KAPPA     # R$ bi (negativo)
    saving_flow = -dGe            # economia orçamentária anual (R$ bi, +)
    revloss_flow = -dTax          # perda de arrecadação anual (R$ bi, +)
    out = {}
    for n in horizons:
        w = C.DISC_SOCIAL ** np.arange(n)
        pv_save = float(np.sum(w * saving_flow[:n]))
        pv_loss = float(np.sum(w * revloss_flow[:n]))
        out[n] = dict(pv_save=pv_save, pv_loss=pv_loss,
                      net=pv_save - pv_loss,
                      eroded=100 * pv_loss / pv_save if pv_save > 0 else np.nan)
    # ano de break-even (perda acumulada de arrecadação supera economia acum.)
    cum_save = np.cumsum(saving_flow)
    cum_loss_ = np.cumsum(revloss_flow)
    be = np.where(cum_loss_ > cum_save)[0]
    out["breakeven"] = int(be[0]) if len(be) else None
    out["saving_flow"] = saving_flow
    out["revloss_flow"] = revloss_flow
    return out


scen_metrics = {}
for cut in CUTS:
    d, lv = results[cut]
    scen_metrics[cut] = dict(
        H_LP=float(d["CapitalHumano"].iloc[-1]),
        Y_LP=float(d["Produto"].iloc[-1]),
        S_LP=float(d["Salario"].iloc[-1]),
        T_LP=float(d["Arrecadacao"].iloc[-1]),
        C_LP=float(d["Consumo"].iloc[-1]),
        loss_Y10=cum_loss(d["Produto"], 10), loss_Y20=cum_loss(d["Produto"], 20),
        loss_S10=cum_loss(d["Salario"], 10), loss_T10=cum_loss(d["Arrecadacao"], 10),
        cev=welfare_cev(lv), mult=multipliers(lv), fisc=fiscal_accounts(lv),
    )
    m = scen_metrics[cut]
    print(f"  {int(cut*100):>2d}%: CEV={m['cev']:.3f}%  mult_cum10={m['mult']['cum10']:.2f}  "
          f"fisc10 net=R$ {C.brl(m['fisc'][10]['net'])} bi  eroded10={m['fisc'][10]['eroded']:.1f}%")

# perdas detalhadas do cenário 1% (para tabela clássica)
losses = {v: {"h10": cum_loss(irf[v], 10), "h20": cum_loss(irf[v], 20)} for v in VARS}

# ===========================================================================
# (3) PROJEÇÃO REGIONAL + DESIGUALDADE
# ===========================================================================
el = pd.read_csv(os.path.join(C.OUT, "elasticidades_uf.csv"))
painel = pd.read_csv(os.path.join(C.OUT, "painel_calibrado.csv"))
# salário-base por UF = média do último ano observado
base_wage = (painel[painel["ano"] == painel["ano"].max()]
             .set_index("uf")["salario_medio"])
el = el.set_index("uf")
el["pop"] = pd.Series(C.POP_UF)
el["base_wage"] = base_wage
el["regiao"] = pd.Series(C.REGIAO)
norm = el["elast_eb"].mean()

# resposta salarial agregada do DSGE (perda acumulada 10a) e razão por UF
agg_sal_loss10 = cum_loss(irf["Salario"], 10)
agg_sal_loss20 = cum_loss(irf["Salario"], 20)
el["perda_salarial_10a"] = agg_sal_loss10 * el["elast_eb"] / norm
el["perda_salarial_20a"] = agg_sal_loss20 * el["elast_eb"] / norm
el = el.sort_values("perda_salarial_10a", ascending=False)
el.reset_index().to_csv(os.path.join(C.OUT, "projecao_regional.csv"), index=False)

# --- Desigualdade: aplica a queda salarial de LONGO PRAZO por UF -------------
def inequality_after(cut):
    S_LP = scen_metrics[cut]["S_LP"]          # % (negativo) resposta agregada LP
    pct = S_LP * (el["elast_eb"] / norm) / 100.0
    w_before = el["base_wage"].values
    w_after = w_before * (1 + pct.values)
    pop = el["pop"].values
    g0 = C.gini(w_before, pop); g1 = C.gini(w_after, pop)
    # razão salário médio Sudeste / Nordeste (ponderado por população)
    def wmean(mask, w):
        return np.average(w[mask], weights=pop[mask])
    se = (el["regiao"] == "SE").values; ne = (el["regiao"] == "NE").values
    ratio0 = wmean(se, w_before) / wmean(ne, w_before)
    ratio1 = wmean(se, w_after) / wmean(ne, w_after)
    return dict(g0=g0, g1=g1, ratio0=ratio0, ratio1=ratio1,
                w_before=w_before, w_after=w_after, pop=pop)


ineq = {cut: inequality_after(cut) for cut in CUTS}
for cut in CUTS:
    q = ineq[cut]
    print(f"  Desigualdade {int(cut*100):>2d}%: Gini {q['g0']:.4f} -> {q['g1']:.4f} "
          f"(+{100*(q['g1']/q['g0']-1):.2f}%) | SE/NE {q['ratio0']:.3f} -> {q['ratio1']:.3f}")

# ===========================================================================
# FIGURAS
# ===========================================================================
# --- F1: IRF permanente vs transitório (cenário 1%) ------------------------
vars6 = [("Produto", C.VAR_COL["Produto"]), ("Investimento", C.VAR_COL["Investimento"]),
         ("Salario", C.VAR_COL["Salario"]), ("Consumo", C.VAR_COL["Consumo"]),
         ("Arrecadacao", C.VAR_COL["Arrecadacao"]), ("CapitalHumano", C.VAR_COL["CapitalHumano"])]
fig, axes = plt.subplots(2, 3, figsize=(11, 6.4))
for ax, (v, c) in zip(axes.ravel(), vars6):
    ax.plot(irf["t"][:35], irf[v][:35], color=c, lw=2.2, label="Corte permanente")
    ax.plot(devs_tr["t"][:35], devs_tr[v][:35], color=c, lw=1.5, ls="--",
            alpha=0.8, label="Corte transitório (AR1)")
    ax.axhline(0, color="black", ls=":", lw=0.7)
    ax.set_title(v if v != "CapitalHumano" else "Capital humano")
    ax.set_xlabel("Anos"); ax.set_ylabel("Desvio % do EE")
axes.ravel()[0].legend(fontsize=7, loc="upper right")
fig.suptitle("Resposta a corte de 1% no gasto público em educação — DSGE fiscal com capital humano",
             fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "irf_corte_educacao.png")); plt.close()

# --- F2: NOVA — múltiplos cortes (Produto, Salário, Arrecadação, Capital humano) ---
fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
panel_vars = [("Produto", "Produto"), ("Salario", "Salário"),
              ("Arrecadacao", "Arrecadação"), ("CapitalHumano", "Capital humano")]
for ax, (v, tit) in zip(axes.ravel(), panel_vars):
    for cut in CUTS:
        d, _ = results[cut]
        ax.plot(d["t"][:40], d[v][:40], color=C.CUT_COL[cut], lw=2,
                label=f"Corte {int(cut*100)}%")
    ax.axhline(0, color="black", ls=":", lw=0.7)
    ax.set_title(tit); ax.set_xlabel("Anos"); ax.set_ylabel("Desvio % do EE")
axes.ravel()[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Trajetórias por magnitude do corte permanente em educação", fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "irf_multi_corte.png")); plt.close()

# --- F3: NOVA — perda de longo prazo por magnitude do corte ----------------
fig, ax = plt.subplots(figsize=(7.6, 4.6))
xlab = [f"{int(c*100)}%" for c in CUTS]
x = np.arange(len(CUTS)); wbar = 0.2
for i, (key, lab, col) in enumerate([
        ("H_LP", "Capital humano", C.VAR_COL["CapitalHumano"]),
        ("Y_LP", "Produto", C.VAR_COL["Produto"]),
        ("S_LP", "Salário", C.VAR_COL["Salario"]),
        ("T_LP", "Arrecadação", C.VAR_COL["Arrecadacao"])]):
    vals = [scen_metrics[c][key] for c in CUTS]
    ax.bar(x + (i - 1.5) * wbar, vals, wbar, label=lab, color=col,
           edgecolor="black", linewidth=0.3)
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(xlab)
ax.set_xlabel("Magnitude do corte permanente")
ax.set_ylabel("Desvio de longo prazo (% do EE)")
ax.set_title("Perda estrutural de longo prazo por magnitude do corte")
ax.legend(fontsize=8, ncol=2)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "perdas_por_corte.png")); plt.close()

# --- F4: NOVA — custo fiscal líquido no tempo (corte 10%) ------------------
fisc10 = scen_metrics[0.10]["fisc"]
HZ = 40
yrs = np.arange(HZ)
cum_save = np.cumsum(fisc10["saving_flow"][:HZ])
cum_loss_r = np.cumsum(fisc10["revloss_flow"][:HZ])
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.plot(yrs, cum_save, color="#2166ac", lw=2.2, label="Economia orçamentária acumulada")
ax.plot(yrs, cum_loss_r, color="#b2182b", lw=2.2, label="Perda de arrecadação acumulada")
ax.fill_between(yrs, cum_save, cum_loss_r, where=(cum_save >= cum_loss_r),
                color="#2166ac", alpha=0.12)
ax.fill_between(yrs, cum_save, cum_loss_r, where=(cum_save < cum_loss_r),
                color="#b2182b", alpha=0.18)
be = scen_metrics[0.10]["fisc"]["breakeven"]
if be is not None:
    ax.axvline(be, color="black", ls="--", lw=1)
    ax.text(be + 0.3, ax.get_ylim()[1] * 0.15, f"Break-even: ano {be}", fontsize=9)
ax.set_xlabel("Anos"); ax.set_ylabel("R\\$ bilhões (valor presente)")
ax.set_title("Custo fiscal líquido de um corte permanente de 10% em educação")
ax.legend(fontsize=9, loc="upper left")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "custo_fiscal.png")); plt.close()

# --- F5: NOVA — economia vs perda de arrecadação (VP 20a) por corte --------
fig, ax = plt.subplots(figsize=(7.8, 4.6))
x = np.arange(len(CUTS)); wbar = 0.28
pv_save = [scen_metrics[c]["fisc"][20]["pv_save"] for c in CUTS]
pv_loss = [scen_metrics[c]["fisc"][20]["pv_loss"] for c in CUTS]
net = [scen_metrics[c]["fisc"][20]["net"] for c in CUTS]
ax.bar(x - wbar, pv_save, wbar, label="Economia orçamentária (VP 20a)",
       color="#2166ac", edgecolor="black", linewidth=0.3)
ax.bar(x, pv_loss, wbar, label="Perda de arrecadação (VP 20a)",
       color="#b2182b", edgecolor="black", linewidth=0.3)
ax.bar(x + wbar, net, wbar, label="Ganho fiscal líquido",
       color="#4d9221", edgecolor="black", linewidth=0.3)
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([f"{int(c*100)}%" for c in CUTS])
ax.set_xlabel("Magnitude do corte permanente")
ax.set_ylabel("R\\$ bilhões (valor presente, 20 anos)")
ax.set_title("A austeridade educacional se autofinancia? Economia vs. perda de base tributária")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "custo_fiscal_barras.png")); plt.close()

# --- F6: projeção regional (perda salarial por UF) -------------------------
elp = el.reset_index().sort_values("perda_salarial_10a")
fig, ax = plt.subplots(figsize=(7.4, 7.6))
ax.barh(np.arange(len(elp)), elp["perda_salarial_10a"],
        color=[C.REG_COL[r] for r in elp["regiao"]], edgecolor="black", linewidth=0.4)
ax.set_yticks(np.arange(len(elp))); ax.set_yticklabels(elp["uf"], fontsize=8)
ax.set_xlabel("Perda salarial acumulada em 10 anos (pontos %-ano, descontada)")
ax.set_title("Prejuízo salarial regional do corte permanente de 1% em educação")
ax.legend(handles=[Patch(facecolor=C.REG_COL[r], label=r) for r in ["N", "NE", "CO", "SE", "S"]],
          fontsize=8, loc="lower right", title="Macrorregião")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "projecao_regional.png")); plt.close()

# --- F7: heatmap por UF ----------------------------------------------------
elh = el.reset_index().sort_values("perda_salarial_10a", ascending=False)
fig, ax = plt.subplots(figsize=(11, 2.5))
im = ax.imshow(elh["perda_salarial_10a"].values.reshape(1, -1), aspect="auto", cmap="OrRd")
ax.set_xticks(np.arange(len(elh))); ax.set_xticklabels(elh["uf"], fontsize=7, rotation=90)
ax.set_yticks([])
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.set_label("Perda salarial 10a", fontsize=8)
ax.set_title("Intensidade do prejuízo salarial por UF (corte permanente de 1% em educação)")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "mapa_calor_uf.png")); plt.close()

# --- F8: NOVA — desigualdade regional (Gini e razão SE/NE) -----------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
x = np.arange(len(CUTS)); wbar = 0.35
g_before = [ineq[c]["g0"] for c in CUTS]; g_after = [ineq[c]["g1"] for c in CUTS]
axes[0].bar(x - wbar / 2, g_before, wbar, label="Antes do corte", color="#5aae61",
            edgecolor="black", linewidth=0.3)
axes[0].bar(x + wbar / 2, g_after, wbar, label="Depois do corte", color="#c2a5cf",
            edgecolor="black", linewidth=0.3)
axes[0].set_xticks(x); axes[0].set_xticklabels([f"{int(c*100)}%" for c in CUTS])
axes[0].set_ylim(min(g_before) * 0.98, max(g_after) * 1.02)
axes[0].set_xlabel("Magnitude do corte"); axes[0].set_ylabel("Gini dos salários entre UF")
axes[0].set_title("Desigualdade salarial regional (Gini)")
axes[0].legend(fontsize=8)
r_before = [ineq[c]["ratio0"] for c in CUTS]; r_after = [ineq[c]["ratio1"] for c in CUTS]
axes[1].bar(x - wbar / 2, r_before, wbar, label="Antes", color="#5aae61",
            edgecolor="black", linewidth=0.3)
axes[1].bar(x + wbar / 2, r_after, wbar, label="Depois", color="#c2a5cf",
            edgecolor="black", linewidth=0.3)
axes[1].set_xticks(x); axes[1].set_xticklabels([f"{int(c*100)}%" for c in CUTS])
axes[1].set_ylim(min(r_before) * 0.98, max(r_after) * 1.02)
axes[1].set_xlabel("Magnitude do corte"); axes[1].set_ylabel("Razão salário SE / NE")
axes[1].set_title("Distância salarial Sudeste-Nordeste")
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "desigualdade_gini.png")); plt.close()

# --- F9: NOVA — curva de Lorenz antes/depois (corte 25%) -------------------
q = ineq[0.25]
p0, L0 = C.lorenz(q["w_before"], q["pop"])
p1, L1 = C.lorenz(q["w_after"], q["pop"])
fig, ax = plt.subplots(figsize=(5.6, 5.4))
ax.plot([0, 1], [0, 1], color="black", lw=1, ls=":", label="Igualdade perfeita")
ax.plot(p0, L0, color="#5aae61", lw=2.2, label=f"Antes (Gini={q['g0']:.3f})")
ax.plot(p1, L1, color="#762a83", lw=2.2, label=f"Depois (Gini={q['g1']:.3f})")
ax.fill_between(p0, p0, L0, color="#5aae61", alpha=0.08)
ax.set_xlabel("Proporção acumulada da população (UF ordenadas por salário)")
ax.set_ylabel("Proporção acumulada da massa salarial")
ax.set_title("Curva de Lorenz dos salários entre UF (corte de 25%)")
ax.legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "lorenz.png")); plt.close()

# --- F10: NOVA — bem-estar e multiplicador ---------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
cev = [scen_metrics[c]["cev"] for c in CUTS]
axes[0].bar([f"{int(c*100)}%" for c in CUTS], cev,
            color=[C.CUT_COL[c] for c in CUTS], edgecolor="black", linewidth=0.3)
axes[0].set_xlabel("Magnitude do corte"); axes[0].set_ylabel("Perda de bem-estar (% cons. equiv.)")
axes[0].set_title("Custo de bem-estar do corte em educação")
for i, c in enumerate(CUTS):
    axes[0].text(i, cev[i], f"{cev[i]:.2f}%", ha="center", va="bottom", fontsize=8)
hs = [1, 5, 10, 20, 30]
mult10 = [scen_metrics[0.10]["mult"][f"cum{n}"] for n in hs]
axes[1].plot(hs, mult10, marker="o", color="#1f78b4", lw=2.2)
axes[1].axhline(1, color="black", ls=":", lw=0.8)
axes[1].set_xlabel("Horizonte (anos)"); axes[1].set_ylabel("Multiplicador acumulado (VP)")
axes[1].set_title("Multiplicador do gasto educacional sobre o produto (corte 10%)")
for xh, ym in zip(hs, mult10):
    axes[1].text(xh, ym, f"{ym:.2f}", ha="center", va="bottom", fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "bemestar_multiplicador.png")); plt.close()

# --- F11: NOVA — trajetória de 30 anos (corte 10% + teto EC-95) ------------
d10, _ = results[0.10]
fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
for ax, (v, tit) in zip(axes.ravel(), panel_vars):
    ax.plot(d10["t"][:30], d10[v][:30], color=C.CUT_COL[0.10], lw=2.4,
            label="Corte permanente 10%")
    ax.plot(devs_ec["t"][:30], devs_ec[v][:30], color="#000000", lw=1.8, ls="--",
            label="Teto de gasto (EC 95, -10% em 8 anos)")
    ax.axhline(0, color="black", ls=":", lw=0.7)
    ax.set_title(tit); ax.set_xlabel("Anos"); ax.set_ylabel("Desvio % do EE")
axes.ravel()[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Trajetórias de longo prazo (30 anos): corte abrupto vs. teto de gasto gradual",
             fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "trajetoria_30anos.png")); plt.close()

# ===========================================================================
# TABELAS LaTeX
# ===========================================================================
# T1: parâmetros
with open(os.path.join(C.TAB, "dsge_params.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{cll}\n\\toprule\n")
    f.write("Parâmetro & Valor & Interpretação \\\\\n\\midrule\n")
    rows = [("$\\alpha$", "0,40", "Participação do capital físico na produção"),
            ("$\\beta$", "0,96", "Fator de desconto intertemporal (anual)"),
            ("$\\delta_k$", "0,10", "Depreciação do capital físico"),
            ("$\\delta_h$", "0,04", "Depreciação do capital humano"),
            ("$\\tau$", "0,33", "Carga tributária efetiva (\\% do PIB)"),
            ("$\\theta$", "1,75", "Peso do lazer na utilidade"),
            ("$\\phi$", "0,30", "Elasticidade do capital humano ao gasto público"),
            ("$g_e$", "0,05", "Gasto público em educação / PIB no EE"),
            ("$\\rho_G$", "0,90", "Persistência do choque fiscal (transitório)"),
            ("$A$", C.num(A, 3), "Produtividade da tecnologia educacional (calibrada)")]
    for p, v, i in rows:
        f.write(f"{p} & {v} & {i} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T2: IRF detalhada (corte permanente 1%)
with open(os.path.join(C.TAB, "dsge_irf.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
    f.write("Variável & Impacto (t=0) & 5 anos & Perda acum. 10a & Perda acum. 20a \\\\\n\\midrule\n")
    for v in VARS:
        f.write(f"{v} & {C.num(irf[v].iloc[0])} & {C.num(irf[v].iloc[5])} & "
                f"{C.num(losses[v]['h10'])} & {C.num(losses[v]['h20'])} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T3: NOVA — cenários de corte (LP + bem-estar + multiplicador)
with open(os.path.join(C.TAB, "cenarios_corte.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
    f.write("Corte & $\\Delta H_{LP}$ & $\\Delta Y_{LP}$ & $\\Delta S_{LP}$ & "
            "$\\Delta T_{LP}$ & Bem-estar & Mult.\\ cum.\\ 10a \\\\\n")
    f.write(" & (\\%) & (\\%) & (\\%) & (\\%) & (\\% cons.) & \\\\\n\\midrule\n")
    for c in CUTS:
        m = scen_metrics[c]
        f.write(f"{int(c*100)}\\% & {C.num(m['H_LP'])} & {C.num(m['Y_LP'])} & "
                f"{C.num(m['S_LP'])} & {C.num(m['T_LP'])} & {C.num(m['cev'])} & "
                f"{C.num(m['mult']['cum10'], 2)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T4: NOVA — custo fiscal em R$ bilhões
with open(os.path.join(C.TAB, "custo_fiscal.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrrr}\n\\toprule\n")
    f.write("Corte & Horizonte & Economia (VP) & Perda arrec.\\ (VP) & "
            "Líquido (VP) & \\% erodido \\\\\n")
    f.write(" & (anos) & R\\$ bi & R\\$ bi & R\\$ bi & \\\\\n\\midrule\n")
    for c in CUTS:
        for n in (10, 20, 30):
            fi = scen_metrics[c]["fisc"][n]
            f.write(f"{int(c*100)}\\% & {n} & {C.brl(fi['pv_save'])} & "
                    f"{C.brl(fi['pv_loss'])} & {C.brl(fi['net'])} & "
                    f"{C.num(fi['eroded'], 1)}\\% \\\\\n")
        f.write("\\midrule\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T5: NOVA — desigualdade regional
with open(os.path.join(C.TAB, "desigualdade.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrrr}\n\\toprule\n")
    f.write("Corte & Gini antes & Gini depois & $\\Delta$ Gini (\\%) & "
            "SE/NE antes & SE/NE depois \\\\\n\\midrule\n")
    for c in CUTS:
        p = ineq[c]
        f.write(f"{int(c*100)}\\% & {C.num(p['g0'], 4)} & {C.num(p['g1'], 4)} & "
                f"{C.num(100*(p['g1']/p['g0']-1), 2)} & {C.num(p['ratio0'], 3)} & "
                f"{C.num(p['ratio1'], 3)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T6: NOVA — multiplicadores por horizonte
with open(os.path.join(C.TAB, "multiplicador.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
    f.write("Corte & Impacto & Cum.\\ 1a & Cum.\\ 5a & Cum.\\ 10a & "
            "Cum.\\ 20a & Cum.\\ 30a \\\\\n\\midrule\n")
    for c in CUTS:
        m = scen_metrics[c]["mult"]
        f.write(f"{int(c*100)}\\% & {C.num(m['impacto'], 2)} & {C.num(m['cum1'], 2)} & "
                f"{C.num(m['cum5'], 2)} & {C.num(m['cum10'], 2)} & "
                f"{C.num(m['cum20'], 2)} & {C.num(m['cum30'], 2)} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# T7: projeção regional (12 UF mais afetadas)
with open(os.path.join(C.TAB, "projecao_regional.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{llrrr}\n\\toprule\n")
    f.write("UF & Região & Elast.\\ EB & Perda sal.\\ 10a & Perda sal.\\ 20a \\\\\n\\midrule\n")
    top = el.reset_index().sort_values("perda_salarial_10a", ascending=False).head(12)
    for _, r in top.iterrows():
        f.write(f"{r['uf']} & {r['regiao']} & {C.num(r['elast_eb'])} & "
                f"{C.num(r['perda_salarial_10a'])} & {C.num(r['perda_salarial_20a'])} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# dumps auxiliares
with open(os.path.join(C.OUT, "ss.json"), "w", encoding="utf-8") as f:
    json.dump({k: float(v) for k, v in ss.items()}, f, indent=2)
with open(os.path.join(C.OUT, "dsge_metrics.json"), "w", encoding="utf-8") as f:
    dump = {str(int(c*100)): {k: (v if not isinstance(v, dict) else
            {kk: (vv if not isinstance(vv, np.ndarray) else None)
             for kk, vv in v.items()})
            for k, v in scen_metrics[c].items() if k != "fisc"} for c in CUTS}
    json.dump(dump, f, indent=2, default=float)

print("\nDSGE concluído. Figuras em paper/figures, tabelas em paper/tables.")
print(f"  Fator de escala R$: 1 unidade-modelo = R$ {C.brl(KAPPA)} bi | "
      f"PIB âncora = R$ {C.brl(C.PIB_BR_BI)} bi")
