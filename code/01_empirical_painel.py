#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLOCO EMPÍRICO — Painel UF x ano (2010-2022)
============================================
Efeitos fixos bidirecionais, erros-padrão robustos (clustered e Driscoll-Kraay),
teste de Hausman, bootstrap de elasticidades por estado, elasticidades por
macrorregião e ENCOLHIMENTO empírico-bayesiano (James-Stein) das elasticidades
estaduais em direção à elasticidade regional.

Os dados são CALIBRADOS/SIMULADOS de forma realista (elasticidade verdadeira
positiva e heterogênea por região). A natureza calibrada é declarada no artigo;
o roteiro de coleta oficial (BCB/SGS, SIDRA/IBGE, SICONFI/FINBRA) está em
`code/script_original.R` e no porte R, ativável com `use_api <- TRUE`.

Saídas: CSVs em out/, figuras em paper/figures/, tabelas LaTeX em paper/tables/.
"""
import os
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import chi2

from linearmodels.panel import PanelOLS, RandomEffects, PooledOLS
import statsmodels.api as sm

import _common as C

warnings.filterwarnings("ignore")
C.set_style()
RNG = np.random.default_rng(2024)

# ===========================================================================
# 1. GERAÇÃO DA BASE CALIBRADA (27 UF x 13 anos = 351 obs)
# ===========================================================================
anos = np.arange(2010, 2023)
nT = len(anos)

# Elasticidade VERDADEIRA por macrorregião: regiões mais pobres (N/NE) dependem
# MAIS do gasto público em educação (menos substitutos privados de formação de
# capital humano) => maior elasticidade salário-educação.
base_elast = dict(zip(["N", "NE", "CO", "SE", "S"],
                      [0.20, 0.24, 0.16, 0.12, 0.15]))
true_elast = {uf: float(np.clip(base_elast[C.REGIAO[uf]] + RNG.normal(0, 0.06),
                                -0.05, 0.45)) for uf in C.UFS}

# Prêmio salarial e de renda por macrorregião (reflete o padrão brasileiro
# real: Sudeste/Sul mais ricos, Norte/Nordeste mais pobres). Entra apenas no
# NÍVEL (intercepto), absorvido pelos efeitos fixos de UF — NÃO afeta a
# estimativa da elasticidade (identificada within-UF).
reg_wage_premium = {"N": -0.15, "NE": -0.22, "CO": 0.08, "SE": 0.22, "S": 0.14}
reg_pib_premium = {"N": -0.12, "NE": -0.20, "CO": 0.16, "SE": 0.20, "S": 0.14}

# Constantes de centragem: a elasticidade e o controle de PIB afetam apenas os
# DESVIOS em torno da média, de modo que o nível salarial regional fique a cargo
# dos prêmios (evita que a maior elasticidade do NE eleve espuriamente seu nível).
LG0 = np.log(700.0)          # gasto educ. pc médio de referência
LP0 = np.log(28000.0)        # PIB pc médio de referência

pib_level = {uf: np.exp(reg_pib_premium[C.REGIAO[uf]] + RNG.normal(0, 0.10))
             for uf in C.UFS}
gasto_level = {uf: RNG.uniform(0.80, 1.30) for uf in C.UFS}
mu_uf = {uf: reg_wage_premium[C.REGIAO[uf]] + RNG.normal(0, 0.06) for uf in C.UFS}
lam_t = {t: 0.015 * (t - 2010) + RNG.normal(0, 0.01) for t in anos}

rows = []
for uf in C.UFS:
    g0 = 600 * gasto_level[uf]                     # gasto educ. pc inicial (R$/hab)
    for k, t in enumerate(anos):
        gasto = g0 * (1.03) ** k * np.exp(RNG.normal(0, 0.16))
        pib = 28000 * pib_level[uf] * (1.018) ** k * np.exp(RNG.normal(0, 0.05))
        lsal = (7.7 + mu_uf[uf] + lam_t[t]
                + true_elast[uf] * (np.log(gasto) - LG0)
                + 0.35 * (np.log(pib) - LP0)
                + RNG.normal(0, 0.05))
        rows.append((uf, t, C.REGIAO[uf], np.exp(lsal), gasto, pib))

painel = pd.DataFrame(rows, columns=["uf", "ano", "regiao", "salario_medio",
                                     "gasto_edu_pc", "pib_pc"])
painel.to_csv(os.path.join(C.OUT, "painel_calibrado.csv"), index=False)
print(f"Painel: {painel.shape} | obs = {len(painel)}")

# ===========================================================================
# 2. EFEITOS FIXOS BIDIRECIONAIS + INFERÊNCIA ROBUSTA + HAUSMAN
# ===========================================================================
df = painel.copy()
df["lsal"] = np.log(df["salario_medio"])
df["lgasto"] = np.log(df["gasto_edu_pc"])
df["lpib"] = np.log(df["pib_pc"])
df = df.set_index(["uf", "ano"])
exog = sm.add_constant(df[["lgasto", "lpib"]])

mod_fe = PanelOLS(df["lsal"], exog, entity_effects=True, time_effects=True,
                  drop_absorbed=True)
res_fe_cl = mod_fe.fit(cov_type="clustered", cluster_entity=True)
res_fe_dk = mod_fe.fit(cov_type="kernel", kernel="bartlett")
res_re = RandomEffects(df["lsal"], exog).fit()
res_pool = PooledOLS(df["lsal"], exog).fit(cov_type="clustered",
                                           cluster_entity=True)

# tendência linear no lugar de efeitos de ano
trend = pd.Series(df.reset_index()["ano"].values - 2010, index=df.index,
                  name="trend")
res_fe_poly = PanelOLS(df["lsal"], sm.add_constant(
    pd.concat([df[["lgasto", "lpib"]], trend], axis=1)),
    entity_effects=True).fit(cov_type="clustered", cluster_entity=True)

# Hausman FE vs RE
b_fe = res_fe_cl.params[["lgasto", "lpib"]].values
b_re = res_re.params[["lgasto", "lpib"]].values
dV = (res_fe_cl.cov.loc[["lgasto", "lpib"], ["lgasto", "lpib"]].values
      - res_re.cov.loc[["lgasto", "lpib"], ["lgasto", "lpib"]].values)
diff = b_fe - b_re
haus = float(diff @ np.linalg.pinv(dV) @ diff)
haus_p = float(1 - chi2.cdf(haus, df=2))

elast_fe = float(res_fe_cl.params["lgasto"])
print(f"Elast. FE (clustered): {elast_fe:.4f} "
      f"(se {res_fe_cl.std_errors['lgasto']:.4f}, "
      f"p {res_fe_cl.pvalues['lgasto']:.3g})")
print(f"Driscoll-Kraay se: {res_fe_dk.std_errors['lgasto']:.4f}")
print(f"Hausman: {haus:.2f} (p {haus_p:.3g}) | "
      f"R2within {res_fe_cl.rsquared_within:.3f}")

# ===========================================================================
# 3. ELASTICIDADES POR MACRORREGIÃO (FE por região)
# ===========================================================================
region_el = {}
for R in ["N", "NE", "CO", "SE", "S"]:
    sub = painel[painel["regiao"] == R].copy()
    sub["lsal"] = np.log(sub["salario_medio"])
    sub["lgasto"] = np.log(sub["gasto_edu_pc"])
    sub["lpib"] = np.log(sub["pib_pc"])
    sub["trend"] = sub["ano"] - 2010
    sub = sub.set_index(["uf", "ano"])
    r_ = PanelOLS(sub["lsal"], sm.add_constant(sub[["lgasto", "lpib", "trend"]]),
                  entity_effects=True, drop_absorbed=True).fit(
        cov_type="clustered", cluster_entity=True)
    region_el[R] = {"beta": float(r_.params["lgasto"]),
                    "se": float(r_.std_errors["lgasto"]),
                    "p": float(r_.pvalues["lgasto"]), "n": int(r_.nobs)}
pd.DataFrame([{"regiao": R, **v} for R, v in region_el.items()]).to_csv(
    os.path.join(C.OUT, "region_elast.csv"), index=False)
print("Elast. por macrorregião:",
      {R: round(v["beta"], 3) for R, v in region_el.items()})

# ===========================================================================
# 4. BOOTSTRAP DE ELASTICIDADES POR UF + ENCOLHIMENTO EMPÍRICO-BAYESIANO
# ===========================================================================
def elast_uf(sub):
    X = sm.add_constant(np.column_stack([np.log(sub["gasto_edu_pc"]),
                                         np.log(sub["pib_pc"]),
                                         sub["ano"].values - 2010]))
    y = np.log(sub["salario_medio"].values)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta[1]

n_boot = 2000
recs = []
for uf in C.UFS:
    sub = painel[painel["uf"] == uf].reset_index(drop=True)
    point = elast_uf(sub)
    idx = np.arange(len(sub))
    bs = np.array([elast_uf(sub.iloc[RNG.choice(idx, len(sub), replace=True)])
                   for _ in range(n_boot)])
    recs.append({"uf": uf, "regiao": C.REGIAO[uf], "est": float(np.mean(bs)),
                 "point": float(point), "lower": float(np.percentile(bs, 2.5)),
                 "upper": float(np.percentile(bs, 97.5)),
                 "true": true_elast[uf], "boot": bs})
boot_df = pd.DataFrame(recs)

# --- Encolhimento (James-Stein empírico-bayesiano) para a região do estado ----
reg_beta = {R: region_el[R]["beta"] for R in region_el}
boot_df["beta_reg"] = boot_df["regiao"].map(reg_beta)
boot_df["sd"] = (boot_df["upper"] - boot_df["lower"]) / (2 * 1.959964)
boot_df["var"] = boot_df["sd"] ** 2
within = boot_df.groupby("regiao").apply(
    lambda g: max(g["est"].var(ddof=1) - g["var"].mean(), 1e-4)).to_dict()
boot_df["tau2"] = boot_df["regiao"].map(within)
boot_df["w"] = boot_df["tau2"] / (boot_df["tau2"] + boot_df["var"])
boot_df["elast_eb"] = (boot_df["w"] * boot_df["est"]
                       + (1 - boot_df["w"]) * boot_df["beta_reg"]).clip(lower=0.0)

boot_save = boot_df[["uf", "regiao", "est", "lower", "upper", "true",
                     "beta_reg", "w", "elast_eb"]].copy()
boot_save.to_csv(os.path.join(C.OUT, "elasticidades_uf.csv"), index=False)
sig = int(((boot_df["lower"] > 0) | (boot_df["upper"] < 0)).sum())
print(f"UF com IC95% que não cruza zero: {sig}/27 | "
      f"elast média entre-UF: {boot_df['est'].mean():.4f}")

summary = {
    "elast_fe": elast_fe, "se_cl": float(res_fe_cl.std_errors["lgasto"]),
    "se_dk": float(res_fe_dk.std_errors["lgasto"]),
    "p_fe": float(res_fe_cl.pvalues["lgasto"]),
    "lpib_fe": float(res_fe_cl.params["lpib"]),
    "r2within": float(res_fe_cl.rsquared_within),
    "nobs": int(res_fe_cl.nobs), "haus": haus, "haus_p": haus_p,
    "elast_pool": float(res_pool.params["lgasto"]),
    "sig_uf": sig, "elast_eb_mean": float(boot_df["elast_eb"].mean()),
}
with open(os.path.join(C.OUT, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

# ===========================================================================
# 5. FIGURAS
# ===========================================================================
# (a) Elasticidades por UF com IC bootstrap
fig, ax = plt.subplots(figsize=(7.2, 7.6))
bd = boot_df.sort_values("est").reset_index(drop=True)
ypos = np.arange(len(bd))
ax.barh(ypos, bd["est"], color=[C.REG_COL[r] for r in bd["regiao"]],
        alpha=0.85, edgecolor="black", linewidth=0.4)
ax.errorbar(bd["est"], ypos, xerr=[bd["est"] - bd["lower"], bd["upper"] - bd["est"]],
            fmt="none", ecolor="black", elinewidth=0.7, capsize=2)
ax.axvline(0, color="black", lw=0.8)
ax.axvline(elast_fe, color="red", ls="--", lw=1.2)
ax.set_yticks(ypos); ax.set_yticklabels(bd["uf"], fontsize=8)
ax.set_xlabel("Elasticidade salário-educação")
ax.set_title("Elasticidade salário-educação por UF (IC bootstrap 95%)")
leg = [Patch(facecolor=C.REG_COL[r], label=r) for r in ["N", "NE", "CO", "SE", "S"]]
leg.append(plt.Line2D([0], [0], color="red", ls="--", label=f"FE={elast_fe:.3f}"))
ax.legend(handles=leg, fontsize=7, loc="lower right", ncol=2)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "elasticidade_uf_ic.png")); plt.close()

# (b) Histograma bootstrap agregado
allboot = np.concatenate([r["boot"] for r in recs])
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.hist(allboot, bins=60, color="#4c72b0", edgecolor="white", alpha=0.9)
ax.axvline(elast_fe, color="red", ls="--", lw=1.4, label=f"Elast. média FE = {elast_fe:.3f}")
ax.axvline(0, color="black", lw=0.8)
ax.set_xlim(-0.6, 0.8); ax.set_xlabel("Elasticidade"); ax.set_ylabel("Frequência")
ax.set_title("Distribuição bootstrap das elasticidades salário-educação")
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "hist_elasticidades.png")); plt.close()

# (c) Dispersão por macrorregião
fig, axes = plt.subplots(2, 3, figsize=(11, 6.6))
axes = axes.ravel()
for i, R in enumerate(["N", "NE", "CO", "SE", "S"]):
    sub = painel[painel["regiao"] == R]
    ax = axes[i]
    ax.scatter(np.log(sub["gasto_edu_pc"]), np.log(sub["salario_medio"]),
               s=12, alpha=0.45, color=C.REG_COL[R])
    b = np.polyfit(np.log(sub["gasto_edu_pc"]), np.log(sub["salario_medio"]), 1)
    xs = np.linspace(np.log(sub["gasto_edu_pc"]).min(),
                     np.log(sub["gasto_edu_pc"]).max(), 50)
    ax.plot(xs, np.polyval(b, xs), color="black", lw=1.6)
    ax.set_title(f"Região {R} (incl.={b[0]:.2f})", fontsize=10)
    ax.set_xlabel("log(Gasto Educ. pc)"); ax.set_ylabel("log(Salário)")
axes[5].axis("off")
fig.suptitle("Relação salário-educação por macrorregião", fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "dispersao_uf.png")); plt.close()

# (d) NOVA — evolução temporal do painel (médias regionais)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for R in ["N", "NE", "CO", "SE", "S"]:
    sub = painel[painel["regiao"] == R].groupby("ano").mean(numeric_only=True)
    axes[0].plot(sub.index, sub["salario_medio"], marker="o", ms=3,
                 color=C.REG_COL[R], label=C.REGIAO_NOME[R])
    axes[1].plot(sub.index, sub["gasto_edu_pc"], marker="o", ms=3,
                 color=C.REG_COL[R], label=C.REGIAO_NOME[R])
axes[0].set_title("Salário médio real por macrorregião")
axes[0].set_xlabel("Ano"); axes[0].set_ylabel("R\\$ (deflacionado)")
axes[1].set_title("Gasto educ. per capita por macrorregião")
axes[1].set_xlabel("Ano"); axes[1].set_ylabel("R\\$/hab (deflacionado)")
axes[1].legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "evolucao_painel.png")); plt.close()

# (e) NOVA — encolhimento empírico-bayesiano: bruto vs EB
fig, ax = plt.subplots(figsize=(7.4, 7.6))
be = boot_df.sort_values("elast_eb").reset_index(drop=True)
yy = np.arange(len(be))
ax.hlines(yy, be["est"], be["elast_eb"], color="gray", lw=0.8, alpha=0.6)
ax.scatter(be["est"], yy, s=26, color="#bbbbbb", edgecolor="black",
           linewidth=0.3, label="Bruta (bootstrap)", zorder=3)
ax.scatter(be["elast_eb"], yy, s=30, color=[C.REG_COL[r] for r in be["regiao"]],
           edgecolor="black", linewidth=0.3, label="Encolhida (EB)", zorder=4)
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(yy); ax.set_yticklabels(be["uf"], fontsize=8)
ax.set_xlabel("Elasticidade salário-educação")
ax.set_title("Encolhimento empírico-bayesiano das elasticidades estaduais")
ax.legend(fontsize=8, loc="lower right")
plt.tight_layout(); plt.savefig(os.path.join(C.FIG, "encolhimento_eb.png")); plt.close()

# ===========================================================================
# 6. TABELAS LaTeX
# ===========================================================================
# Descritivas
desc = painel[["salario_medio", "gasto_edu_pc", "pib_pc"]].describe().T[
    ["mean", "std", "min", "max"]]
desc.index = ["Salário médio (R\\$)", "Gasto educ. pc (R\\$)", "PIB pc (R\\$)"]
with open(os.path.join(C.TAB, "descritivas.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
    f.write("Variável & Média & Desvio & Mínimo & Máximo \\\\\n\\midrule\n")
    for idx, row in desc.iterrows():
        f.write(f"{idx} & {C.brl(row['mean'])} & {C.brl(row['std'])} & "
                f"{C.brl(row['min'])} & {C.brl(row['max'])} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# Modelo FE (4 colunas: pooled, FE clustered, FE Driscoll-Kraay, FE+tendência)
with open(os.path.join(C.TAB, "modelo_fe.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lcccc}\n\\toprule\n")
    f.write(" & (1) Pooled & (2) FE 2-vias & (3) FE 2-vias & (4) FE entidade \\\\\n")
    f.write(" & OLS & Clustered & Driscoll-Kraay & + tendência \\\\\n\\midrule\n")
    def line(nm, key):
        vals, ses = [], []
        for r in (res_pool, res_fe_cl, res_fe_dk, res_fe_poly):
            vals.append(f"{r.params[key]:.3f}{C.star(r.pvalues[key])}")
            ses.append(f"({r.std_errors[key]:.3f})")
        f.write(f"{nm} & " + " & ".join(vals) + " \\\\\n")
        f.write(" & " + " & ".join(ses) + " \\\\\n")
    line("$\\log$(Gasto educ. pc)", "lgasto")
    line("$\\log$(PIB pc)", "lpib")
    f.write("\\midrule\n")
    f.write("Efeitos fixos UF & Não & Sim & Sim & Sim \\\\\n")
    f.write("Efeitos fixos ano & Não & Sim & Sim & Não \\\\\n")
    f.write(f"Observações & {int(res_pool.nobs)} & {int(res_fe_cl.nobs)} & "
            f"{int(res_fe_dk.nobs)} & {int(res_fe_poly.nobs)} \\\\\n")
    f.write(f"$R^2$ (within) & -- & {res_fe_cl.rsquared_within:.3f} & "
            f"{res_fe_dk.rsquared_within:.3f} & {res_fe_poly.rsquared_within:.3f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# Elasticidades regionais
with open(os.path.join(C.TAB, "region_elast.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lcccc}\n\\toprule\n")
    f.write("Macrorregião & Elasticidade & Erro-padrão & $p$-valor & $N$ \\\\\n\\midrule\n")
    for R in ["N", "NE", "CO", "SE", "S"]:
        v = region_el[R]
        f.write(f"{C.REGIAO_NOME[R]} & {v['beta']:.3f}{C.star(v['p'])} & "
                f"{v['se']:.3f} & {v['p']:.3f} & {v['n']} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# Elasticidades por UF (com EB) — tabela completa
with open(os.path.join(C.TAB, "elasticidades_uf.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{llrrrcr}\n\\toprule\n")
    f.write("UF & Reg. & Elast. & IC inf. & IC sup. & Sig.\\,95\\% & Elast. EB \\\\\n\\midrule\n")
    for _, r in boot_df.sort_values("est", ascending=False).iterrows():
        s = "Sim" if (r["lower"] > 0 or r["upper"] < 0) else "--"
        f.write(f"{r['uf']} & {r['regiao']} & {r['est']:.3f} & {r['lower']:.3f} & "
                f"{r['upper']:.3f} & {s} & {r['elast_eb']:.3f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# NOVA — matriz de correlação (nível log)
corr = df[["lsal", "lgasto", "lpib"]].corr()
labels = {"lsal": "$\\log$ Salário", "lgasto": "$\\log$ Gasto educ.",
          "lpib": "$\\log$ PIB pc"}
with open(os.path.join(C.TAB, "correlacao.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lccc}\n\\toprule\n")
    f.write(" & " + " & ".join(labels[c] for c in corr.columns) + " \\\\\n\\midrule\n")
    for i in corr.index:
        f.write(labels[i] + " & " + " & ".join(f"{corr.loc[i, j]:.3f}"
                for j in corr.columns) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

print("\nBloco empírico concluído.")
print(f"  Elast. média entre-UF: {boot_df['est'].mean():.4f} | "
      f"EB média: {boot_df['elast_eb'].mean():.4f}")
