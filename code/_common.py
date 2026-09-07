#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_common.py — configurações, constantes e utilitários compartilhados pelos
blocos empírico (01) e DSGE (02). Centraliza: lista de UF, macrorregiões,
populações estaduais (Censo 2022, IBGE, em milhares), âncora macro em R$,
estilo dos gráficos e funções auxiliares.

Projeto: "Gasto Público em Educação, Salários Regionais e Dinâmica Fiscal
no Brasil: Evidências em Painel e um DSGE Fiscal com Capital Humano".
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Caminhos (relativos à raiz do repositório)
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "out")
FIG = os.path.join(ROOT, "paper", "figures")
TAB = os.path.join(ROOT, "paper", "tables")
for d in (OUT, FIG, TAB):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# Geografia
# ---------------------------------------------------------------------------
UFS = ["AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", "MT", "MS",
       "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", "RS", "RO", "RR", "SC",
       "SP", "SE", "TO"]

REGIAO = {
    "AC": "N", "AP": "N", "AM": "N", "PA": "N", "RO": "N", "RR": "N", "TO": "N",
    "AL": "NE", "BA": "NE", "CE": "NE", "MA": "NE", "PB": "NE", "PE": "NE",
    "PI": "NE", "RN": "NE", "SE": "NE",
    "DF": "CO", "GO": "CO", "MT": "CO", "MS": "CO",
    "ES": "SE", "MG": "SE", "RJ": "SE", "SP": "SE",
    "PR": "S", "RS": "S", "SC": "S",
}

REGIAO_NOME = {"N": "Norte", "NE": "Nordeste", "CO": "Centro-Oeste",
               "SE": "Sudeste", "S": "Sul"}

# População residente por UF — Censo Demográfico 2022 (IBGE), em milhares
POP_UF = {
    "AC": 830, "AL": 3128, "AP": 733, "AM": 3941, "BA": 14136, "CE": 8794,
    "DF": 2817, "ES": 3833, "GO": 7056, "MA": 6776, "MT": 3658, "MS": 2757,
    "MG": 20539, "PA": 8121, "PB": 3974, "PR": 11444, "PE": 9058, "PI": 3271,
    "RJ": 16055, "RN": 3303, "RS": 10882, "RO": 1581, "RR": 636, "SC": 7610,
    "SP": 44411, "SE": 2210, "TO": 1511,
}

# ---------------------------------------------------------------------------
# Âncora macroeconômica para conversão em R$ (aprox. 2023)
# ---------------------------------------------------------------------------
PIB_BR_BI = 10900.0      # PIB nominal Brasil, R$ bilhões (~R$ 10,9 trilhões)
DISC = 0.96              # fator de desconto anual para perdas acumuladas
DISC_SOCIAL = 0.96       # idem para valor presente fiscal

# ---------------------------------------------------------------------------
# Paleta e estilo
# ---------------------------------------------------------------------------
REG_COL = {"N": "#1b9e77", "NE": "#d95f02", "CO": "#7570b3",
           "SE": "#e7298a", "S": "#66a61e"}

VAR_COL = {
    "Produto": "#1f77b4", "Investimento": "#2ca02c", "Salario": "#e377c2",
    "Consumo": "#bcbd22", "Arrecadacao": "#d62728", "CapitalHumano": "#9467bd",
    "GastoEduc": "#8c564b", "Emprego": "#17becf",
}

CUT_COL = {0.01: "#ffb14e", 0.05: "#fa8775", 0.10: "#ea5f94", 0.25: "#9d02d7"}


def set_style():
    plt.rcParams.update({
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.dpi": 160, "savefig.dpi": 200, "axes.grid": True,
        "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False,
        "font.family": "DejaVu Sans", "legend.frameon": False,
    })


# ---------------------------------------------------------------------------
# Utilitários numéricos e de formatação
# ---------------------------------------------------------------------------
def cum_disc(series, n, disc=DISC):
    """Soma descontada dos n primeiros desvios (perda reportada positiva)."""
    series = np.asarray(series, dtype=float)
    w = disc ** np.arange(n)
    return -float(np.sum(series[:n] * w))


def gini(x, w=None):
    """Índice de Gini (0-1), opcionalmente ponderado por w (população)."""
    x = np.asarray(x, dtype=float)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w, dtype=float)
    order = np.argsort(x)
    x, w = x[order], w[order]
    cw = np.cumsum(w)
    cxw = np.cumsum(x * w)
    # Gini ponderado (fórmula trapezoidal de Lorenz)
    num = np.sum(cxw[:-1] * w[1:] - cxw[1:] * w[:-1])
    g = (np.sum(w) * cxw[-1] - 2 * np.sum((cw - w) * x * w) - cxw[-1] * np.sum(w))
    # implementação estável via curva de Lorenz
    p = cw / cw[-1]
    L = cxw / cxw[-1]
    p = np.insert(p, 0, 0.0)
    L = np.insert(L, 0, 0.0)
    return float(1 - np.sum((p[1:] - p[:-1]) * (L[1:] + L[:-1])))


def lorenz(x, w=None):
    """Pontos (p, L) da curva de Lorenz ponderada, iniciando em (0,0)."""
    x = np.asarray(x, dtype=float)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w, dtype=float)
    order = np.argsort(x)
    x, w = x[order], w[order]
    cw = np.cumsum(w)
    cxw = np.cumsum(x * w)
    p = np.insert(cw / cw[-1], 0, 0.0)
    L = np.insert(cxw / cxw[-1], 0, 0.0)
    return p, L


def brl(v, casas=1):
    """Formata número em R$ bilhões, padrão brasileiro (vírgula decimal)."""
    s = f"{v:,.{casas}f}"
    return s.replace(",", "X").replace(".", ",").replace("X", ".")


def num(v, casas=3):
    """Número com vírgula decimal (padrão brasileiro) para LaTeX."""
    return f"{v:.{casas}f}".replace(".", ",")


def star(p):
    return "***" if p < 0.01 else ("**" if p < 0.05 else ("*" if p < 0.1 else ""))
