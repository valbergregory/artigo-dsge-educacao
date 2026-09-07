#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py — executa o pipeline completo (bloco empírico -> bloco DSGE).
Uso:  python run_all.py
Gera CSVs em out/, figuras em paper/figures/ e tabelas LaTeX em paper/tables/.
Depois, compile o artigo em paper/main.tex (pdflatex + bibtex).
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")

for script in ("01_empirical_painel.py", "02_dsge_capital_humano.py",
               "05_extensoes.py"):
    print(f"\n{'='*70}\n>>> {script}\n{'='*70}")
    r = subprocess.run([sys.executable, script], cwd=CODE)
    if r.returncode != 0:
        sys.exit(f"Falha em {script} (código {r.returncode})")

print("\nPipeline concluído. Para gerar o PDF:")
print("  cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main")
