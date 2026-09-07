# Gasto Público em Educação, Salários Regionais e Dinâmica Fiscal no Brasil

**Evidências em painel e um DSGE fiscal com capital humano.**

Este repositório reúne o código, os dados calibrados, as figuras, as tabelas e o
artigo (LaTeX/Overleaf) de um estudo que avalia o impacto do gasto público em
educação sobre os salários médios estaduais e sobre a dinâmica fiscal no Brasil,
com ênfase na **dimensão regional** e nos efeitos de **cortes/desinvestimentos**
(por exemplo, a limitação de repasses a universidades e institutos federais).

O trabalho combina:

1. **Bloco empírico** — painel das 27 unidades federativas (2010–2022), efeitos
   fixos bidirecionais, erros-padrão de Driscoll–Kraay, teste de Hausman,
   *bootstrap* de elasticidades por UF, elasticidades por macrorregião e
   **encolhimento empírico-bayesiano** (James–Stein).
2. **Bloco DSGE** — modelo neoclássico com capital físico e humano, tributação
   distorciva e gasto público em educação como insumo produtivo da acumulação de
   capital humano (Uzawa–Lucas / Glomm–Ravikumar com governo), resolvido por
   **previsão perfeita não linear** (*stacked-time* Newton) com **estado
   estacionário terminal endógeno**.

## Principais resultados

- Elasticidade salário-educação média de **0,17** (significativa a 1%), com forte
  heterogeneidade regional (≈0,21 no Norte/Nordeste; ≈0 no Sudeste).
- Um corte permanente de **10%** reduz o capital humano de longo prazo em **3,1%**
  e o produto em **3,4%**.
- **Efeito Laffer dinâmico da educação**: o corte corrói a base tributária; em
  valor presente de 30 anos, até **78%** da economia orçamentária é revertida pela
  perda de arrecadação (*break-even* fiscal por volta do ano 33).
- Multiplicador acumulado do gasto educacional sobre o produto **> 1** (≈1,2).
- Perda de bem-estar de até **2,4%** do consumo permanente (corte de 25%).
- A austeridade educacional **amplia a desigualdade regional**: o Gini salarial
  entre estados sobe até **+17,7%** e a razão Sudeste/Nordeste aumenta.

## Artigos

O material está organizado em **três formatos do mesmo trabalho**, para diferentes usos:

| Pasta | Conteúdo | Uso |
|---|---|---|
| `paper/` | **Artigo integrado** (26 pág.): painel + DSGE + custo fiscal + desigualdade. | Versão completa, referência. |
| `anpec/` | Versão **ANPEC** condensada (20 pág., Times 12, A4), nas variantes `main_full` (identificada) e `main_blind` (avaliação às cegas). | Submissão ao ANPEC. |
| `paperA/` | **Paper A — "A austeridade educacional se autofinancia?"**: DSGE fiscal, cenários de corte, efeito Laffer, multiplicador, bem-estar. | Periódico (Setor Público/Macro). |
| `paperB/` | **Paper B — "Cortes em educação e desigualdade regional"**: painel, heterogeneidade, projeção regional, Gini/Lorenz/Theil. | Periódico (Economia Regional). |

Cada paper traz também uma versão **`main_guiado.tex`/`.pdf`** com caixas de orientação (o que escrever/expandir em cada seção); remova-as antes de submeter (`\renewcommand{\guia}[1]{}`).

## Estrutura do repositório

```
.
├── README.md
├── LICENSE                   # MIT (código) + CC-BY 4.0 (texto)
├── requirements.txt          # dependências Python
├── run_all.py                # executa o pipeline Python completo (01 → 02 → 05)
├── code/
│   ├── _common.py            # constantes, populações UF, âncora R$, utilitários
│   ├── 01_empirical_painel.py    # bloco empírico (Python, validado)
│   ├── 02_dsge_capital_humano.py # bloco DSGE (Python, validado)
│   ├── 05_extensoes.py           # extensões: sensibilidade Laffer (φ); Theil e trajetória do Gini
│   ├── 03_empirical_painel.R     # porte R do bloco empírico (RStudio)
│   ├── 04_dsge_capital_humano.R  # porte R do bloco DSGE (RStudio)
│   └── script_original.R         # versão R original com coleta via API (legado)
├── paper/  anpec/  paperA/  paperB/   # ver seção "Artigos" acima
│   ├── main.tex / referencias.bib
│   ├── figures/              # figuras (.png) geradas pelo código
│   └── tables/               # tabelas (.tex) geradas pelo código
└── out/                      # CSVs/JSON intermediários (regerados)
```

## Como reproduzir

### Python (pipeline validado)

```bash
pip install -r requirements.txt
python run_all.py                       # gera figuras e tabelas
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

### R (RStudio)

Abra `code/03_empirical_painel.R` e `code/04_dsge_capital_humano.R` no RStudio,
defina o diretório de trabalho para a pasta `code/` (*Session → Set Working
Directory → To Source File Location*) e execute. Os pacotes necessários
(`plm`, `lmtest`, `sandwich`, `boot`, `nleqslv`, `ggplot2`, `dplyr`, `tidyr`)
são instalados automaticamente.

> **Observação:** os portes em R reproduzem a lógica do pipeline Python (a
> referência numérica validada). Rode-os localmente no RStudio para regenerar as
> versões `*_R.png` / `*_R.tex`.

## Dados oficiais via API (Governo Federal)

Por padrão o código usa uma **base calibrada** realista (declarada abertamente no
artigo), para reprodutibilidade imediata. Para coletar séries **oficiais**,
defina `use_api <- TRUE` em `code/03_empirical_painel.R` (ou consulte
`script_original.R`). As fontes previstas:

| Variável | Fonte | Acesso |
|---|---|---|
| IPCA (deflator) | Banco Central — **SGS** série 433 | `api.bcb.gov.br/dados/serie/bcdata.sgs.433` |
| PIB estadual / rendimento | **SIDRA/IBGE** | pacote `sidrar` (tabelas 5938, 6387) |
| Finanças estaduais (receita/despesa) | **SICONFI/FINBRA** (Tesouro) | API SICONFI |
| Gasto federal em educação | **SIOP / Tesouro** | execução orçamentária (função 12) |

Todos os valores monetários são **deflacionados pelo IPCA** para preços constantes.

## Modelo DSGE (resumo)

- Produção: `Y = K^α (H·L)^{1-α}`.
- Capital humano: `H_{t+1} = (1-δ_h)H_t + A·(G^e_t)^φ` (retornos decrescentes ao
  gasto público; elasticidade de longo prazo `φ`).
- Governo com orçamento equilibrado e tributação distorciva; arrecadação
  `T = τY`.
- Solução por previsão perfeita não linear com EE terminal endógeno (solução
  fechada para nível fixo de gasto — ver apêndice do artigo).

## Como citar

> Santos, V. G. B. C. B. (2026). *Gasto Público em Educação, Salários Regionais e
> Dinâmica Fiscal no Brasil: Evidências em Painel e um DSGE Fiscal com Capital
> Humano*. Repositório: https://github.com/valbergregory/artigo-dsge-educacao

## Licença e uso de IA

Código sob **licença MIT**; texto dos artigos, figuras e tabelas sob **CC-BY 4.0**
(ver `LICENSE`). O artigo foi organizado com auxílio de IA generativa a partir de
código, modelos e instruções do autor; a concepção, a validação e a
responsabilidade científica são integralmente do autor.
