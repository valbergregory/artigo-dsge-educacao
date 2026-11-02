# Gasto Público em Educação, Salários Regionais e Dinâmica Fiscal no Brasil

Evidências em painel e um DSGE fiscal com capital humano.

Este repositório contém o código que estima o efeito do gasto público em
educação sobre os salários médios estaduais e sobre a dinâmica fiscal no
Brasil, com atenção à heterogeneidade entre macrorregiões e a cenários de
corte do gasto. O trabalho tem dois blocos:

1. **Bloco empírico** — painel das 27 unidades da federação, 2010–2022:
   efeitos fixos bidirecionais, erros-padrão agrupados e de Driscoll–Kraay,
   teste de Hausman, *bootstrap* de elasticidades por UF, elasticidades por
   macrorregião e encolhimento empírico-bayesiano (James–Stein).
2. **Bloco DSGE** — modelo neoclássico com capital físico e humano,
   tributação distorciva e gasto público em educação como insumo da
   acumulação de capital humano (Uzawa–Lucas / Glomm–Ravikumar com governo),
   resolvido por previsão perfeita não linear (*stacked-time* Newton) com
   estado estacionário terminal endógeno. Cenários de corte de 1, 5, 10 e
   25% do gasto e de teto de gasto; efeito Laffer, bem-estar (equivalente de
   consumo), multiplicadores e desigualdade regional (Gini, Lorenz, Theil).

Perguntas de pesquisa: qual é a elasticidade dos salários estaduais ao gasto
público em educação e como ela varia entre regiões; e quanto um corte nesse
gasto custa, ao longo do tempo, em produto, salários e arrecadação.

## Dados e fontes

Por padrão o código gera uma **base calibrada** (painel UF × ano simulado
com elasticidade positiva e heterogênea por região, semente fixa). Essa
natureza é declarada no código e no artigo. A coleta de séries oficiais está
implementada no porte R (`code/03_empirical_painel.R`, `use_api <- TRUE`) e
no `code/script_original.R`, com as seguintes fontes:

| Variável | Fonte | Acesso |
|---|---|---|
| IPCA (deflator, base 2010) | Banco Central do Brasil — SGS, série 433 | `api.bcb.gov.br/dados/serie/bcdata.sgs.433` |
| PIB estadual a preços correntes | IBGE — SIDRA, tabela 5938 (variável 37) | pacote `sidrar` |
| População residente | IBGE — SIDRA, tabela 6579 (variável 9324) | pacote `sidrar` |
| Rendimento médio domiciliar per capita (PNAD Contínua) | IBGE — SIDRA, tabela 7453 (variável 10809) | pacote `sidrar` |
| Despesa estadual em educação (função 12, liquidada) | Tesouro Nacional — SICONFI, RREO Anexo 02 | `apidatalake.tesouro.gov.br/ords/siconfi/tt/rreo` |

Os códigos de tabela e variável do SIDRA/SICONFI podem mudar entre edições e
devem ser conferidos antes do uso. Valores monetários são deflacionados pelo
IPCA. A população por UF do Censo 2022 (IBGE) e a âncora de PIB nominal usada
para converter os resultados do DSGE em R$ estão em `code/_common.py`.

Não há credenciais necessárias: as APIs usadas são públicas.

## Estrutura do repositório

```
.
├── README.md
├── LICENSE                        # MIT (código)
├── LICENSING.md                   # MIT para código; CC-BY 4.0 para texto, figuras e tabelas
├── requirements.txt               # dependências Python
├── run_all.py                     # executa 01 -> 02 -> 05
├── code/
│   ├── _common.py                 # constantes (UF, regiões, população), âncora em R$, estilo, utilitários
│   ├── 01_empirical_painel.py     # bloco empírico (Python)
│   ├── 02_dsge_capital_humano.py  # bloco DSGE (Python)
│   ├── 05_extensoes.py            # sensibilidade do efeito Laffer a phi; trajetória do Gini e decomposição de Theil
│   ├── 03_empirical_painel.R      # porte R do bloco empírico, com coleta oficial opcional via API
│   ├── 04_dsge_capital_humano.R   # porte R do bloco DSGE
│   └── script_original.R          # versão R original (coleta via rbcb/sidrar/SICONFI)
├── out/                           # CSV/JSON intermediários gerados pelo código
└── .github/workflows/             # notificação da página pessoal após push
```

Conteúdo de `out/`:

| Arquivo | Conteúdo |
|---|---|
| `painel_calibrado.csv` | painel UF × ano (salário médio, gasto em educação per capita, PIB per capita) |
| `region_elast.csv` | elasticidade por macrorregião |
| `elasticidades_uf.csv` | elasticidades por UF (bootstrap e encolhimento) |
| `summary.json` | resumo do bloco empírico |
| `ss.json` | estado estacionário do DSGE |
| `irf_dsge.csv`, `irf_dsge_transitorio.csv` | trajetórias após corte permanente e transitório |
| `projecao_regional.csv` | projeção regional dos efeitos |
| `dsge_metrics.json` | métricas dos cenários do DSGE |

Figuras (`.png`) e tabelas (`.tex`) são gravadas em `paper/figures/` e
`paper/tables/`; essas pastas e o texto do artigo não são versionados
(ver `.gitignore`).

## Como reproduzir

Python:

```bash
pip install -r requirements.txt
python run_all.py
```

R (RStudio): abra `code/03_empirical_painel.R` e `code/04_dsge_capital_humano.R`,
defina o diretório de trabalho como a pasta `code/` e execute. Os pacotes
(`plm`, `lmtest`, `sandwich`, `boot`, `nleqslv`, `ggplot2`, `dplyr`, `tidyr`;
`jsonlite` e `sidrar` com `use_api <- TRUE`) são instalados pelo script se
ausentes.

## Modelo DSGE (resumo)

- Produção: `Y = K^α (H·L)^{1-α}`.
- Capital humano: `H_{t+1} = (1-δ_h) H_t + A (G^e_t)^φ`.
- Governo com orçamento equilibrado e tributação distorciva: `T = τY`.
- Calibração anual em `code/02_dsge_capital_humano.py` (`P = dict(...)`).

## Licença

Código: MIT ([LICENSE](LICENSE)). Texto, figuras e tabelas: CC-BY 4.0
([LICENSING.md](LICENSING.md)).

## Como citar

> Santos, V. G. B. C. B. (2026). *Gasto Público em Educação, Salários Regionais e
> Dinâmica Fiscal no Brasil: Evidências em Painel e um DSGE Fiscal com Capital
> Humano*. Repositório: https://github.com/valbergregory/artigo-dsge-educacao
