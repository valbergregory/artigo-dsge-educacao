# =============================================================================
# BLOCO EMPÍRICO (porte R do 01_empirical_painel.py) — RStudio
# Painel UF x ano (2010-2022): efeitos fixos bidirecionais, Driscoll-Kraay,
# Hausman, bootstrap de elasticidades por UF, elasticidades por macrorregião e
# encolhimento empírico-bayesiano (James-Stein).
#
# use_api <- TRUE  ativa a coleta de séries OFICIAIS via API do Governo Federal
#   - IPCA (deflator): Banco Central / SGS série 433
#   - PIB e rendimento estadual: SIDRA/IBGE (pacote sidrar)
#   - Finanças estaduais: SICONFI/Tesouro (FINBRA)  |  Gasto federal: SIOP
# Por padrão (use_api <- FALSE) gera-se uma base CALIBRADA realista, idêntica
# em espírito à do porte Python, para reprodutibilidade imediata.
# =============================================================================

## ---- 0. Pacotes -----------------------------------------------------------
pkgs <- c("plm", "lmtest", "sandwich", "boot", "ggplot2", "dplyr", "tidyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(2024)
use_api <- FALSE                      # <- mude para TRUE para coletar dados oficiais

# Raiz do repositório: assume que este script está em <root>/code/.
# Em RStudio, defina o diretório de trabalho para a pasta 'code' (Session ->
# Set Working Directory -> To Source File Location) ou ajuste ROOT manualmente.
get_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  base <- if (length(f)) dirname(normalizePath(f)) else getwd()
  normalizePath(file.path(base, ".."), mustWork = FALSE)
}
ROOT <- get_root()
OUT <- file.path(ROOT, "out"); FIG <- file.path(ROOT, "paper", "figures")
TAB <- file.path(ROOT, "paper", "tables")
for (d in c(OUT, FIG, TAB)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ufs <- c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT","MS","MG","PA",
         "PB","PR","PE","PI","RJ","RN","RS","RO","RR","SC","SP","SE","TO")
regiao <- c(AC="N",AP="N",AM="N",PA="N",RO="N",RR="N",TO="N",
            AL="NE",BA="NE",CE="NE",MA="NE",PB="NE",PE="NE",PI="NE",RN="NE",SE="NE",
            DF="CO",GO="CO",MT="CO",MS="CO",ES="SE",MG="SE",RJ="SE",SP="SE",
            PR="S",RS="S",SC="S")
reg_col <- c(N="#1b9e77", NE="#d95f02", CO="#7570b3", SE="#e7298a", S="#66a61e")
anos <- 2010:2022

## ---- 1. COLETA DE DADOS OFICIAIS VIA API (use_api <- TRUE) -----------------
# Códigos de tabela SIDRA/SICONFI que PODEM variar entre edições das pesquisas.
# Confira em https://sidra.ibge.gov.br e https://apidatalake.tesouro.gov.br .
TAB_PIB        <- 5938   # Contas Regionais — PIB a preços correntes (var 37)
VAR_PIB        <- 37
TAB_POP        <- 6579   # Estimativas de população residente (var 9324)
VAR_POP        <- 9324
TAB_RENDIMENTO <- 7453   # PNADC anual — rendimento médio mensal real domiciliar
VAR_RENDIMENTO <- 10809  #   per capita (var 10809). *** CONFIRA no SIDRA ***
ANOS_PIB       <- 2010:2021   # Contas Regionais saem com defasagem (~2 anos)

# Código IBGE (2 dígitos) de cada UF — usado pelo SICONFI (id_ente).
uf_cod <- c(RO=11,AC=12,AM=13,RR=14,PA=15,AP=16,TO=17,MA=21,PI=22,CE=23,RN=24,
            PB=25,PE=26,AL=27,SE=28,BA=29,MG=31,ES=32,RJ=33,SP=35,PR=41,SC=42,
            RS=43,MS=50,MT=51,GO=52,DF=53)
cod_uf <- setNames(names(uf_cod), as.character(unname(uf_cod)))

# Deflator IPCA pelo SGS do Banco Central (série 433 = IPCA % mensal).
get_ipca_index <- function() {
  url <- paste0("https://api.bcb.gov.br/dados/serie/bcdata.sgs.433/dados",
                "?formato=json&dataInicial=01/01/2009&dataFinal=31/12/2023")
  js <- jsonlite::fromJSON(url)
  js$valor <- as.numeric(gsub(",", ".", js$valor))
  js$data  <- as.Date(js$data, "%d/%m/%Y")
  js$ano   <- as.integer(format(js$data, "%Y"))
  fac <- cumprod(1 + js$valor / 100)          # índice mensal acumulado
  idx <- tapply(fac, js$ano, mean)            # média anual
  idx / idx["2010"]                           # base 2010 = 1
}

# SIDRA/IBGE via pacote sidrar. Retorna data.frame (uf, ano, <col>).
get_sidra_uf <- function(tabela, variavel, periodo, col) {
  d <- sidrar::get_sidra(x = tabela, variable = variavel,
                         period = as.character(periodo), geo = "State")
  cod <- as.character(d[["Unidade da Federação (Código)"]])
  out <- data.frame(uf = cod_uf[cod], ano = as.integer(d[["Ano"]]),
                    valor = suppressWarnings(as.numeric(d[["Valor"]])),
                    stringsAsFactors = FALSE)
  names(out)[3] <- col
  out[!is.na(out$uf), ]
}

# Gasto estadual em EDUCAÇÃO (função 12) via SICONFI/Tesouro — RREO Anexo 02,
# despesas liquidadas. 27 UF x anos => muitas chamadas; usa Sys.sleep leve.
# *** A estrutura de conta/coluna do RREO pode exigir ajuste fino (confira). ***
get_gasto_educ_uf <- function(periodo) {
  base <- "https://apidatalake.tesouro.gov.br/ords/siconfi/tt/rreo"
  recs <- list(); k <- 1
  for (u in ufs) for (a in periodo) {
    url <- sprintf("%s?an_exercicio=%d&nr_periodo=6&co_tipo_demonstrativo=RREO&no_anexo=RREO-Anexo%%2002&co_esfera=E&id_ente=%d",
                   base, a, uf_cod[u])
    val <- tryCatch({
      js <- jsonlite::fromJSON(url)$items
      sel <- js[grepl("Educa", js$conta, ignore.case = TRUE) &
                grepl("LIQUIDADAS", toupper(js$coluna)), ]
      sum(suppressWarnings(as.numeric(sel$valor)), na.rm = TRUE)
    }, error = function(e) NA_real_)
    recs[[k]] <- data.frame(uf = u, ano = a, gasto_edu_total = val); k <- k + 1
    Sys.sleep(0.15)
  }
  do.call(rbind, recs)
}

# Monta o painel REAL (uf, ano, regiao, salario_medio, gasto_edu_pc, pib_pc),
# deflacionado pelo IPCA. Qualquer falha propaga => fallback à base calibrada.
build_panel_api <- function() {
  message("  [1/5] IPCA (BCB/SGS 433)...");        ipca <- get_ipca_index()
  message("  [2/5] PIB estadual (SIDRA ", TAB_PIB, ")...")
  pib <- get_sidra_uf(TAB_PIB, VAR_PIB, ANOS_PIB, "pib_mil")      # R$ mil correntes
  message("  [3/5] População (SIDRA ", TAB_POP, ")...")
  pop <- get_sidra_uf(TAB_POP, VAR_POP, anos, "pop")
  message("  [4/5] Rendimento (SIDRA ", TAB_RENDIMENTO, ")...")
  ren <- get_sidra_uf(TAB_RENDIMENTO, VAR_RENDIMENTO, anos, "rend_nom")
  message("  [5/5] Gasto em educação (SICONFI/RREO)... pode demorar")
  edu <- get_gasto_educ_uf(anos)                                 # R$ correntes

  d <- Reduce(function(x, y) merge(x, y, by = c("uf","ano"), all = TRUE),
              list(pib, pop, ren, edu))
  d$regiao <- regiao[d$uf]
  defl <- as.numeric(ipca[as.character(d$ano)])                  # deflator base 2010
  # valores reais (preços de 2010): divide o nominal pelo índice do ano
  d$pib_pc        <- (d$pib_mil * 1000 / d$pop) / defl           # R$/hab reais
  d$gasto_edu_pc  <- (d$gasto_edu_total / d$pop) / defl          # R$/hab reais
  d$salario_medio <- d$rend_nom / defl                          # R$ reais
  out <- d[, c("uf","ano","regiao","salario_medio","gasto_edu_pc","pib_pc")]
  out <- out[stats::complete.cases(out) & out$salario_medio > 0 &
             out$gasto_edu_pc > 0 & out$pib_pc > 0, ]
  if (nrow(out) < 50) stop("Painel oficial incompleto (", nrow(out), " obs).")
  message(sprintf("  Painel oficial montado: %d obs (%d UF, anos %d-%d).",
                  nrow(out), length(unique(out$uf)), min(out$ano), max(out$ano)))
  out[order(out$uf, out$ano), ]
}

painel <- NULL
if (use_api) {
  message("Coletando séries oficiais via API do Governo Federal...")
  for (p in c("jsonlite","sidrar"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  painel <- tryCatch(build_panel_api(), error = function(e) {
    message("!! Coleta oficial falhou: ", conditionMessage(e),
            "\n   Verifique conexão e os códigos de tabela (TAB_RENDIMENTO/SICONFI).",
            "\n   Usando a base CALIBRADA como fallback.")
    NULL
  })
  if (!is.null(painel))
    write.csv(painel, file.path(OUT, "painel_real_R.csv"), row.names = FALSE)
}

## ---- 2. Base CALIBRADA realista (usada se use_api=FALSE ou se a API falhar) --
if (is.null(painel)) {
base_elast <- c(N=0.20, NE=0.24, CO=0.16, SE=0.12, S=0.15)
true_elast <- sapply(ufs, function(u)
  max(min(base_elast[regiao[u]] + rnorm(1, 0, 0.06), 0.45), -0.05))
names(true_elast) <- ufs

reg_wage_premium <- c(N=-0.15, NE=-0.22, CO=0.08, SE=0.22, S=0.14)
reg_pib_premium  <- c(N=-0.12, NE=-0.20, CO=0.16, SE=0.20, S=0.14)
LG0 <- log(700); LP0 <- log(28000)

pib_level   <- sapply(ufs, function(u) exp(reg_pib_premium[regiao[u]] + rnorm(1,0,0.10)))
gasto_level <- runif(length(ufs), 0.80, 1.30); names(gasto_level) <- ufs
mu_uf <- sapply(ufs, function(u) reg_wage_premium[regiao[u]] + rnorm(1,0,0.06))
lam_t <- 0.015 * (anos - 2010) + rnorm(length(anos), 0, 0.01); names(lam_t) <- anos

rows <- list(); k <- 1
for (u in ufs) {
  g0 <- 600 * gasto_level[u]
  for (i in seq_along(anos)) {
    t <- anos[i]
    gasto <- g0 * 1.03^(i-1) * exp(rnorm(1,0,0.16))
    pib   <- 28000 * pib_level[u] * 1.018^(i-1) * exp(rnorm(1,0,0.05))
    lsal  <- 7.7 + mu_uf[u] + lam_t[as.character(t)] +
             true_elast[u]*(log(gasto)-LG0) + 0.35*(log(pib)-LP0) + rnorm(1,0,0.05)
    rows[[k]] <- data.frame(uf=u, ano=t, regiao=regiao[u],
                            salario_medio=exp(lsal), gasto_edu_pc=gasto, pib_pc=pib)
    k <- k + 1
  }
}
painel <- do.call(rbind, rows)
write.csv(painel, file.path(OUT, "painel_calibrado_R.csv"), row.names = FALSE)
fonte_dados <- "calibrada"
}  # ---- fim do bloco da base calibrada ----

# 'true_elast' (elasticidade verdadeira) existe apenas na base calibrada; com
# dados reais preenche-se com NA (coluna meramente diagnóstica no bootstrap).
if (!exists("true_elast")) true_elast <- setNames(rep(NA_real_, length(ufs)), ufs)
if (!exists("fonte_dados")) fonte_dados <- "oficial (API)"
cat(sprintf("Painel [%s]: %d obs\n", fonte_dados, nrow(painel)))

## ---- 3. Efeitos fixos bidirecionais + Driscoll-Kraay + Hausman ------------
painel$lsal   <- log(painel$salario_medio)
painel$lgasto <- log(painel$gasto_edu_pc)
painel$lpib   <- log(painel$pib_pc)
pd <- pdata.frame(painel, index = c("uf", "ano"))

fe  <- plm(lsal ~ lgasto + lpib, data = pd, model = "within", effect = "twoways")
re  <- plm(lsal ~ lgasto + lpib, data = pd, model = "random")
pol <- plm(lsal ~ lgasto + lpib, data = pd, model = "pooling")

se_cl <- coeftest(fe, vcov = vcovHC(fe, type = "HC0", cluster = "group"))
se_dk <- coeftest(fe, vcov = vcovSCC(fe, type = "HC0"))    # Driscoll-Kraay
haus  <- phtest(fe, re)
cat(sprintf("Elast. FE (clustered): %.4f (se %.4f)\n",
            se_cl["lgasto","Estimate"], se_cl["lgasto","Std. Error"]))
cat(sprintf("Driscoll-Kraay se: %.4f | Hausman p=%.3g\n",
            se_dk["lgasto","Std. Error"], haus$p.value))
elast_fe <- unname(coef(fe)["lgasto"])

## ---- 4. Elasticidades por macrorregião ------------------------------------
region_el <- lapply(c("N","NE","CO","SE","S"), function(R) {
  sub <- subset(painel, regiao == R); sub$trend <- sub$ano - 2010
  pdr <- pdata.frame(sub, index = c("uf","ano"))
  m <- plm(lsal ~ lgasto + lpib + trend, data = pdr, model = "within")
  s <- coeftest(m, vcov = vcovHC(m, cluster = "group"))
  data.frame(regiao=R, beta=coef(m)["lgasto"], se=s["lgasto","Std. Error"],
             p=s["lgasto","Pr(>|t|)"], n=nrow(sub))
})
region_el <- do.call(rbind, region_el)
write.csv(region_el, file.path(OUT, "region_elast_R.csv"), row.names = FALSE)
print(region_el)

## ---- 5. Bootstrap por UF + encolhimento EB --------------------------------
elast_uf_fun <- function(d) {
  X <- cbind(1, log(d$gasto_edu_pc), log(d$pib_pc), d$ano - 2010)
  coef(lm.fit(X, log(d$salario_medio)))[2]
}
boot_stat <- function(d, idx) elast_uf_fun(d[idx, ])
recs <- lapply(ufs, function(u) {
  sub <- subset(painel, uf == u)
  bo <- boot(sub, boot_stat, R = 2000)
  ci <- boot.ci(bo, type = "perc")$percent[4:5]
  data.frame(uf=u, regiao=regiao[u], est=mean(bo$t, na.rm=TRUE),
             lower=ci[1], upper=ci[2], true=true_elast[u])
})
boot_df <- do.call(rbind, recs)
boot_df$beta_reg <- region_el$beta[match(boot_df$regiao, region_el$regiao)]
boot_df$sd  <- (boot_df$upper - boot_df$lower) / (2 * 1.959964)
boot_df$var <- boot_df$sd^2
within <- tapply(seq_len(nrow(boot_df)), boot_df$regiao, function(ix)
  max(var(boot_df$est[ix]) - mean(boot_df$var[ix]), 1e-4))
boot_df$tau2 <- within[boot_df$regiao]
boot_df$w    <- boot_df$tau2 / (boot_df$tau2 + boot_df$var)
boot_df$elast_eb <- pmax(boot_df$w*boot_df$est + (1-boot_df$w)*boot_df$beta_reg, 0)
write.csv(boot_df[c("uf","regiao","est","lower","upper","true","beta_reg","w","elast_eb")],
          file.path(OUT, "elasticidades_uf_R.csv"), row.names = FALSE)
cat(sprintf("Elast média entre-UF: %.4f | EB média: %.4f\n",
            mean(boot_df$est), mean(boot_df$elast_eb)))

## ---- 6. Figuras (ggplot2) -------------------------------------------------
bd <- boot_df[order(boot_df$est), ]; bd$uf <- factor(bd$uf, levels = bd$uf)
g1 <- ggplot(bd, aes(est, uf, fill = regiao)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  geom_vline(xintercept = elast_fe, linetype = "dashed", color = "red") +
  scale_fill_manual(values = reg_col) +
  labs(x = "Elasticidade salário-educação", y = NULL,
       title = "Elasticidade por UF (IC bootstrap 95%)") + theme_minimal()
ggsave(file.path(FIG, "elasticidade_uf_ic_R.png"), g1, width = 7.2, height = 7.6, dpi = 200)

g2 <- ggplot(painel, aes(log(gasto_edu_pc), log(salario_medio), color = regiao)) +
  geom_point(alpha = 0.4, size = 1) + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~regiao, scales = "free") + scale_color_manual(values = reg_col) +
  labs(x = "log(Gasto Educ. pc)", y = "log(Salário)",
       title = "Relação salário-educação por macrorregião") + theme_minimal()
ggsave(file.path(FIG, "dispersao_uf_R.png"), g2, width = 11, height = 6.6, dpi = 200)

## ---- 7. Tabelas LaTeX -----------------------------------------------------
star <- function(p) ifelse(p<0.01,"***",ifelse(p<0.05,"**",ifelse(p<0.1,"*","")))
fmt  <- function(x, d=3) formatC(x, format="f", digits=d, decimal.mark=",")

con <- file(file.path(TAB, "region_elast_R.tex"), "w", encoding = "UTF-8")
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
  "Macrorregião & Elasticidade & Erro-padrão & $p$-valor & $N$ \\\\", "\\midrule"), con)
nomes <- c(N="Norte",NE="Nordeste",CO="Centro-Oeste",SE="Sudeste",S="Sul")
for (i in seq_len(nrow(region_el))) {
  r <- region_el[i, ]
  writeLines(sprintf("%s & %s%s & %s & %s & %d \\\\", nomes[r$regiao],
    fmt(r$beta), star(r$p), fmt(r$se), fmt(r$p), r$n), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), con); close(con)

cat("\nBloco empírico (R) concluído. Figuras e tabelas gerados.\n")
