# =============================================================================
# BLOCO DSGE (porte R do 02_dsge_capital_humano.py) — RStudio
# Modelo neoclássico com capital físico e humano, tributação distorciva e gasto
# público em educação como insumo produtivo da acumulação de capital humano.
# Solução por previsão perfeita não linear (stacked-time Newton) com estado
# estacionário terminal endógeno. Cenários de corte 1/5/10/25% + teto (EC 95),
# custo fiscal em R$, efeito Laffer, bem-estar, multiplicadores e desigualdade.
#
# Referência numérica validada: o porte Python (02_...py). Este porte R
# reproduz a mesma lógica para uso no RStudio / publicação no GitHub.
# =============================================================================

pkgs <- c("nleqslv", "ggplot2", "dplyr", "tidyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

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

## ---- Parâmetros e estado estacionário -------------------------------------
P <- list(alpha=0.40, beta=0.96, dk=0.10, dh=0.04, tau=0.33, theta=1.75,
          phi=0.30, gsh=0.05, rhoG=0.90)
a<-P$alpha; b<-P$beta; dk<-P$dk; dh<-P$dh; tau<-P$tau; th<-P$theta; phi<-P$phi
PIB_BR_BI <- 10900; DISC <- 0.96

steady_state <- function(P) {
  Rtil <- (1/b - 1 + dk)/(1 - tau); KY <- a/Rtil
  coef <- (a/Rtil)^(a/(1-a)); cshare <- 1 - dk*KY - P$gsh
  M <- (1-tau)*(1-a)/cshare; l <- M/(th+M)
  Y <- coef*l; K <- KY*Y; Ge <- P$gsh*Y; Cc <- cshare*Y
  A <- dh/(Ge^phi)
  list(H=1, K=K, l=l, Y=Y, Ge=Ge, C=Cc, A=A, I=dk*K, S=(1-a)*Y/l, T=tau*Y,
       r=a*Y/K, KY=KY)
}
ss <- steady_state(P); A <- ss$A
Kss<-ss$K; Css<-ss$C; lss<-ss$l
Rtil <- (1/b-1+dk)/(1-tau); KY <- a/Rtil; Bcoef <- (a/Rtil)^(a/(1-a))
KAPPA <- PIB_BR_BI/ss$Y

steady_state_for_Ge <- function(Ge) {              # EE terminal (nível fixo)
  H <- A*Ge^phi/dh
  l <- ((1-tau)*(1-a)*Bcoef*H + th*Ge) /
       (Bcoef*H*(th*(1-dk*KY) + (1-tau)*(1-a)))
  Y <- Bcoef*H*l; K <- KY*Y
  list(H=H, K=K, l=l, Y=Y, Ge=Ge, C=Y*(1-dk*KY)-Ge, I=dk*K, S=(1-a)*Y/l, T=tau*Y)
}

TT <- 220
Yfun <- function(K,H,l) K^a*(H*l)^(1-a)
human_capital_path <- function(Ge) {
  H <- numeric(TT+2); H[1] <- ss$H
  for (t in 1:(TT+1)) H[t+1] <- (1-dh)*H[t] + A*(Ge[min(t,TT+1)]^phi)
  H
}

solve_path <- function(Ge) {
  Hp <- human_capital_path(Ge)
  tss <- steady_state_for_Ge(Ge[TT+1]); Cterm<-tss$C; lterm<-tss$l
  resid <- function(x) {
    Cc<-x[1:TT]; l<-x[(TT+1):(2*TT)]; Kn<-x[(2*TT+1):(3*TT)]
    K<-c(Kss, Kn); Y<-Yfun(K[1:TT], Hp[1:TT], l)
    r1 <- th/(1-l) - (1-tau)*(1-a)*Y/(Cc*l)
    r2 <- Cc + Kn - (1-dk)*K[1:TT] + Ge[1:TT] - Y
    l_n<-c(l[-1], lterm); C_n<-c(Cc[-1], Cterm)
    Yn<-Yfun(Kn, Hp[2:(TT+1)], l_n); Rn<-(1-tau)*a*Yn/Kn + (1-dk)
    r3 <- 1/Cc - b/C_n*Rn
    c(r1, r2, r3)
  }
  x0 <- c(rep(Cterm,TT), rep(lterm,TT), rep(tss$K,TT))
  sol <- nleqslv::nleqslv(x0, resid, method="Newton",
                          control=list(maxit=500, ftol=1e-11))
  Cc<-sol$x[1:TT]; l<-sol$x[(TT+1):(2*TT)]; Kn<-sol$x[(2*TT+1):(3*TT)]
  K<-c(Kss,Kn); Y<-Yfun(K[1:TT], Hp[1:TT], l)
  I<-Kn-(1-dk)*K[1:TT]; S<-(1-a)*Y/l; Tax<-tau*Y
  dev <- function(x,x0) 100*(x/x0-1)
  devs <- data.frame(t=0:(TT-1),
    Produto=dev(Y,ss$Y), Investimento=dev(I,ss$I), Salario=dev(S,ss$S),
    Consumo=dev(Cc,ss$C), Arrecadacao=dev(Tax,ss$T),
    CapitalHumano=dev(Hp[1:TT],ss$H), GastoEduc=dev(Ge[1:TT],ss$Ge))
  levels <- data.frame(t=0:(TT-1), Y=Y, I=I, S=S, C=Cc, Tax=Tax,
                       H=Hp[1:TT], Ge=Ge[1:TT], l=l, K=K[1:TT])
  list(dev=devs, lev=levels, ok=sol$termcd %in% c(1,2),
       res=max(abs(sol$fvec)))
}

path_permanent <- function(cut) rep(ss$Ge*(1-cut), TT+1)
path_ceiling <- function(total=0.10, ramp=8)
  ss$Ge*(1 - total*pmin((0:TT)/ramp, 1))

## ---- Cenários e métricas ---------------------------------------------------
CUTS <- c(0.01, 0.05, 0.10, 0.25)
results <- lapply(CUTS, function(cut) solve_path(path_permanent(cut)))
names(results) <- as.character(CUTS)
for (cut in CUTS) {
  d <- results[[as.character(cut)]]$dev
  cat(sprintf("Corte %2d%%: H_LP=%.3f%% Y_LP=%.3f%%\n", cut*100,
              tail(d$CapitalHumano,1), tail(d$Produto,1)))
}
ec <- solve_path(path_ceiling(0.10, 8))

cum_disc <- function(x, n) -sum(x[1:n]*DISC^(0:(n-1)))
welfare_cev <- function(lev) {
  u <- log(lev$C) + th*log(1-lev$l); disc <- b^(0:(TT-1))
  W <- sum(disc*u) + b^TT * tail(u,1)/(1-b)
  W_ss <- (log(Css)+th*log(1-lss))/(1-b)
  -(exp((W-W_ss)*(1-b))-1)*100
}
multiplier_cum <- function(lev, n) {
  dY<-lev$Y-ss$Y; dGe<-lev$Ge-ss$Ge; w<-DISC^(0:(n-1))
  sum(w*dY[1:n])/sum(w*dGe[1:n])
}
fiscal <- function(lev, n) {
  dGe<-(lev$Ge-ss$Ge)*KAPPA; dTax<-(lev$Tax-ss$T)*KAPPA
  sav<--dGe; los<--dTax; w<-DISC^(0:(n-1))
  ps<-sum(w*sav[1:n]); pl<-sum(w*los[1:n])
  list(pv_save=ps, pv_loss=pl, net=ps-pl, eroded=100*pl/ps)
}

## ---- Figura: múltiplos cortes (ggplot2) -----------------------------------
long <- do.call(rbind, lapply(CUTS, function(cut) {
  d <- results[[as.character(cut)]]$dev
  d <- d[d$t < 40, c("t","Produto","Salario","Arrecadacao","CapitalHumano")]
  d <- tidyr::pivot_longer(d, -t, names_to="Variavel", values_to="Desvio")
  d$corte <- paste0(cut*100, "%"); d
}))
gA <- ggplot(long, aes(t, Desvio, color=corte)) + geom_line(linewidth=0.9) +
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~Variavel, scales="free_y") +
  labs(x="Anos", y="Desvio % do EE", color="Corte",
       title="Trajetórias por magnitude do corte permanente") + theme_minimal()
ggsave(file.path(FIG, "irf_multi_corte_R.png"), gA, width=11, height=7, dpi=200)

## ---- Tabela: cenários -----------------------------------------------------
fmt <- function(x,d=3) formatC(x, format="f", digits=d, decimal.mark=",")
con <- file(file.path(TAB, "cenarios_corte_R.tex"), "w", encoding="UTF-8")
writeLines(c("\\begin{tabular}{lrrrrr}", "\\toprule",
  "Corte & $\\Delta H_{LP}$ & $\\Delta Y_{LP}$ & Bem-estar & Mult.\\ 10a & \\% erodido 20a \\\\",
  "\\midrule"), con)
for (cut in CUTS) {
  d<-results[[as.character(cut)]]$dev; lv<-results[[as.character(cut)]]$lev
  fi<-fiscal(lv,20)
  writeLines(sprintf("%d\\%% & %s & %s & %s & %s & %s\\%% \\\\", cut*100,
    fmt(tail(d$CapitalHumano,1)), fmt(tail(d$Produto,1)),
    fmt(welfare_cev(lv)), fmt(multiplier_cum(lv,10),2), fmt(fi$eroded,1)), con)
}
writeLines(c("\\bottomrule","\\end{tabular}"), con); close(con)

cat("\nBloco DSGE (R) concluído. Figuras e tabelas gerados.\n")
