# =========================================================
# 0) Paquetes y setup
# =========================================================
# install.packages(c("zoo","ggplot2","tseries","vars","urca","forecast","dplyr","tidyr"))
library(zoo)
library(ggplot2)
library(tseries)
library(vars)
library(urca)
library(forecast)
library(dplyr)
library(tidyr)

set.seed(42)

# =========================================================
# 1) Calendario y variables (72 trimestres 2007Q1–2024Q4)
# =========================================================
T  <- 72
fq <- 4
fechas <- as.yearqtr(seq(from = as.Date("2007-01-01"),
                         by = "quarter", length.out = T))
vars <- c("L_PP","L_PIB","L_GG","L_GC","L_GK","L_IF","L_IVA","L_ICE","L_IR")
k <- length(vars)

# =========================================================
# 2) DGP: I(1) con 2 relaciones de cointegración plausibles
# =========================================================
# β (columnas = vectores cointegrados)
beta <- matrix(0, nrow = k, ncol = 2, dimnames = list(vars, c("CI_gasto","CI_tributos")))
beta["L_GG","CI_gasto"] <- 1;   beta["L_GC","CI_gasto"] <- -0.7; beta["L_GK","CI_gasto"] <- -0.3
beta["L_IF","CI_tributos"] <- 1; beta["L_IVA","CI_tributos"] <- -0.5; beta["L_IR","CI_tributos"] <- -0.3; beta["L_ICE","CI_tributos"] <- -0.2

# α (velocidades de ajuste)
alpha <- matrix(0, nrow = k, ncol = 2, dimnames = list(vars, c("CI_gasto","CI_tributos")))
alpha["L_GG","CI_gasto"]  <- -0.25
alpha["L_GC","CI_gasto"]  <- -0.10
alpha["L_GK","CI_gasto"]  <- -0.05
alpha["L_IF","CI_tributos"]  <- -0.20
alpha["L_IVA","CI_tributos"] <- -0.10
alpha["L_IR","CI_tributos"]  <- -0.08
alpha["L_ICE","CI_tributos"] <- -0.05
alpha["L_PIB","CI_gasto"]    <- -0.02
alpha["L_PIB","CI_tributos"] <- -0.01
alpha["L_PP",] <- c(0,0)  # petróleo no corrige directamente

# Γ1 (dinámica de corto plazo sobre ΔY_{t-1})
Gamma1 <- matrix(0, nrow = k, ncol = k, dimnames = list(vars, vars))
Gamma1["L_PIB","L_PP"] <- 0.10
Gamma1["L_IF","L_PP"]  <- 0.08
Gamma1["L_IF","L_IVA"] <- 0.10; Gamma1["L_IF","L_IR"] <- 0.07; Gamma1["L_IF","L_ICE"] <- 0.05
Gamma1["L_GG","L_PIB"] <- 0.06
Gamma1["L_GC","L_GG"]  <- 0.10; Gamma1["L_GK","L_GG"] <- 0.06

# Σ (varianzas-covarianzas de shocks)
Sigma <- diag(c(0.20, 0.18, 0.15, 0.12, 0.12, 0.18, 0.15, 0.12, 0.12))
dimnames(Sigma) <- list(vars, vars)
Sigma["L_PIB","L_PP"] <- Sigma["L_PP","L_PIB"] <- 0.05
Sigma["L_IF","L_PP"]  <- Sigma["L_PP","L_IF"]  <- 0.04
Sigma["L_IF","L_IVA"] <- Sigma["L_IVA","L_IF"] <- 0.06
Sigma["L_GG","L_PIB"] <- Sigma["L_PIB","L_GG"] <- 0.04
Sigma <- (Sigma + t(Sigma))/2
C <- t(chol(Sigma))

# =========================================================
# 3) Simulación VECM (ΔY_t = Γ1 ΔY_{t-1} + α β' Y_{t-1} + ε_t)
# =========================================================
Y  <- matrix(0, nrow = k, ncol = T, dimnames = list(vars, NULL))
dY <- matrix(0, nrow = k, ncol = T, dimnames = list(vars, NULL))
Y[,1] <- c(4.6, 8.5, 9.2, 8.9, 7.8, 8.6, 7.9, 6.5, 7.6)

for (t in 2:T) {
  eps_t  <- C %*% rnorm(k)
  EC_lag <- t(beta) %*% Y[, t-1]           # (2x9)*(9x1) = (2x1)
  dY[,t] <- Gamma1 %*% dY[,t-1] + alpha %*% EC_lag + eps_t
  Y[,t]  <- Y[,t-1] + dY[,t]
}
