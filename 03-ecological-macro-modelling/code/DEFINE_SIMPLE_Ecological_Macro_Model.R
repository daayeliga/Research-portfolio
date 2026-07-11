# =============================================================
# DEFINE-SIMPLE-MODIFIED Ecological Macroeconomic Model
# with climate damages and Task 2 policy-mix scenarios
# Advanced Macroeconomics, MSc Development Economics, SOAS
# Author: David Ayuuya Ayeliga
# =============================================================

# ---- WORKING DIRECTORY + PACKAGES ----
# NOTE: update this path to your own local directory before running
setwd("C:/Users/ayeli/OneDrive/Desktop/Second Semester/Macro/Assessment")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
p_load(ggplot2, reshape2, ggthemes, showtext, sysfonts, showtextdb, xtable, readr, dplyr)

# Check current folder
getwd()

rm(list = ls(all = TRUE))
T <- 48
Scenarios <- 3

# =============================================================
# EQUATIONS
# =============================================================
Equations <- function(t){

  # --- Households ---
  Y_H[t]   <<- W[t] + DP[t] + BP[t] + int_D * D[t-1]
  W[t]     <<- s_W * Y[t]
  CO_PRI[t] <<- c_1 * Y_H[t-1] + c_2 * D[t-1]
  D[t]     <<- D[t-1] + Y_H[t] - CO_PRI[t]

  # --- Firms ---
  Y[t]     <<- CO_PRI[t] + I[t] + CO_GOV[t]
  TP_G[t]  <<- Y[t] - W[t] - int_C * L_C[t-1] - int_G * L_G[t-1] - delta[t] * K[t-1]
  TP[t]    <<- TP_G[t] - TAX_C[t]
  RP[t]    <<- s_F * TP[t]
  DP[t]    <<- TP[t] - RP[t]

  # Investment depends on lagged climate damages
  I[t]     <<- ((alpha_0 + alpha_1 * r[t-1]) * K[t-1] + delta[t] * K[t-1]) * (1 - D_T[t-1])
  r[t]     <<- TP[t] / K[t]
  I_G[t]   <<- beta[t] * I[t]
  beta[t]  <<- beta_0 - beta_1 * (int_G - int_C) + beta_2 * tucf[t-1]
  tucf[t]  <<- pucf + tau_C * omega
  I_C[t]   <<- I[t] - I_G[t]

  K_G[t]   <<- K_G[t-1] + I_G[t] - delta[t] * K_G[t-1]
  K_C[t]   <<- K_C[t-1] + I_C[t] - delta[t] * K_C[t-1]
  K[t]     <<- K_C[t] + K_G[t]

  L_G[t]   <<- L_G[t-1] + I_G[t] - beta[t] * RP[t] - delta[t] * K_G[t-1]
  L_C[t]   <<- L_C[t-1] + I_C[t] + I_G[t] - RP[t] - (L_G[t] - L_G[t-1]) - delta[t] * K[t-1]
  L[t]     <<- L_C[t] + L_G[t]

  # Depreciation depends on lagged climate damage
  delta[t] <<- delta_0 + (1 - delta_0) * (1 - ad_K) * D_T[t-1]

  # --- Banks ---
  BP[t]     <<- int_C * L_C[t-1] + int_G * L_G[t-1] + int_S * SEC[t-1] - int_D * D[t-1]
  SEC_red[t] <<- D[t] - L[t]

  # --- Government ---
  SEC[t]   <<- SEC[t-1] + int_S * SEC[t-1] - TAX_C[t] + CO_GOV[t]
  CO_GOV[t] <<- gov_C * Y[t-1]
  TAX_C[t] <<- tau_C * EMIS_F[t-1]

  # --- Emissions ---
  EMIS_F[t] <<- CI[t] * Y[t]
  CI[t]     <<- CI_max - ((CI_max - CI_min) / (1 + ci_1 * exp(-ci_2 * (K_G[t-1] / K_C[t-1]))))

  # Cumulative emissions and temperature
  CO2_CUM[t] <<- CO2_CUM[t-1] + EMIS_F[t]
  TEMP[t]    <<- (1 / (1 - f_nc)) * TCRE * CO2_CUM[t]

  # Climate damage rate
  D_T[t] <<- (1 - (1 / (1 + eta_1 * TEMP[t] + eta_2 * TEMP[t]^2 + eta_3 * TEMP[t]^6.754)))

  # --- Auxiliary equations ---
  Y_P[t]     <<- v * K[t]
  u[t]       <<- Y[t] / Y_P[t]
  g_Y[t]     <<- (Y[t] - Y[t-1]) / Y[t-1]
  lev[t]     <<- L[t] / K[t]
  SEC_Y[t]   <<- SEC[t] / Y[t]
  g_D[t]     <<- (D[t] - D[t-1]) / D[t-1]
  g_L[t]     <<- (L[t] - L[t-1]) / L[t-1]
  g_LG[t]    <<- (L_G[t] - L_G[t-1]) / L_G[t-1]
  g_LC[t]    <<- (L_C[t] - L_C[t-1]) / L_C[t-1]
  g_K[t]     <<- (K[t] - K[t-1]) / K[t-1]
  g_KG[t]    <<- (K_G[t] - K_G[t-1]) / K_G[t-1]
  g_KC[t]    <<- (K_C[t] - K_C[t-1]) / K_C[t-1]
  g_YH[t]    <<- (Y_H[t] - Y_H[t-1]) / Y_H[t-1]
  g_SEC[t]   <<- (SEC[t] - SEC[t-1]) / SEC[t-1]
  g_EMISF[t] <<- (EMIS_F[t] - EMIS_F[t-1]) / EMIS_F[t-1]
}

# =============================================================
# PARAMETERS AND INITIAL VALUES
# =============================================================
Parameters_InitialValues <- function(){

  c_2     <<- 0.05      # Propensity to consume out of deposits
  int_C   <<- 0.09      # Interest rate on conventional loans
  int_G   <<- 0.09      # Interest rate on green loans
  int_D   <<- 0.04      # Interest rate on deposits
  s_W     <<- 0.55      # Wage share
  alpha_1 <<- 0.1       # Sensitivity of investment rate to profit rate
  beta_1  <<- -1        # Sensitivity of green investment share to interest rate differential
  beta_2  <<- 20        # Sensitivity of green investment share to carbon tax
  delta_0 <<- 0.05      # Autonomous depreciation rate
  int_S   <<- 0.0254    # Interest rate on government securities
  pucf    <<- 0.021     # Pre-tax unit cost of producing fossil energy (USD trillion/EJ)
  omega   <<- 0.072     # CO2 intensity of fossil energy (GtCO2/EJ)
  ci_2    <<- 3.579244  # Parameter linking K_G/K_C to carbon intensity
  CI_max  <<- 0.6       # Max carbon intensity
  CI_min  <<- 0.05      # Min carbon intensity

  # Climate/temperature parameters (Matthews et al. 2021)
  f_nc <<- 0.14         # Non-CO2 fraction of total anthropogenic forcing
  TCRE <<- 0.44 / 1000  # TCRE (degC/GtCO2)

  # Damage function parameters (Dietz and Stern 2015)
  eta_1 <<- 0
  eta_2 <<- 0.00284
  eta_3 <<- 0.0000819

  # Adaptation parameter
  ad_K <<- 0.8          # Fraction of gross damages to capital avoided through adaptation

  # Initial growth and macro ratios
  g_Y[1]    <<- 0.028      # Growth rate of output
  I[1]      <<- 0.26 * Y[1]  # Investment
  L[1]      <<- 0.96 * Y[1]  # Total loans
  u[1]      <<- 0.73       # Capacity utilisation rate
  Y[1]      <<- 106.17     # Output (USD trillion)
  I_G[1]    <<- 1.7        # Green investment (USD trillion)
  SEC_Y[1]  <<- 0.94       # Government debt-to-GDP ratio
  TAX_C[1]  <<- 0.104      # Carbon tax revenues (USD trillion)
  EMIS_F[1] <<- 37.8       # Fossil CO2 emissions (GtCO2)

  g_D[1]     <<- g_Y[1]
  g_L[1]     <<- g_Y[1]
  g_K[1]     <<- g_Y[1]
  g_YH[1]    <<- g_Y[1]
  g_LG[1]    <<- g_Y[1]
  g_LC[1]    <<- g_Y[1]
  g_KG[1]    <<- g_Y[1]
  g_KC[1]    <<- g_Y[1]
  g_SEC[1]   <<- g_Y[1]
  g_EMISF[1] <<- g_Y[1]

  # ------------------- Model-constrained -------------------
  Y_H[1]  <<- W[1] + DP[1] + BP[1] + int_D * (D[1] / (1 + g_D[1]))
  W[1]    <<- s_W * Y[1]
  c_1     <<- (CO_PRI[1] - c_2 * D[1] / (1 + g_D[1])) / (Y_H[1] / (1 + g_YH[1]))
  CO_PRI[1] <<- Y[1] - I[1] - CO_GOV[1]
  TP_G[1] <<- Y[1] - W[1] - int_C * (L_C[1] / (1 + g_LC[1])) - int_G * (L_G[1] / (1 + g_LG[1])) - delta[1] * K[1] / (1 + g_K[1])
  TP[1]   <<- TP_G[1] - TAX_C[1]
  s_F     <<- RP[1] / TP[1]
  DP[1]   <<- TP[1] - RP[1]
  r[1]    <<- TP[1] / K[1]
  beta[1] <<- I_G[1] / I[1]
  beta_0  <<- beta[1] + beta_1 * (int_G - int_C) - beta_2 * tucf[1]
  I_C[1]  <<- I[1] - I_G[1]
  K_G[1]  <<- ((1 + g_KG[1]) / (g_KG[1] + delta[1])) * I_G[1]
  K_C[1]  <<- ((1 + g_KC[1]) / (g_KC[1] + delta[1])) * I_C[1]
  K[1]    <<- K_G[1] + K_C[1]
  L_G[1]  <<- ((1 + g_LG[1]) / g_LG[1]) * (I_G[1] - beta[1] * RP[1] - delta[1] * K_G[1] / (1 + g_KG[1]))
  RP[1]   <<- I[1] - g_L[1] / (1 + g_L[1]) * L[1] - delta[1] * K[1] / (1 + g_K[1])

  # Temperature and damages initial conditions
  TEMP[1]    <<- 1.4  # Atmospheric temperature (degC)
  CO2_CUM[1] <<- TEMP[1] * (1 - f_nc) / (TCRE)  # Cumulative CO2 emissions (GtCO2)
  D_T[1]     <<- (1 - (1 / (1 + eta_1 * TEMP[1] + eta_2 * TEMP[1]^2 + eta_3 * TEMP[1]^6.754)))

  # Depreciation rate from climate damage
  delta[1] <<- delta_0 + (1 - delta_0) * (1 - ad_K) * D_T[1]

  L_C[1]     <<- L[1] - L_G[1]
  BP[1]      <<- int_C * (L_C[1] / (1 + g_LC[1])) + int_G * (L_G[1] / (1 + g_LG[1])) + int_S * (SEC[1] / (1 + g_SEC[1])) - int_D * (D[1] / (1 + g_D[1]))
  D[1]       <<- L[1] + SEC[1]
  SEC_red[1] <<- SEC[1]
  gov_C      <<- (((g_SEC[1] - int_S) / (1 + g_SEC[1])) * SEC[1] + TAX_C[1]) / (Y[1] / (1 + g_Y[1]))
  tau_C      <<- TAX_C[1] / (EMIS_F[1] / (1 + g_EMISF[1]))
  CO_GOV[1]  <<- gov_C * (Y[1] / (1 + g_Y[1]))
  CI[1]      <<- EMIS_F[1] / Y[1]
  ci_1       <<- ((CI_max - CI_min) / (CI_max - CI[1]) - 1) * exp(ci_2 * ((K_G[1] / (1 + g_KG[1])) / ((K_C[1]) / (1 + g_KC[1]))))
  v          <<- Y_P[1] / K[1]
  Y_P[1]     <<- Y[1] / u[1]
  lev[1]     <<- L[1] / K[1]
  SEC[1]     <<- SEC_Y[1] * Y[1]
  tucf[1]    <<- pucf + tau_C * omega

  # Alpha_0 consistent with damage-adjusted investment
  alpha_0 <<- I[1] / ((1 - D_T[1]) * (K[1] / (1 + g_K[1]))) - alpha_1 * r[1] - delta[1]
}

# =============================================================
# PARAMETERS FOR SCENARIOS (TASK 2 POLICY MIX)
# =============================================================

# Scenario 2: Higher carbon tax and higher government spending
Parameters_scenario2 <- function(){
  tau_C <<- 0.08   # Carbon tax (trillion USD per GtCO2)
  gov_C <<- 0.02   # Government consumption-to-GDP ratio
}

# Scenario 3: Degrowth and higher carbon tax
Parameters_scenario3 <- function(){
  c_1     <<- 0.8   # Lower propensity to consume out of disposable income
  alpha_0 <<- 0.02  # Lower autonomous investment spending
  tau_C   <<- 0.08  # Higher carbon tax
}

# =============================================================
# SYMBOLS AND VECTORS FOR ENDOGENOUS VARIABLES/PARAMETERS
# =============================================================

# List of Greek letters
greek_letters <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
                    "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho",
                    "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega")

# Symbols and vectors for endogenous variables
code_body  <- body(Equations)
code_lines <- as.character(code_body)
Symbols_var <- regmatches(code_lines, regexpr("^[^\\[]+", code_lines))
Symbols_var <- sub("\\[t\\] <<-.*", "", Symbols_var)
Symbols_var <- Symbols_var[!Symbols_var %in% c("{", "}") & Symbols_var != ""]
for (v in 1:length(Symbols_var)){
  assign(noquote(paste(Symbols_var[v])), vector(length = T))
}
Symbols_var_sort <- sort(Symbols_var)
is_greek <- grepl(paste(greek_letters, collapse = "|"), Symbols_var_sort, ignore.case = TRUE)
Symbols_var_sort <- c(sort(Symbols_var_sort[!is_greek]), Symbols_var_sort[is_greek])

# Symbols for parameters
code_body  <- body(Parameters_InitialValues)
code_lines <- as.character(code_body)
Symbols_par <- regmatches(code_lines, regexpr("^[^\\[]+<<-", code_lines))
Symbols_par <- sub("<<-$", "", Symbols_par)
Symbols_par <- Symbols_par[Symbols_par != ""]

Symbols_par_sort <- sort(Symbols_par)
is_greek <- grepl(paste(greek_letters, collapse = "|"), Symbols_par_sort, ignore.case = TRUE)
Symbols_par_sort <- c(sort(Symbols_par_sort[!is_greek]), Symbols_par_sort[is_greek])

# =============================================================
# RUN THE MODEL AND SAVE SCENARIO RESULTS
# =============================================================
for (j in 1:Scenarios){

  for (iterations in 1:10){ Parameters_InitialValues() }

  for (i in 2:T){
    if (j == 2 & i >= 5) { Parameters_scenario2() }
    if (j == 3 & i >= 5) { Parameters_scenario3() }
    for (iterations in 1:5){ Equations(t = i) }
  }

  # Tables with variables in same order as equations
  tablename <- paste("Variables", j, sep = "")
  values <- setNames(data.frame(matrix(ncol = length(Symbols_var), nrow = T)), Symbols_var)
  for (v in 1:length(Symbols_var)){
    values[, v] <- eval(parse(text = (paste(noquote(Symbols_var[v])))))
  }
  assign(tablename, round(values, 4))

  # Tables with variables sorted alphabetically
  tablename <- paste("Variables_sorted", j, sep = "")
  values <- setNames(data.frame(matrix(ncol = length(Symbols_var_sort), nrow = T)), Symbols_var_sort)
  for (v in 1:length(Symbols_var_sort)){
    values[, v] <- eval(parse(text = (paste(noquote(Symbols_var_sort[v])))))
  }
  assign(tablename, round(values, 6))

  # Table with parameters sorted alphabetically
  tablename <- paste("Parameters_sorted", j, sep = "")
  values <- setNames(data.frame(matrix(ncol = length(Symbols_par_sort), nrow = 1)), Symbols_par_sort)
  for (v in 1:length(Symbols_par_sort)){
    values[, v] <- eval(parse(text = (paste(noquote(Symbols_par_sort[v])))))
  }
  assign(tablename, round(values, 4))
  rm(values)

  # Save as CSV
  write.csv(get(paste0("Variables", j)), paste0("Variables", j, ".csv"))
  write.csv(get(paste0("Variables_sorted", j)), paste0("Variables-sorted", j, ".csv"))
  write.csv(get(paste0("Parameters_sorted", j)), paste0("Parameters-sorted", j, ".csv"))
}

# =============================================================
# CHART INFO
# =============================================================
library(ggplot2)
library(reshape2)
library(ggthemes)
library(showtext)
library(sysfonts)
library(showtextdb)

showtext_auto()

Colour1 <- "#000000"
Colour2 <- "#009E73"
Colour3 <- "#D55E00"
Linetype1 <- "solid"
Linetype2 <- "solid"
Linetype3 <- "solid"

pch1 <- 26
pch2 <- 26
pch3 <- 26

Size1 <- 0.8
Size2 <- 0.8
Size3 <- 0.8

Scenario1 <- "Baseline"
Scenario2 <- "Higher carbon tax + higher gov. spending"
Scenario3 <- "Degrowth + higher carbon tax"

Periods <- 1:T

Chart3 <- function(fig_title, variable, label, min, max, sce1, sce2, sce3, legend_x, legend_y){

  eval(parse(text = (paste0("Line1=Variables", sce1, "[,c('", variable, "')]"))))
  eval(parse(text = (paste0("Line2=Variables", sce2, "[,c('", variable, "')]"))))
  eval(parse(text = (paste0("Line3=Variables", sce3, "[,c('", variable, "')]"))))

  dataframe <<- data.frame(Periods, Line1, Line2, Line3)
  data_long <<- melt(dataframe, id = "Periods")

  plot <- ggplot(data = data_long, aes(x = Periods, y = value, colour = variable,
                                        linetype = variable, size = variable, shape = variable)) +
    geom_line(aes(linetype = variable, size = variable)) +
    geom_point(aes(shape = variable)) +
    scale_color_manual(name = NULL,
                        labels = c(eval(parse(text = (paste0("Scenario", sce1)))),
                                   eval(parse(text = (paste0("Scenario", sce2)))),
                                   eval(parse(text = (paste0("Scenario", sce3))))),
                        values = c(Colour1, Colour2, Colour3)) +
    scale_linetype_manual(name = NULL,
                          labels = c(eval(parse(text = (paste0("Scenario", sce1)))),
                                     eval(parse(text = (paste0("Scenario", sce2)))),
                                     eval(parse(text = (paste0("Scenario", sce3))))),
                          values = c(Linetype1, Linetype2, Linetype3)) +
    scale_size_manual(name = NULL,
                       labels = c(eval(parse(text = (paste0("Scenario", sce1)))),
                                  eval(parse(text = (paste0("Scenario", sce2)))),
                                  eval(parse(text = (paste0("Scenario", sce3))))),
                       values = c(Size1, Size2, Size3)) +
    scale_shape_manual(name = NULL,
                        labels = c(eval(parse(text = (paste0("Scenario", sce1)))),
                                   eval(parse(text = (paste0("Scenario", sce2)))),
                                   eval(parse(text = (paste0("Scenario", sce3))))),
                        values = as.numeric(c(pch1, pch2, pch3))) +
    theme_bw() +
    theme(axis.title = element_text(size = 31, colour = "black"),
          axis.text = element_text(size = 30, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = c(legend_x, legend_y),
          legend.title = element_blank(),
          legend.text = element_text(size = 30, colour = "black"),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0, "cm"),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.border = element_rect(color = "black", size = 0.8),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey")) +
    labs(y = label, x = "Year") +
    geom_rect(aes(xmin = 1, xmax = T, ymin = min, ymax = max), fill = "transparent") +
    coord_cartesian(ylim = c(min, max)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(T - 40, T - 30, T - 20, T - 10, T),
                        labels = c(2030, 2040, 2050, 2060, 2070)) +
    scale_y_continuous(expand = c(0, 0))

  ggsave(fig_title, plot, width = 6, height = 4, dpi = 300)
}

# =============================================================
# CREATE CHARTS (including temperature, depreciation, damage)
# =============================================================
Chart3(fig_title = "Fig1-Leverage.jpeg", variable = 'lev', label = "Firm leverage",
       min = 0, max = 0.6, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig2-Growth.jpeg", variable = 'g_Y', label = "Growth rate of output",
       min = -0.02, max = 0.06, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig3-Output.jpeg", variable = 'Y', label = "Output (trillion USD)",
       min = 100, max = 400, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig4-Profit-rate.jpeg", variable = 'r', label = "Rate of profit",
       min = 0, max = 0.08, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig5-Share-of-green-investment.jpeg", variable = 'beta', label = "Share of green investment",
       min = 0, max = 0.4, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig6-Carbon-intensity.jpeg", variable = 'CI', label = expression("CO"[2] * " intensity"),
       min = 0, max = 0.5, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig7-Emissions.jpeg", variable = 'EMIS_F', label = expression("Fossil CO"[2] * " emissions (Gt)"),
       min = 0, max = 150, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig8-Government-debt.jpeg", variable = 'SEC_Y', label = expression("Government debt-to-GDP ratio"),
       min = 0, max = 1.4, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig9-Temperature.jpeg", variable = 'TEMP', label = "Atmospheric temperature (degC)",
       min = 1, max = 4, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig10-Depreciation.jpeg", variable = 'delta', label = "Depreciation rate",
       min = 0, max = 0.1, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)

Chart3(fig_title = "Fig11-Damage-rate.jpeg", variable = 'D_T', label = "Climate damage rate",
       min = 0, max = 0.4, sce1 = 1, sce2 = 2, sce3 = 3, legend_x = 0.25, legend_y = 0.92)
