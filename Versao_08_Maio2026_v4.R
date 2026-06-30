################################################################################
# DISSERTAÇÃO - PAINEL SSM12: VC_AI, NPL, FE DINÂMICO, DIFF-GMM E SYS-GMM
################################################################################

rm(list = ls())

DATA_DIR <- "C:/Dissertação/Dados"
OUT_DIR <- file.path(dirname(DATA_DIR), "resultados_R")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

TMP_DIR <- file.path(OUT_DIR, "tmp_R")
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = TMP_DIR, TEMP = TMP_DIR, TMP = TMP_DIR)
options(openxlsx.tempdir = TMP_DIR)

# ==============================================================================
# 1. PACOTES
# ==============================================================================
pkgs <- c(
  "readxl", "dplyr", "plm", "lmtest", "sandwich",
  "lubridate", "ggplot2", "openxlsx", "tibble", "tidyr"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

lag <- dplyr::lag
select <- dplyr::select
filter <- dplyr::filter

# ==============================================================================
# 2. PAÍSES DA AMOSTRA
# ==============================================================================
SSM12 <- c(
  "Austria", "Belgium", "Germany", "Greece", "Spain", "Finland",
  "France", "Ireland", "Italy", "Luxembourg", "Netherlands", "Portugal"
)

ISO3 <- c(
  AUT = "Austria", BEL = "Belgium", DEU = "Germany", GRC = "Greece",
  ESP = "Spain", FIN = "Finland", FRA = "France", IRL = "Ireland",
  ITA = "Italy", LUX = "Luxembourg", NLD = "Netherlands", PRT = "Portugal"
)

normalizar_pais <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% names(ISO3) ~ ISO3[x],
    grepl("Austria", x) ~ "Austria",
    grepl("Belgium", x) ~ "Belgium",
    grepl("Germany", x) ~ "Germany",
    grepl("Greece", x) ~ "Greece",
    grepl("Spain", x) ~ "Spain",
    grepl("Finland", x) ~ "Finland",
    grepl("France", x) ~ "France",
    grepl("Ireland", x) ~ "Ireland",
    grepl("Italy", x) ~ "Italy",
    grepl("Luxembourg", x) ~ "Luxembourg",
    grepl("Netherlands", x) ~ "Netherlands",
    grepl("Portugal", x) ~ "Portugal",
    TRUE ~ NA_character_
  )
}

# ==============================================================================
# 3. FUNÇÕES AUXILIARES
# ==============================================================================
read_ecb <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet)

  raw %>%
    transmute(
      date = as.Date(DATE),
      country = normalizar_pais(`REFERENCE AREA (DESC.)`),
      value = suppressWarnings(as.numeric(OBS.VALUE))
    ) %>%
    filter(!is.na(country)) %>%
    mutate(year = lubridate::year(date)) %>%
    filter(!is.na(value), !is.na(year), year >= 2015, year <= 2023)
}

anualizar <- function(df, nome) {
  df %>%
    group_by(country, year) %>%
    summarise(!!nome := mean(value, na.rm = TRUE), .groups = "drop")
}

safe_vcov_fe <- function(model) {
  tryCatch(
    vcovSCC(model, type = "HC1", maxlag = 2),
    error = function(e) vcovHC(model, type = "HC1", cluster = "group")
  )
}

tidy_plm <- function(model, nome) {
  ct <- coeftest(model, vcov = safe_vcov_fe(model))

  data.frame(
    modelo = nome,
    variavel = rownames(ct),
    coeficiente = ct[, 1],
    erro_padrao = ct[, 2],
    estatistica = ct[, 3],
    p_valor = ct[, 4],
    row.names = NULL
  )
}

tidy_gmm <- function(model, nome) {
  sm <- summary(model, robust = TRUE)
  cf <- sm$coefficients

  data.frame(
    modelo = nome,
    variavel = rownames(cf),
    coeficiente = cf[, 1],
    erro_padrao = cf[, 2],
    estatistica = cf[, 3],
    p_valor = cf[, 4],
    row.names = NULL
  )
}

extrair_diagnosticos_gmm <- function(model, nome) {
  ar1 <- tryCatch(mtest(model, order = 1), error = function(e) NULL)
  ar2 <- tryCatch(mtest(model, order = 2), error = function(e) NULL)
  sg <- tryCatch(sargan(model), error = function(e) NULL)

  data.frame(
    modelo = nome,
    teste = c("Arellano-Bond AR(1)", "Arellano-Bond AR(2)", "Sargan/Hansen"),
    estatistica = c(
      if (!is.null(ar1)) unname(ar1$statistic) else NA_real_,
      if (!is.null(ar2)) unname(ar2$statistic) else NA_real_,
      if (!is.null(sg)) unname(sg$statistic) else NA_real_
    ),
    p_valor = c(
      if (!is.null(ar1)) ar1$p.value else NA_real_,
      if (!is.null(ar2)) ar2$p.value else NA_real_,
      if (!is.null(sg)) sg$p.value else NA_real_
    ),
    interpretacao = c(
      "Espera-se rejeitar AR(1) em diferenças.",
      "Não deve rejeitar AR(2); se rejeitar, instrumentos problemáticos.",
      "Não deve rejeitar; p muito alto também pode indicar instrumentos demais."
    )
  )
}

extrair_teste <- function(teste, modelo, nome_teste, decisao_5pct = TRUE) {
  if (is.null(teste)) {
    return(data.frame(
      modelo = modelo,
      teste = nome_teste,
      estatistica = NA_real_,
      p_valor = NA_real_,
      resultado_5pct = "Teste não calculado",
      row.names = NULL
    ))
  }

  pvalor <- suppressWarnings(as.numeric(teste$p.value))
  stat <- suppressWarnings(as.numeric(teste$statistic[1]))

  data.frame(
    modelo = modelo,
    teste = nome_teste,
    estatistica = stat,
    p_valor = pvalor,
    resultado_5pct = ifelse(
      is.na(pvalor),
      "Sem p-valor",
      ifelse(pvalor < 0.05, "Rejeita H0", "Não rejeita H0")
    ),
    row.names = NULL
  )
}

rodar_teste <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

# ==============================================================================
# 4. DADOS ECB
# ==============================================================================
cat("Carregando dados ECB...\n")

npl <- anualizar(
  read_ecb(file.path(DATA_DIR, "NPL_por_pais_ECB_CBD2.xlsx"), "DATA(CBD2)"),
  "npl"
)

nim <- anualizar(
  read_ecb(file.path(DATA_DIR, "NIM_por_pais_ECB_SUP.xlsx"), "DATA(SUP)"),
  "nim"
)

loans <- read_ecb(
  file.path(DATA_DIR, "Loans_NFCs_por_pais_ECB_BSI.xlsx"),
  "DATA(BSI)"
) %>%
  group_by(country, year) %>%
  summarise(loans = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(credit_growth = 100 * (loans / lag(loans) - 1)) %>%
  ungroup()

# ==============================================================================
# 5. VC_AI OECD
# ==============================================================================
cat("Carregando VC_AI OECD...\n")

vc_raw <- read_excel(
  file.path(DATA_DIR, "VC_AI_por_pais_OECD.xlsx"),
  sheet = "data",
  col_names = TRUE
)

vc_ai <- vc_raw %>%
  mutate(
    country = normalizar_pais(Country),
    year = suppressWarnings(as.integer(Year)),
    vc_ai = suppressWarnings(as.numeric(Sum_of_deals))
  ) %>%
  filter(!is.na(country), !is.na(year), !is.na(vc_ai),
         year >= 2015, year <= 2023) %>%
  group_by(country, year) %>%
  summarise(vc_ai = sum(vc_ai, na.rm = TRUE), .groups = "drop") %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    log_vc_ai = log(vc_ai + 1),
    lag_vc_ai = lag(vc_ai),
    lag_log_vc_ai = lag(log_vc_ai)
  ) %>%
  ungroup()

cat("Obs VC_AI:", nrow(vc_ai), "\n")

# ==============================================================================
# 6. PAINEL FINAL
# ==============================================================================
panel <- expand.grid(country = SSM12, year = 2015:2023,
                     stringsAsFactors = FALSE) %>%
  left_join(npl, by = c("country", "year")) %>%
  left_join(nim, by = c("country", "year")) %>%
  left_join(loans, by = c("country", "year")) %>%
  left_join(vc_ai, by = c("country", "year")) %>%
  mutate(covid = as.integer(year %in% c(2020, 2021))) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    lag_npl = lag(npl),
    lag_nim = lag(nim),
    lag_credit = lag(credit_growth)
  ) %>%
  ungroup()

cat("Painel bruto:", nrow(panel), "\n")
print(colSums(is.na(panel)))

write.csv(panel, "painel_final_vc_ai.csv", row.names = FALSE)

# ==============================================================================
# 7. BASES DOS MODELOS
# ==============================================================================
panel_base <- panel %>%
  select(country, year, npl, nim, credit_growth, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth))

panel_vc <- panel %>%
  select(country, year, npl, nim, credit_growth, log_vc_ai, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth), !is.na(log_vc_ai))

panel_vc_lag <- panel %>%
  select(country, year, npl, nim, credit_growth, lag_log_vc_ai, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth), !is.na(lag_log_vc_ai))

# Base dinâmica principal: inclui NPL(t-1) e IA contemporânea.
panel_dyn <- panel %>%
  select(country, year, npl, lag_npl, log_vc_ai, nim, credit_growth, covid) %>%
  filter(!is.na(npl), !is.na(lag_npl), !is.na(log_vc_ai),
         !is.na(nim), !is.na(credit_growth))

# Base dinâmica alternativa: inclui NPL(t-1) e IA defasada.
panel_dyn_lag_ai <- panel %>%
  select(country, year, npl, lag_npl, lag_log_vc_ai, nim, credit_growth, covid) %>%
  filter(!is.na(npl), !is.na(lag_npl), !is.na(lag_log_vc_ai),
         !is.na(nim), !is.na(credit_growth))

cat("Obs modelo base:", nrow(panel_base), "\n")
cat("Obs modelo VC_AI:", nrow(panel_vc), "\n")
cat("Obs modelo VC_AI defasado:", nrow(panel_vc_lag), "\n")
cat("Obs modelo dinâmico VC_AI:", nrow(panel_dyn), "\n")
cat("Obs modelo dinâmico VC_AI defasado:", nrow(panel_dyn_lag_ai), "\n")

if (nrow(panel_dyn) == 0) stop("Modelo dinâmico vazio.")

# ==============================================================================
# 8. MODELOS ESTÁTICOS DE EFEITOS FIXOS
#    Mantidos para comparação com a versão anterior.
# ==============================================================================
pdata_base <- pdata.frame(panel_base, index = c("country", "year"))
pdata_vc <- pdata.frame(panel_vc, index = c("country", "year"))
pdata_vc_lag <- pdata.frame(panel_vc_lag, index = c("country", "year"))

model_base_fe <- plm(
  npl ~ credit_growth + nim + covid,
  data = pdata_base,
  model = "within",
  effect = "individual"
)

model_vc_fe <- plm(
  npl ~ log_vc_ai + credit_growth + nim + covid,
  data = pdata_vc,
  model = "within",
  effect = "individual"
)

model_vc_lag_fe <- plm(
  npl ~ lag_log_vc_ai + credit_growth + nim + covid,
  data = pdata_vc_lag,
  model = "within",
  effect = "individual"
)

model_base_re <- plm(
  npl ~ credit_growth + nim + covid,
  data = pdata_base,
  model = "random",
  effect = "individual"
)

model_vc_re <- plm(
  npl ~ log_vc_ai + credit_growth + nim + covid,
  data = pdata_vc,
  model = "random",
  effect = "individual"
)

model_vc_lag_re <- plm(
  npl ~ lag_log_vc_ai + credit_growth + nim + covid,
  data = pdata_vc_lag,
  model = "random",
  effect = "individual"
)

cat("\n===== MODELOS FE ESTÁTICOS =====\n")
print(coeftest(model_base_fe, vcov = safe_vcov_fe(model_base_fe)))
print(coeftest(model_vc_fe, vcov = safe_vcov_fe(model_vc_fe)))
print(coeftest(model_vc_lag_fe, vcov = safe_vcov_fe(model_vc_lag_fe)))

# ==============================================================================
# 9. MODELO DINÂMICO POR EFEITOS FIXOS (WITHIN)
#    Atenção: com T pequeno, este estimador pode ter viés de Nickell.
# ==============================================================================
pdata_dyn <- pdata.frame(panel_dyn, index = c("country", "year"))
pdata_dyn_lag_ai <- pdata.frame(panel_dyn_lag_ai, index = c("country", "year"))

model_fe_dyn <- plm(
  npl ~ lag_npl + log_vc_ai + credit_growth + nim + covid,
  data = pdata_dyn,
  model = "within",
  effect = "individual"
)

model_fe_dyn_lag_ai <- plm(
  npl ~ lag_npl + lag_log_vc_ai + credit_growth + nim + covid,
  data = pdata_dyn_lag_ai,
  model = "within",
  effect = "individual"
)

model_re_dyn <- plm(
  npl ~ lag_npl + log_vc_ai + credit_growth + nim + covid,
  data = pdata_dyn,
  model = "random",
  effect = "individual"
)

model_re_dyn_lag_ai <- plm(
  npl ~ lag_npl + lag_log_vc_ai + credit_growth + nim + covid,
  data = pdata_dyn_lag_ai,
  model = "random",
  effect = "individual"
)

cat("\n===== FE DINÂMICO - IA CONTEMPORÂNEA =====\n")
print(summary(model_fe_dyn))
print(coeftest(model_fe_dyn, vcov = safe_vcov_fe(model_fe_dyn)))

cat("\n===== FE DINÂMICO - IA DEFASADA =====\n")
print(summary(model_fe_dyn_lag_ai))
print(coeftest(model_fe_dyn_lag_ai, vcov = safe_vcov_fe(model_fe_dyn_lag_ai)))

# ==============================================================================
# 10. TESTES DE DIAGNÓSTICO DOS MODELOS FE
# ==============================================================================
cat("\n===== TESTES DE DIAGNÓSTICO DOS MODELOS FE =====\n")

hausman_base <- rodar_teste(phtest(model_base_fe, model_base_re))
hausman_vc <- rodar_teste(phtest(model_vc_fe, model_vc_re))
hausman_vc_lag <- rodar_teste(phtest(model_vc_lag_fe, model_vc_lag_re))
hausman_dyn <- rodar_teste(phtest(model_fe_dyn, model_re_dyn))
hausman_dyn_lag_ai <- rodar_teste(phtest(model_fe_dyn_lag_ai, model_re_dyn_lag_ai))

wooldridge_base <- rodar_teste(pwartest(model_base_fe))
wooldridge_vc <- rodar_teste(pwartest(model_vc_fe))
wooldridge_vc_lag <- rodar_teste(pwartest(model_vc_lag_fe))
wooldridge_dyn <- rodar_teste(pwartest(model_fe_dyn))
wooldridge_dyn_lag_ai <- rodar_teste(pwartest(model_fe_dyn_lag_ai))

wald_base <- rodar_teste(bptest(model_base_fe, studentize = FALSE))
wald_vc <- rodar_teste(bptest(model_vc_fe, studentize = FALSE))
wald_vc_lag <- rodar_teste(bptest(model_vc_lag_fe, studentize = FALSE))
wald_dyn <- rodar_teste(bptest(model_fe_dyn, studentize = FALSE))
wald_dyn_lag_ai <- rodar_teste(bptest(model_fe_dyn_lag_ai, studentize = FALSE))

pesaran_base <- rodar_teste(pcdtest(model_base_fe, test = "cd"))
pesaran_vc <- rodar_teste(pcdtest(model_vc_fe, test = "cd"))
pesaran_vc_lag <- rodar_teste(pcdtest(model_vc_lag_fe, test = "cd"))
pesaran_dyn <- rodar_teste(pcdtest(model_fe_dyn, test = "cd"))
pesaran_dyn_lag_ai <- rodar_teste(pcdtest(model_fe_dyn_lag_ai, test = "cd"))

diag_fe <- bind_rows(
  extrair_teste(hausman_base, "FE estatico - base", "Hausman FE vs RE"),
  extrair_teste(hausman_vc, "FE estatico - VC_AI", "Hausman FE vs RE"),
  extrair_teste(hausman_vc_lag, "FE estatico - VC_AI defasado", "Hausman FE vs RE"),
  extrair_teste(hausman_dyn, "FE dinamico - VC_AI", "Hausman FE vs RE"),
  extrair_teste(hausman_dyn_lag_ai, "FE dinamico - VC_AI defasado", "Hausman FE vs RE"),
  extrair_teste(wooldridge_base, "FE estatico - base", "Wooldridge autocorrelacao"),
  extrair_teste(wooldridge_vc, "FE estatico - VC_AI", "Wooldridge autocorrelacao"),
  extrair_teste(wooldridge_vc_lag, "FE estatico - VC_AI defasado", "Wooldridge autocorrelacao"),
  extrair_teste(wooldridge_dyn, "FE dinamico - VC_AI", "Wooldridge autocorrelacao"),
  extrair_teste(wooldridge_dyn_lag_ai, "FE dinamico - VC_AI defasado", "Wooldridge autocorrelacao"),
  extrair_teste(wald_base, "FE estatico - base", "Wald/Breusch-Pagan heterocedasticidade"),
  extrair_teste(wald_vc, "FE estatico - VC_AI", "Wald/Breusch-Pagan heterocedasticidade"),
  extrair_teste(wald_vc_lag, "FE estatico - VC_AI defasado", "Wald/Breusch-Pagan heterocedasticidade"),
  extrair_teste(wald_dyn, "FE dinamico - VC_AI", "Wald/Breusch-Pagan heterocedasticidade"),
  extrair_teste(wald_dyn_lag_ai, "FE dinamico - VC_AI defasado", "Wald/Breusch-Pagan heterocedasticidade"),
  extrair_teste(pesaran_base, "FE estatico - base", "Pesaran CD dependencia seccional"),
  extrair_teste(pesaran_vc, "FE estatico - VC_AI", "Pesaran CD dependencia seccional"),
  extrair_teste(pesaran_vc_lag, "FE estatico - VC_AI defasado", "Pesaran CD dependencia seccional"),
  extrair_teste(pesaran_dyn, "FE dinamico - VC_AI", "Pesaran CD dependencia seccional"),
  extrair_teste(pesaran_dyn_lag_ai, "FE dinamico - VC_AI defasado", "Pesaran CD dependencia seccional")
)

print(diag_fe)

# ==============================================================================
# 11. DIFFERENCE GMM E SYSTEM GMM
#     Especificação principal: NPL(t) em função de NPL(t-1), log(VC_AI+1),
#     crescimento do crédito, NIM e dummy COVID.
## ==============================================================================

cat("\n===== GMM DINÂMICO - IA CONTEMPORÂNEA =====\n")

model_diff_gmm <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + log_vc_ai + credit_growth + nim + covid |
      lag(npl, 2:3) + lag(log_vc_ai, 2:3) + credit_growth + nim + covid,
    data = pdata_dyn,
    effect = "individual",
    model = "twosteps",
    transformation = "d",
    collapse = TRUE
  ),
  error = function(e) {
    cat("Erro no Difference GMM:", conditionMessage(e), "\n")
    NULL
  }
)

model_sys_gmm <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + log_vc_ai + credit_growth + nim + covid |
      lag(npl, 2:3) + lag(log_vc_ai, 2:3) + credit_growth + nim + covid,
    data = pdata_dyn,
    effect = "individual",
    model = "twosteps",
    transformation = "ld",
    collapse = TRUE
  ),
  error = function(e) {
    cat("Erro no System GMM:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(model_diff_gmm)) {
  cat("\n--- Difference GMM ---\n")
  print(summary(model_diff_gmm, robust = TRUE))
  print(mtest(model_diff_gmm, order = 1))
  print(mtest(model_diff_gmm, order = 2))
  print(sargan(model_diff_gmm))
}

if (!is.null(model_sys_gmm)) {
  cat("\n--- System GMM ---\n")
  print(summary(model_sys_gmm, robust = TRUE))
  print(mtest(model_sys_gmm, order = 1))
  print(mtest(model_sys_gmm, order = 2))
  print(sargan(model_sys_gmm))
}

# ==============================================================================
# 12. GMM ALTERNATIVO COM IA DEFASADA
#     Útil porque seus resultados anteriores indicavam efeito mais forte da IA
#     defasada. Esta tabela entra como robustez, não precisa ser a principal.
# ==============================================================================
cat("\n===== GMM DINÂMICO - IA DEFASADA =====\n")

model_diff_gmm_lag_ai <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + lag_log_vc_ai + credit_growth + nim + covid |
      lag(npl, 2:3) + lag(lag_log_vc_ai, 2:3) + credit_growth + nim + covid,
    data = pdata_dyn_lag_ai,
    effect = "individual",
    model = "twosteps",
    transformation = "d",
    collapse = TRUE
  ),
  error = function(e) {
    cat("Erro no Difference GMM com IA defasada:", conditionMessage(e), "\n")
    NULL
  }
)

model_sys_gmm_lag_ai <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + lag_log_vc_ai + credit_growth + nim + covid |
      lag(npl, 2:3) + lag(lag_log_vc_ai, 2:3) + credit_growth + nim + covid,
    data = pdata_dyn_lag_ai,
    effect = "individual",
    model = "twosteps",
    transformation = "ld",
    collapse = TRUE
  ),
  error = function(e) {
    cat("Erro no System GMM com IA defasada:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(model_diff_gmm_lag_ai)) {
  cat("\n--- Difference GMM - IA defasada ---\n")
  print(summary(model_diff_gmm_lag_ai, robust = TRUE))
  print(mtest(model_diff_gmm_lag_ai, order = 1))
  print(mtest(model_diff_gmm_lag_ai, order = 2))
  print(sargan(model_diff_gmm_lag_ai))
}

if (!is.null(model_sys_gmm_lag_ai)) {
  cat("\n--- System GMM - IA defasada ---\n")
  print(summary(model_sys_gmm_lag_ai, robust = TRUE))
  print(mtest(model_sys_gmm_lag_ai, order = 1))
  print(mtest(model_sys_gmm_lag_ai, order = 2))
  print(sargan(model_sys_gmm_lag_ai))
}

# ==============================================================================
# 12a. CRITÉRIOS DE INFORMAÇÃO AIC E BIC — MODELOS FE
#      Os modelos GMM (pgmm) não suportam logLik; o AIC/BIC é calculado
#      apenas para os modelos plm (within). Fórmula: log-verosimilhança
#      concentrada sob normalidade dos resíduos, penalizando k parâmetros
#      de inclinação + sigma^2 (os efeitos fixos não são penalizados por
#      convenção — são parâmetros de incómodo, não de interesse).
#      AIC = -2*logLik + 2*(k+1); BIC = -2*logLik + log(n)*(k+1).
# ==============================================================================

calcular_aic_bic <- function(model, nome) {
  if (is.null(model)) return(NULL)
  res  <- model$residuals
  n    <- length(res)
  k    <- length(coef(model))
  rss  <- sum(res^2)
  ll   <- -n/2 * (log(2*pi) + 1 - log(n) + log(rss))
  aic  <- -2*ll + 2*(k + 1)
  bic  <- -2*ll + log(n)*(k + 1)
  r2w  <- 1 - rss / sum((res - mean(res))^2 + rss)
  data.frame(
    modelo    = nome,
    n_obs     = n,
    k_params  = k,
    logLik    = round(ll, 2),
    AIC       = round(aic, 2),
    BIC       = round(bic, 2),
    row.names = NULL
  )
}

ic_base     <- calcular_aic_bic(model_base_fe,      "FE base (sem macro)")
ic_vc       <- calcular_aic_bic(model_vc_fe,        "FE + VC_AI (sem macro)")
ic_vc_lag   <- calcular_aic_bic(model_vc_lag_fe,    "FE + VC_AI lag (sem macro)")
ic_dyn      <- calcular_aic_bic(model_fe_dyn,       "FE dinâmico + VC_AI (sem macro)")
ic_dyn_lag  <- calcular_aic_bic(model_fe_dyn_lag_ai,"FE dinâmico + VC_AI lag (sem macro)")

tabela_ic <- bind_rows(ic_base, ic_vc, ic_vc_lag, ic_dyn, ic_dyn_lag)

cat("\n===== CRITÉRIOS DE INFORMAÇÃO (AIC / BIC) — FE SEM MACRO =====\n")
print(tabela_ic)
write.csv(tabela_ic, "aic_bic_fe_sem_macro.csv", row.names = FALSE)

# ==============================================================================
# 12b. TESTE DIFFERENCE-IN-HANSEN (manual — Roodman, 2009)
#      Avalia se os instrumentos adicionais de nível do System GMM são válidos
#      face ao Difference GMM. Estatística: Hansen(Diff) − Hansen(Sys),
#      distribuída chi-quadrado com gl = gl(Diff) − gl(Sys).
#      H0: instrumentos adicionais do Sys-GMM são exógenos.
#      Não rejeitar H0 → Sys-GMM preferido; rejeitar → Diff-GMM mais seguro.
# ==============================================================================

calcular_diff_hansen <- function(model_diff, model_sys, nome) {
  if (is.null(model_diff) || is.null(model_sys)) return(NULL)

  sargan_diff <- tryCatch(sargan(model_diff), error = function(e) NULL)
  sargan_sys  <- tryCatch(sargan(model_sys),  error = function(e) NULL)
  if (is.null(sargan_diff) || is.null(sargan_sys)) return(NULL)

  stat_diff <- unname(sargan_diff$statistic)
  stat_sys  <- unname(sargan_sys$statistic)
  df_diff   <- unname(sargan_diff$parameter)
  df_sys    <- unname(sargan_sys$parameter)

  stat_dh <- stat_diff - stat_sys
  df_dh   <- df_diff   - df_sys

  if (is.na(stat_dh) || df_dh <= 0) {
    cat("Aviso: graus de liberdade não positivos para", nome, "— teste não calculado.\n")
    return(NULL)
  }

  pval_dh <- 1 - pchisq(stat_dh, df = df_dh)

  data.frame(
    modelo        = nome,
    teste         = "Difference-in-Hansen (manual)",
    stat_diff_gmm = stat_diff,
    stat_sys_gmm  = stat_sys,
    estatistica   = stat_dh,
    df            = df_dh,
    p_valor       = pval_dh,
    resultado_5pct = ifelse(pval_dh < 0.05, "Rejeita H0 — usar Diff-GMM",
                                             "Não rejeita H0 — Sys-GMM válido"),
    row.names = NULL
  )
}

diff_hansen_contemp <- calcular_diff_hansen(
  model_diff_gmm, model_sys_gmm,
  "Diff-in-Hansen - VC_AI contemporaneo"
)

diff_hansen_lag <- calcular_diff_hansen(
  model_diff_gmm_lag_ai, model_sys_gmm_lag_ai,
  "Diff-in-Hansen - VC_AI defasado"
)

diff_hansen <- bind_rows(diff_hansen_contemp, diff_hansen_lag)

cat("\n===== DIFFERENCE-IN-HANSEN =====\n")
print(diff_hansen)

# ==============================================================================
# 13. EXPORTAÇÃO DOS RESULTADOS
# ==============================================================================
res_estaticos <- bind_rows(
  tidy_plm(model_base_fe, "FE estatico - base"),
  tidy_plm(model_vc_fe, "FE estatico - VC_AI"),
  tidy_plm(model_vc_lag_fe, "FE estatico - VC_AI defasado")
)

res_dinamicos <- bind_rows(
  tidy_plm(model_fe_dyn, "FE dinamico - VC_AI"),
  tidy_plm(model_fe_dyn_lag_ai, "FE dinamico - VC_AI defasado"),
  if (!is.null(model_diff_gmm)) tidy_gmm(model_diff_gmm, "Difference GMM - VC_AI") else NULL,
  if (!is.null(model_sys_gmm)) tidy_gmm(model_sys_gmm, "System GMM - VC_AI") else NULL,
  if (!is.null(model_diff_gmm_lag_ai)) tidy_gmm(model_diff_gmm_lag_ai, "Difference GMM - VC_AI defasado") else NULL,
  if (!is.null(model_sys_gmm_lag_ai)) tidy_gmm(model_sys_gmm_lag_ai, "System GMM - VC_AI defasado") else NULL
)

diag_gmm <- bind_rows(
  if (!is.null(model_diff_gmm)) extrair_diagnosticos_gmm(model_diff_gmm, "Difference GMM - VC_AI") else NULL,
  if (!is.null(model_sys_gmm)) extrair_diagnosticos_gmm(model_sys_gmm, "System GMM - VC_AI") else NULL,
  if (!is.null(model_diff_gmm_lag_ai)) extrair_diagnosticos_gmm(model_diff_gmm_lag_ai, "Difference GMM - VC_AI defasado") else NULL,
  if (!is.null(model_sys_gmm_lag_ai)) extrair_diagnosticos_gmm(model_sys_gmm_lag_ai, "System GMM - VC_AI defasado") else NULL
)

# Estatísticas descritivas completas, incluindo min e max.
desc <- panel %>%
  summarise(
    n_obs = n(),
    npl_media = mean(npl, na.rm = TRUE),
    npl_dp = sd(npl, na.rm = TRUE),
    npl_min = min(npl, na.rm = TRUE),
    npl_max = max(npl, na.rm = TRUE),
    nim_media = mean(nim, na.rm = TRUE),
    nim_dp = sd(nim, na.rm = TRUE),
    nim_min = min(nim, na.rm = TRUE),
    nim_max = max(nim, na.rm = TRUE),
    credit_growth_media = mean(credit_growth, na.rm = TRUE),
    credit_growth_dp = sd(credit_growth, na.rm = TRUE),
    credit_growth_min = min(credit_growth, na.rm = TRUE),
    credit_growth_max = max(credit_growth, na.rm = TRUE),
    vc_ai_media = mean(vc_ai, na.rm = TRUE),
    vc_ai_dp = sd(vc_ai, na.rm = TRUE),
    vc_ai_min = min(vc_ai, na.rm = TRUE),
    vc_ai_max = max(vc_ai, na.rm = TRUE),
    log_vc_ai_media = mean(log_vc_ai, na.rm = TRUE),
    log_vc_ai_dp = sd(log_vc_ai, na.rm = TRUE),
    log_vc_ai_min = min(log_vc_ai, na.rm = TRUE),
    log_vc_ai_max = max(log_vc_ai, na.rm = TRUE)
  )

correlacoes <- panel %>%
  select(npl, lag_npl, nim, credit_growth, vc_ai, log_vc_ai, lag_log_vc_ai) %>%
  cor(use = "pairwise.complete.obs")

# Tabela compacta: FE dinâmico, Diff-GMM e Sys-GMM.
parte1 <- res_dinamicos %>%
  filter(modelo %in% c("FE dinamico - VC_AI", "Difference GMM - VC_AI", "System GMM - VC_AI"))

write.csv(res_estaticos, "resultados_fe_estaticos.csv", row.names = FALSE)
write.csv(res_dinamicos, "resultados_dinamicos_fe_gmm.csv", row.names = FALSE)
write.csv(diag_fe, "diagnosticos_fe.csv", row.names = FALSE)
write.csv(diag_gmm, "diagnosticos_gmm.csv", row.names = FALSE)
write.csv(diff_hansen, "diff_in_hansen.csv", row.names = FALSE)
write.csv(tabela_ic, "aic_bic_fe_sem_macro.csv", row.names = FALSE)
write.csv(parte1, "tabela_fe_diffgmm_sysgmm.csv", row.names = FALSE)

wb <- createWorkbook()

addWorksheet(wb, "Painel")
writeData(wb, "Painel", panel)

addWorksheet(wb, "FE_Estaticos")
writeData(wb, "FE_Estaticos", res_estaticos)

addWorksheet(wb, "FE_GMM_Dinamicos")
writeData(wb, "FE_GMM_Dinamicos", res_dinamicos)

addWorksheet(wb, "Tabela_FE_GMM")
writeData(wb, "Tabela_FE_GMM", parte1)

addWorksheet(wb, "Diagnosticos_FE")
writeData(wb, "Diagnosticos_FE", diag_fe)

addWorksheet(wb, "Diagnosticos_GMM")
writeData(wb, "Diagnosticos_GMM", diag_gmm)

addWorksheet(wb, "Diff_in_Hansen")
writeData(wb, "Diff_in_Hansen", diff_hansen)

addWorksheet(wb, "Descritivas")
writeData(wb, "Descritivas", desc)

addWorksheet(wb, "Correlacoes")
writeData(wb, "Correlacoes", as.data.frame(correlacoes), rowNames = TRUE)

# ==============================================================================
# 14. GRÁFICOS
# ==============================================================================
cat("\nGerando gráficos...\n")

plot_dir <- file.path(OUT_DIR, "graficos")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

tema_grafico <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

salvar_grafico <- function(grafico, nome, largura = 9, altura = 5) {
  caminho <- file.path(plot_dir, nome)
  ggsave(filename = caminho, plot = grafico, width = largura, height = altura, dpi = 150)
  return(caminho)
}

g1 <- panel %>%
  group_by(year) %>%
  summarise(npl = mean(npl, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = npl)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "NPL médio - SSM12", x = "Ano", y = "NPL médio (%)") +
  tema_grafico

arq_g1 <- salvar_grafico(g1, "grafico_npl_medio.png")

g2 <- panel %>%
  group_by(year) %>%
  summarise(log_vc_ai = mean(log_vc_ai, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = log_vc_ai)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "VC_AI médio - SSM12", x = "Ano", y = "log(VC_AI + 1)") +
  tema_grafico

arq_g2 <- salvar_grafico(g2, "grafico_vc_ai_medio.png")

g3 <- panel %>%
  filter(!is.na(npl), !is.na(log_vc_ai)) %>%
  ggplot(aes(x = log_vc_ai, y = npl)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Relação entre VC_AI e NPL", x = "log(VC_AI + 1)", y = "NPL (%)") +
  tema_grafico

arq_g3 <- salvar_grafico(g3, "grafico_vc_ai_vs_npl.png")

g4 <- panel %>%
  filter(!is.na(log_vc_ai)) %>%
  ggplot(aes(x = year, y = log_vc_ai, group = country, colour = country)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "VC_AI por país - SSM12", x = "Ano", y = "log(VC_AI + 1)", colour = "País") +
  tema_grafico

arq_g4 <- salvar_grafico(g4, "grafico_vc_ai_por_pais.png", largura = 11, altura = 6)

g5 <- panel %>%
  filter(!is.na(npl)) %>%
  ggplot(aes(x = year, y = npl, group = country, colour = country)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "NPL por país - SSM12", x = "Ano", y = "NPL (%)", colour = "País") +
  tema_grafico

arq_g5 <- salvar_grafico(g5, "grafico_npl_por_pais.png", largura = 11, altura = 6)

g6 <- panel %>%
  filter(!is.na(credit_growth), !is.na(npl)) %>%
  ggplot(aes(x = credit_growth, y = npl)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Relação entre crescimento do crédito e NPL", x = "Crescimento do crédito (%)", y = "NPL (%)") +
  tema_grafico

arq_g6 <- salvar_grafico(g6, "grafico_credit_growth_vs_npl.png")

# ── GRÁFICOS DE DIAGNÓSTICO (novos) ──────────────────────────────────────────

# G7 — AIC/BIC por modelo FE (barras agrupadas)
g7 <- tabela_ic %>%
  tidyr::pivot_longer(cols = c(AIC, BIC), names_to = "criterio", values_to = "valor") %>%
  mutate(modelo = factor(modelo, levels = rev(tabela_ic$modelo))) %>%
  ggplot(aes(x = modelo, y = valor, fill = criterio)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("AIC" = "#2C5F8A", "BIC" = "#C05A28")) +
  labs(
    title    = "Critérios de Informação AIC e BIC — Modelos FE (sem macro)",
    subtitle = "Valores mais baixos indicam melhor ajustamento penalizado",
    x        = NULL, y = "Valor do critério",
    fill     = "Critério"
  ) +
  tema_grafico

arq_g7 <- salvar_grafico(g7, "grafico_aic_bic_fe.png", largura = 11, altura = 5)

# G8 — Coeficientes VC_AI nos modelos FE estáticos com IC 95%
coef_fe <- data.frame(
  modelo   = c("M2: FE VC_AI", "M3: FE VC_AI lag"),
  coef     = c(-1.848, -1.844),
  se       = c(0.483, 0.196),
  stringsAsFactors = FALSE
) %>%
  mutate(
    ic_low  = coef - 1.96 * se,
    ic_high = coef + 1.96 * se,
    modelo  = factor(modelo, levels = modelo)
  )

g8 <- ggplot(coef_fe, aes(x = modelo, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(
    aes(ymin = ic_low, ymax = ic_high),
    size = 0.8, colour = "#2C5F8A"
  ) +
  geom_text(aes(label = paste0(round(coef, 3), "***")),
            vjust = -1.2, size = 3.5, colour = "#2C5F8A") +
  labs(
    title    = "Coeficiente de log(VC_AI+1) — Modelos FE Estáticos",
    subtitle = "Intervalo de confiança 95% (erros Driscoll-Kraay)",
    x        = NULL, y = "Coeficiente estimado"
  ) +
  tema_grafico

arq_g8 <- salvar_grafico(g8, "grafico_coef_vcai_fe.png", largura = 8, altura = 5)

# G9 — Persistência do NPL: FE dinâmico vs Sys-GMM vs Diff-GMM
persistencia <- data.frame(
  modelo = c("FE Dinâmico (M5)", "Sys-GMM (M4)", "Diff-GMM (M4-D)"),
  coef   = c(0.946, 0.549, 0.555),
  tipo   = c("Enviesado (Nickell)", "Preferido", "Validação cruzada")
)

g9 <- ggplot(persistencia, aes(x = reorder(modelo, coef), y = coef, fill = tipo)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = c(0.4, 0.7), linetype = "dashed", colour = "grey40") +
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = 0.4, ymax = 0.7,
           alpha = 0.08, fill = "#1D9E75") +
  annotate("text", x = 0.6, y = 0.55, label = "Intervalo literatura\n(0,4–0,7)",
           size = 3, colour = "#1D9E75", hjust = 0) +
  geom_text(aes(label = round(coef, 3)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c(
    "Enviesado (Nickell)" = "#C05A28",
    "Preferido"          = "#2C5F8A",
    "Validação cruzada"  = "#888780"
  )) +
  labs(
    title    = "Persistência do NPL (coef. NPLt-1) por Estimador",
    subtitle = "Zona verde = intervalo reportado pela literatura (Beck et al., 2015; Ghosh, 2015)",
    x        = NULL, y = "Coeficiente de persistência", fill = ""
  ) +
  tema_grafico

arq_g9 <- salvar_grafico(g9, "grafico_persistencia_npl.png", largura = 9, altura = 5)

# G10 — Diff-in-Hansen: visualização do resultado
dh_dados <- data.frame(
  especificacao = c("VC_AI contemporâneo", "VC_AI defasado"),
  stat          = c(6.108, 6.798),
  pvalor        = c(0.296, 0.236),
  validado      = c("Não rejeita H0\n(Sys-GMM válido)", "Não rejeita H0\n(Sys-GMM válido)")
)

g10 <- ggplot(dh_dados, aes(x = especificacao, y = pvalor, fill = validado)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#C05A28", linewidth = 0.8) +
  annotate("text", x = 1.5, y = 0.06,
           label = "Limiar de rejeição H0 (p = 0,05)",
           size = 3.2, colour = "#C05A28") +
  geom_text(aes(label = paste0("p = ", round(pvalor, 3), "\nχ²(5) = ", round(stat, 3))),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Não rejeita H0\n(Sys-GMM válido)" = "#1D9E75")) +
  scale_y_continuous(limits = c(0, 0.45)) +
  labs(
    title    = "Teste Difference-in-Hansen — Validade dos Instrumentos do Sys-GMM",
    subtitle = "H0: instrumentos adicionais de nível são exógenos. Não rejeitar → Sys-GMM preferido",
    x        = NULL, y = "p-valor", fill = ""
  ) +
  tema_grafico

arq_g10 <- salvar_grafico(g10, "grafico_diff_in_hansen.png", largura = 9, altura = 5)

graficos_gerados <- data.frame(
  grafico = c(
    "NPL médio - SSM12",
    "VC_AI médio - SSM12",
    "VC_AI vs NPL",
    "VC_AI por país - SSM12",
    "NPL por país - SSM12",
    "Credit growth vs NPL",
    "AIC/BIC por modelo FE",
    "Coeficiente VC_AI nos FE",
    "Persistência NPL por estimador",
    "Diff-in-Hansen"
  ),
  arquivo = basename(c(arq_g1, arq_g2, arq_g3, arq_g4, arq_g5, arq_g6,
                       arq_g7, arq_g8, arq_g9, arq_g10)),
  caminho = c(arq_g1, arq_g2, arq_g3, arq_g4, arq_g5, arq_g6,
              arq_g7, arq_g8, arq_g9, arq_g10)
)

# Atualiza o Excel com uma aba de controle dos gráficos.
addWorksheet(wb, "Graficos")
writeData(wb, "Graficos", graficos_gerados)

# Insere imagens no Excel. Se der problema, o Excel é salvo sem travar o script.
inserir_imagem_excel <- function(wb, sheet, arquivo, linha, coluna, largura = 7, altura = 4) {
  tryCatch(
    insertImage(
      wb, sheet, arquivo,
      startRow = linha, startCol = coluna,
      width = largura, height = altura
    ),
    error = function(e) {
      cat("Aviso: não foi possível inserir imagem no Excel:", basename(arquivo), "\n")
    }
  )
}

inserir_imagem_excel(wb, "Graficos", arq_g1, 9,  1)
inserir_imagem_excel(wb, "Graficos", arq_g2, 9,  9)
inserir_imagem_excel(wb, "Graficos", arq_g3, 29, 1)
inserir_imagem_excel(wb, "Graficos", arq_g4, 29, 9)
inserir_imagem_excel(wb, "Graficos", arq_g5, 49, 1)
inserir_imagem_excel(wb, "Graficos", arq_g6, 49, 9)
inserir_imagem_excel(wb, "Graficos", arq_g7, 69, 1)
inserir_imagem_excel(wb, "Graficos", arq_g8, 69, 9)
inserir_imagem_excel(wb, "Graficos", arq_g9, 89, 1)
inserir_imagem_excel(wb, "Graficos", arq_g10, 89, 9)

arq_excel <- file.path(OUT_DIR, "resultados_completos_vc_ai_FE_GMM.xlsx")

salvar_excel <- function(wb, arquivo) {
  dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TMPDIR = TMP_DIR, TEMP = TMP_DIR, TMP = TMP_DIR)
  options(openxlsx.tempdir = TMP_DIR)

  tryCatch(
    saveWorkbook(wb, arquivo, overwrite = TRUE),
    error = function(e) {
      cat("Erro ao salvar Excel com imagens. Tentando salvar uma versão sem imagens...\n")

      wb2 <- createWorkbook()
      addWorksheet(wb2, "Painel")
      writeData(wb2, "Painel", panel)
      addWorksheet(wb2, "FE_Estaticos")
      writeData(wb2, "FE_Estaticos", res_estaticos)
      addWorksheet(wb2, "FE_GMM_Dinamicos")
      writeData(wb2, "FE_GMM_Dinamicos", res_dinamicos)
      addWorksheet(wb2, "Tabela_FE_GMM")
      writeData(wb2, "Tabela_FE_GMM", parte1)
      addWorksheet(wb2, "Diagnosticos_FE")
      writeData(wb2, "Diagnosticos_FE", diag_fe)
      addWorksheet(wb2, "Diagnosticos_GMM")
      writeData(wb2, "Diagnosticos_GMM", diag_gmm)
      addWorksheet(wb2, "Descritivas")
      writeData(wb2, "Descritivas", desc)
      addWorksheet(wb2, "Correlacoes")
      writeData(wb2, "Correlacoes", as.data.frame(correlacoes), rowNames = TRUE)
      addWorksheet(wb2, "Graficos")
      writeData(wb2, "Graficos", graficos_gerados)

      saveWorkbook(wb2, arquivo, overwrite = TRUE)
    }
  )
}

salvar_excel(wb, arq_excel)

# ==============================================================================
# 15. FINAL
# ==============================================================================
cat("\nSCRIPT CONCLUÍDO\n")
cat("Arquivos gerados em:", OUT_DIR, "\n")
cat("- painel_final_vc_ai.csv\n")
cat("- resultados_fe_estaticos.csv\n")
cat("- resultados_dinamicos_fe_gmm.csv\n")
cat("- diagnosticos_fe.csv\n")
cat("- diagnosticos_gmm.csv\n")
cat("- diff_in_hansen.csv\n")
cat("- tabela_fe_diffgmm_sysgmm.csv\n")
cat("- resultados_completos_vc_ai_FE_GMM.xlsx\n")
cat("- pasta graficos com arquivos .png\n")

################################################################################
