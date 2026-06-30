################################################################################
# DISSERTAÇÃO - PAINEL SSM12: VC_AI, NPL, FE DINÂMICO, DIFF-GMM E SYS-GMM
# VERSÃO COM CONTROLOS MACROECONÓMICOS (PIB, INFLAÇÃO, DESEMPREGO)
# ################################################################################

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

# ------------------------------------------------------------------------------
# Leitor de ficheiros Eurostat exportados em formato "wide" 
# ------------------------------------------------------------------------------
read_eurostat_wide <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)

  # localizar a linha "TIME"
  time_row_idx <- which(raw[[1]] == "TIME")[1]
  if (is.na(time_row_idx)) stop(paste("Linha TIME não encontrada em", path, sheet))

  header <- raw[time_row_idx, ]
  geo_row_idx <- time_row_idx + 1
  data_start_idx <- geo_row_idx + 1

  # identificar colunas de ano (valores numéricos de 4 dígitos no cabeçalho)
  year_cols <- list()
  for (j in 2:ncol(header)) {
    v <- suppressWarnings(as.numeric(as.character(header[[j]])))
    if (!is.na(v) && v >= 2000 && v <= 2030) {
      year_cols[[as.character(v)]] <- j
    }
  }

  dados_pais <- raw[data_start_idx:nrow(raw), ]

  out <- list()
  for (i in seq_len(nrow(dados_pais))) {
    pais_raw <- as.character(dados_pais[[i, 1]])
    pais <- normalizar_pais(pais_raw)
    if (is.na(pais)) next
    for (ano_str in names(year_cols)) {
      col_idx <- year_cols[[ano_str]]
      valor <- suppressWarnings(as.numeric(as.character(dados_pais[[i, col_idx]])))
      out[[length(out) + 1]] <- data.frame(
        country = pais,
        year = as.integer(ano_str),
        value = valor
      )
    }
  }

  bind_rows(out) %>%
    filter(!is.na(value), year >= 2015, year <= 2023) %>%
    group_by(country, year) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
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
# 6. VARIÁVEIS MACROECONÓMICAS (EUROSTAT) — NOVO
#==============================================================================
cat("Carregando variáveis macroeconómicas (Eurostat)...\n")

# PIB per capita em níveis (euros correntes);
gdp_levels <- read_eurostat_wide(
  file.path(DATA_DIR, "GDP_indice_por_pais_Eurostat.xlsx"),
  "Sheet 11"
) %>%
  rename(gdp_pc = value)

gdp_growth_df <- gdp_levels %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(gdp_growth = 100 * (gdp_pc / lag(gdp_pc) - 1)) %>%
  ungroup() %>%
  select(country, year, gdp_growth)

# Inflação: HICP, "Annual average rate of change"
inflation_df <- read_eurostat_wide(
  file.path(DATA_DIR, "HICP_indice_por_pais_Eurostat.xlsx"),
  "Sheet 6"
) %>%
  rename(inflation = value)

# Desemprego: 

unemployment_df <- read_eurostat_wide(
  file.path(DATA_DIR, "Desemprego_taxa_por_pais_Eurostat.xlsx"),
  "Sheet 3"
) %>%
  rename(unemployment = value)

cat("Obs GDP growth:", nrow(gdp_growth_df), "\n")
cat("Obs Inflação:", nrow(inflation_df), "\n")
cat("Obs Desemprego:", nrow(unemployment_df), "\n")

# ==============================================================================
# 7. PAINEL FINAL
# ==============================================================================
panel <- expand.grid(country = SSM12, year = 2015:2023,
                     stringsAsFactors = FALSE) %>%
  left_join(npl, by = c("country", "year")) %>%
  left_join(nim, by = c("country", "year")) %>%
  left_join(loans, by = c("country", "year")) %>%
  left_join(vc_ai, by = c("country", "year")) %>%
  left_join(gdp_growth_df, by = c("country", "year")) %>%
  left_join(inflation_df, by = c("country", "year")) %>%
  left_join(unemployment_df, by = c("country", "year")) %>%
  mutate(covid = as.integer(year %in% c(2020, 2021))) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    lag_npl = lag(npl),
    lag_nim = lag(nim),
    lag_credit = lag(credit_growth),
    lag_gdp_growth = lag(gdp_growth),
    lag_inflation = lag(inflation),
    lag_unemployment = lag(unemployment)
  ) %>%
  ungroup()

cat("Painel bruto:", nrow(panel), "\n")
print(colSums(is.na(panel)))

write.csv(panel, "painel_final_vc_ai_macro.csv", row.names = FALSE)

# ==============================================================================
# 8. BASES DOS MODELOS
# ==============================================================================
panel_base <- panel %>%
  select(country, year, npl, nim, credit_growth, gdp_growth, inflation,
         unemployment, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth),
         !is.na(gdp_growth), !is.na(inflation), !is.na(unemployment))

panel_vc <- panel %>%
  select(country, year, npl, nim, credit_growth, gdp_growth, inflation,
         unemployment, log_vc_ai, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth), !is.na(log_vc_ai),
         !is.na(gdp_growth), !is.na(inflation), !is.na(unemployment))

panel_vc_lag <- panel %>%
  select(country, year, npl, nim, credit_growth, gdp_growth, inflation,
         unemployment, lag_log_vc_ai, covid) %>%
  filter(!is.na(npl), !is.na(nim), !is.na(credit_growth), !is.na(lag_log_vc_ai),
         !is.na(gdp_growth), !is.na(inflation), !is.na(unemployment))

# Base dinâmica principal: inclui NPL(t-1) e IA contemporânea.
panel_dyn <- panel %>%
  select(country, year, npl, lag_npl, log_vc_ai, nim, credit_growth,
         gdp_growth, inflation, unemployment, covid) %>%
  filter(!is.na(npl), !is.na(lag_npl), !is.na(log_vc_ai),
         !is.na(nim), !is.na(credit_growth),
         !is.na(gdp_growth), !is.na(inflation), !is.na(unemployment))

# Base dinâmica alternativa: inclui NPL(t-1) e IA defasada.
panel_dyn_lag_ai <- panel %>%
  select(country, year, npl, lag_npl, lag_log_vc_ai, nim, credit_growth,
         gdp_growth, inflation, unemployment, covid) %>%
  filter(!is.na(npl), !is.na(lag_npl), !is.na(lag_log_vc_ai),
         !is.na(nim), !is.na(credit_growth),
         !is.na(gdp_growth), !is.na(inflation), !is.na(unemployment))

cat("Obs modelo base:", nrow(panel_base), "\n")
cat("Obs modelo VC_AI:", nrow(panel_vc), "\n")
cat("Obs modelo VC_AI defasado:", nrow(panel_vc_lag), "\n")
cat("Obs modelo dinâmico VC_AI:", nrow(panel_dyn), "\n")
cat("Obs modelo dinâmico VC_AI defasado:", nrow(panel_dyn_lag_ai), "\n")

if (nrow(panel_dyn) == 0) stop("Modelo dinâmico vazio.")

# ==============================================================================
# 9. MODELOS ESTÁTICOS DE EFEITOS FIXOS
# ==============================================================================
pdata_base <- pdata.frame(panel_base, index = c("country", "year"))
pdata_vc <- pdata.frame(panel_vc, index = c("country", "year"))
pdata_vc_lag <- pdata.frame(panel_vc_lag, index = c("country", "year"))

model_base_fe <- plm(
  npl ~ credit_growth + nim + gdp_growth + inflation + unemployment + covid,
  data = pdata_base,
  model = "within",
  effect = "individual"
)

model_vc_fe <- plm(
  npl ~ log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_vc,
  model = "within",
  effect = "individual"
)

model_vc_lag_fe <- plm(
  npl ~ lag_log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_vc_lag,
  model = "within",
  effect = "individual"
)

model_base_re <- plm(
  npl ~ credit_growth + nim + gdp_growth + inflation + unemployment + covid,
  data = pdata_base,
  model = "random",
  effect = "individual"
)

model_vc_re <- plm(
  npl ~ log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_vc,
  model = "random",
  effect = "individual"
)

model_vc_lag_re <- plm(
  npl ~ lag_log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_vc_lag,
  model = "random",
  effect = "individual"
)

cat("\n===== MODELOS FE ESTÁTICOS =====\n")
print(coeftest(model_base_fe, vcov = safe_vcov_fe(model_base_fe)))
print(coeftest(model_vc_fe, vcov = safe_vcov_fe(model_vc_fe)))
print(coeftest(model_vc_lag_fe, vcov = safe_vcov_fe(model_vc_lag_fe)))

# ==============================================================================
# 10. MODELO DINÂMICO POR EFEITOS FIXOS (WITHIN)
# ==============================================================================
pdata_dyn <- pdata.frame(panel_dyn, index = c("country", "year"))
pdata_dyn_lag_ai <- pdata.frame(panel_dyn_lag_ai, index = c("country", "year"))

model_fe_dyn <- plm(
  npl ~ lag_npl + log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_dyn,
  model = "within",
  effect = "individual"
)

model_fe_dyn_lag_ai <- plm(
  npl ~ lag_npl + lag_log_vc_ai + credit_growth + nim + gdp_growth +
    inflation + unemployment + covid,
  data = pdata_dyn_lag_ai,
  model = "within",
  effect = "individual"
)

model_re_dyn <- plm(
  npl ~ lag_npl + log_vc_ai + credit_growth + nim + gdp_growth + inflation +
    unemployment + covid,
  data = pdata_dyn,
  model = "random",
  effect = "individual"
)

model_re_dyn_lag_ai <- plm(
  npl ~ lag_npl + lag_log_vc_ai + credit_growth + nim + gdp_growth +
    inflation + unemployment + covid,
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
# 11. TESTES DE DIAGNÓSTICO DOS MODELOS FE
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
# 12. DIFFERENCE GMM E SYSTEM GMM
# ==============================================================================

cat("\n===== GMM DINÂMICO - IA CONTEMPORÂNEA =====\n")

model_diff_gmm <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + log_vc_ai + credit_growth + nim + gdp_growth +
      inflation + unemployment + covid |
      lag(npl, 2:3) + lag(log_vc_ai, 2:3) + credit_growth + nim +
      gdp_growth + inflation + unemployment + covid,
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
    npl ~ lag(npl, 1) + log_vc_ai + credit_growth + nim + gdp_growth +
      inflation + unemployment + covid |
      lag(npl, 2:3) + lag(log_vc_ai, 2:3) + credit_growth + nim +
      gdp_growth + inflation + unemployment + covid,
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
# 13. GMM ALTERNATIVO COM IA DEFASADA (robustez)
# ==============================================================================
cat("\n===== GMM DINÂMICO - IA DEFASADA =====\n")

model_diff_gmm_lag_ai <- tryCatch(
  pgmm(
    npl ~ lag(npl, 1) + lag_log_vc_ai + credit_growth + nim + gdp_growth +
      inflation + unemployment + covid |
      lag(npl, 2:3) + lag(lag_log_vc_ai, 2:3) + credit_growth + nim +
      gdp_growth + inflation + unemployment + covid,
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
    npl ~ lag(npl, 1) + lag_log_vc_ai + credit_growth + nim + gdp_growth +
      inflation + unemployment + covid |
      lag(npl, 2:3) + lag(lag_log_vc_ai, 2:3) + credit_growth + nim +
      gdp_growth + inflation + unemployment + covid,
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
# ==============================================================================
# 13a. CRITÉRIOS DE INFORMAÇÃO AIC E BIC — MODELOS FE COM MACRO
#      Comparação directa com os modelos sem variáveis macroeconómicas.
#      AIC/BIC menores indicam melhor ajustamento penalizado pelo n.º parâmetros.
#      Nota: amostras ligeiramente diferentes (macro exige gdp_growth não NA),
#      pelo que a comparação entre sem/com macro deve ser feita com cautela;
#      o BIC, mais penalizador, é o critério preferido quando N é pequeno.
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

ic_base_m    <- calcular_aic_bic(model_base_fe,       "FE base (com macro)")
ic_vc_m      <- calcular_aic_bic(model_vc_fe,         "FE + VC_AI (com macro)")
ic_vc_lag_m  <- calcular_aic_bic(model_vc_lag_fe,     "FE + VC_AI lag (com macro)")
ic_dyn_m     <- calcular_aic_bic(model_fe_dyn,        "FE dinâmico + VC_AI (com macro)")
ic_dyn_lag_m <- calcular_aic_bic(model_fe_dyn_lag_ai, "FE dinâmico + VC_AI lag (com macro)")

tabela_ic_macro <- bind_rows(ic_base_m, ic_vc_m, ic_vc_lag_m, ic_dyn_m, ic_dyn_lag_m)

cat("\n===== CRITÉRIOS DE INFORMAÇÃO (AIC / BIC) — FE COM MACRO =====\n")
print(tabela_ic_macro)
write.csv(tabela_ic_macro, "aic_bic_fe_com_macro.csv", row.names = FALSE)

# ==============================================================================
# 13b. TESTE DIFFERENCE-IN-HANSEN (manual — Roodman, 2009)
#      Testa se os instrumentos adicionais de nível do Sys-GMM são exógenos.
#      H0: instrumentos adicionais válidos → Sys-GMM preferido.
# ==============================================================================

calcular_diff_hansen <- function(model_diff, model_sys, nome) {
  if (is.null(model_diff) || is.null(model_sys)) return(NULL)

  sargan_diff <- tryCatch(sargan(model_diff), error = function(e) NULL)
  sargan_sys <- tryCatch(sargan(model_sys), error = function(e) NULL)
  if (is.null(sargan_diff) || is.null(sargan_sys)) return(NULL)

  stat_diff <- unname(sargan_diff$statistic)
  stat_sys <- unname(sargan_sys$statistic)
  df_diff <- unname(sargan_diff$parameter)
  df_sys <- unname(sargan_sys$parameter)

  stat_dh <- stat_diff - stat_sys
  df_dh <- df_diff - df_sys

  if (is.na(stat_dh) || df_dh <= 0) return(NULL)

  pval_dh <- 1 - pchisq(stat_dh, df = df_dh)

  data.frame(
    modelo = nome,
    teste = "Difference-in-Hansen (manual)",
    estatistica = stat_dh,
    p_valor = pval_dh,
    df = df_dh,
    resultado_5pct = ifelse(pval_dh < 0.05, "Rejeita H0", "Não rejeita H0"),
    row.names = NULL
  )
}

diff_hansen_contemp <- calcular_diff_hansen(
  model_diff_gmm, model_sys_gmm, "Diff-in-Hansen - VC_AI"
)
diff_hansen_lag <- calcular_diff_hansen(
  model_diff_gmm_lag_ai, model_sys_gmm_lag_ai, "Diff-in-Hansen - VC_AI defasado"
)

diff_hansen <- bind_rows(diff_hansen_contemp, diff_hansen_lag)
print(diff_hansen)

# ==============================================================================
# 14. EXPORTAÇÃO DOS RESULTADOS
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
    gdp_growth_media = mean(gdp_growth, na.rm = TRUE),
    gdp_growth_dp = sd(gdp_growth, na.rm = TRUE),
    gdp_growth_min = min(gdp_growth, na.rm = TRUE),
    gdp_growth_max = max(gdp_growth, na.rm = TRUE),
    inflation_media = mean(inflation, na.rm = TRUE),
    inflation_dp = sd(inflation, na.rm = TRUE),
    inflation_min = min(inflation, na.rm = TRUE),
    inflation_max = max(inflation, na.rm = TRUE),
    unemployment_media = mean(unemployment, na.rm = TRUE),
    unemployment_dp = sd(unemployment, na.rm = TRUE),
    unemployment_min = min(unemployment, na.rm = TRUE),
    unemployment_max = max(unemployment, na.rm = TRUE),
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
  select(npl, lag_npl, nim, credit_growth, gdp_growth, inflation,
         unemployment, vc_ai, log_vc_ai, lag_log_vc_ai) %>%
  cor(use = "pairwise.complete.obs")

# Tabela compacta: FE dinâmico, Diff-GMM e Sys-GMM.
parte1 <- res_dinamicos %>%
  filter(modelo %in% c("FE dinamico - VC_AI", "Difference GMM - VC_AI", "System GMM - VC_AI"))

write.csv(panel, "painel_final_vc_ai_macro.csv", row.names = FALSE)
write.csv(res_estaticos, "resultados_fe_estaticos_macro.csv", row.names = FALSE)
write.csv(res_dinamicos, "resultados_dinamicos_fe_gmm_macro.csv", row.names = FALSE)
write.csv(diag_fe, "diagnosticos_fe_macro.csv", row.names = FALSE)
write.csv(diag_gmm, "diagnosticos_gmm_macro.csv", row.names = FALSE)
write.csv(diff_hansen, "diff_in_hansen_macro.csv", row.names = FALSE)
write.csv(tabela_ic_macro, "aic_bic_fe_com_macro.csv", row.names = FALSE)
write.csv(parte1, "tabela_fe_diffgmm_sysgmm_macro.csv", row.names = FALSE)

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
# 15. GRÁFICOS
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

# NOVO: gráfico de dispersão crescimento do PIB vs NPL
g7 <- panel %>%
  filter(!is.na(gdp_growth), !is.na(npl)) %>%
  ggplot(aes(x = gdp_growth, y = npl)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Relação entre crescimento do PIB e NPL", x = "Crescimento do PIB (%, nominal)", y = "NPL (%)") +
  tema_grafico

arq_g7 <- salvar_grafico(g7, "grafico_gdp_growth_vs_npl.png")

# NOVO: gráfico de dispersão desemprego vs NPL
g8 <- panel %>%
  filter(!is.na(unemployment), !is.na(npl)) %>%
  ggplot(aes(x = unemployment, y = npl)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Relação entre desemprego e NPL", x = "Taxa de desemprego (%)", y = "NPL (%)") +
  tema_grafico

arq_g8 <- salvar_grafico(g8, "grafico_unemployment_vs_npl.png")

# ── GRÁFICOS DE DIAGNÓSTICO (novos) ──────────────────────────────────────────

# G9 — AIC/BIC por modelo FE com macro (barras agrupadas)
g9 <- tabela_ic_macro %>%
  tidyr::pivot_longer(cols = c(AIC, BIC), names_to = "criterio", values_to = "valor") %>%
  mutate(modelo = factor(modelo, levels = rev(tabela_ic_macro$modelo))) %>%
  ggplot(aes(x = modelo, y = valor, fill = criterio)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("AIC" = "#2C5F8A", "BIC" = "#C05A28")) +
  labs(
    title    = "Critérios de Informação AIC e BIC — Modelos FE (com macro)",
    subtitle = "Valores mais baixos indicam melhor ajustamento penalizado",
    x        = NULL, y = "Valor do critério", fill = "Critério"
  ) +
  tema_grafico

arq_g9 <- salvar_grafico(g9, "grafico_aic_bic_fe_macro.png", largura = 11, altura = 5)

# G10 — Comparação coeficientes VC_AI: sem macro vs com macro
coef_comp <- data.frame(
  modelo     = c("FE Estático\n(sem macro)", "FE Estático\n(com macro)",
                 "Sys-GMM\n(sem macro)", "Sys-GMM\n(com macro)"),
  coef       = c(-1.848, 0.658, -0.058, -1.540),
  se         = c(0.483, 0.301, 0.489, 0.829),
  versao     = c("Sem macro", "Com macro", "Sem macro", "Com macro"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    ic_low  = coef - 1.96 * se,
    ic_high = coef + 1.96 * se,
    modelo  = factor(modelo, levels = unique(modelo))
  )

g10 <- ggplot(coef_comp, aes(x = modelo, y = coef, colour = versao)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(
    aes(ymin = ic_low, ymax = ic_high),
    size = 0.8, position = position_dodge(width = 0)
  ) +
  geom_text(aes(label = round(coef, 3)),
            vjust = -1.2, size = 3.2) +
  scale_colour_manual(values = c("Sem macro" = "#2C5F8A", "Com macro" = "#C05A28")) +
  labs(
    title    = "Coeficiente de log(VC_AI+1): Sem vs Com Variáveis Macroeconómicas",
    subtitle = "IC 95%. Reversão de sinal nos FE estáticos evidencia endogeneidade por variáveis de confusão",
    x        = NULL, y = "Coeficiente estimado", colour = "Especificação"
  ) +
  tema_grafico

arq_g10 <- salvar_grafico(g10, "grafico_coef_vcai_comparacao.png", largura = 11, altura = 5)

# G11 — Persistência do NPL: comparação entre estimadores e versões
persistencia_comp <- data.frame(
  modelo = c("FE Din. (sem macro)", "Sys-GMM (sem macro)", "Diff-GMM (sem macro)",
             "FE Din. (com macro)", "Sys-GMM (com macro)", "Diff-GMM (com macro)"),
  coef   = c(0.946, 0.549, 0.555, 0.827, 0.394, 0.332),
  versao = c("Sem macro", "Sem macro", "Sem macro",
             "Com macro", "Com macro", "Com macro")
)

g11 <- ggplot(persistencia_comp, aes(x = reorder(modelo, coef), y = coef, fill = versao)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = c(0.4, 0.7), linetype = "dashed", colour = "grey40") +
  annotate("rect", xmin = 0.4, xmax = 6.6, ymin = 0.4, ymax = 0.7,
           alpha = 0.07, fill = "#1D9E75") +
  annotate("text", x = 0.65, y = 0.55, label = "Intervalo literatura (0,4–0,7)",
           size = 3, colour = "#1D9E75", hjust = 0) +
  geom_text(aes(label = round(coef, 3)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("Sem macro" = "#2C5F8A", "Com macro" = "#C05A28")) +
  coord_flip() +
  labs(
    title    = "Persistência do NPL (coef. NPLt-1) — Todos os Estimadores",
    subtitle = "Zona verde = intervalo da literatura. Sys-GMM sem macro é o único dentro do intervalo esperado",
    x        = NULL, y = "Coeficiente de persistência", fill = "Versão"
  ) +
  tema_grafico

arq_g11 <- salvar_grafico(g11, "grafico_persistencia_comp.png", largura = 11, altura = 5)

# G12 — Diff-in-Hansen (apenas modelo sem macro — único calculável)
dh_dados <- data.frame(
  especificacao = c("VC_AI contemporâneo", "VC_AI defasado"),
  stat          = c(6.108, 6.798),
  pvalor        = c(0.296, 0.236),
  resultado     = c("Válido", "Válido")
)

g12 <- ggplot(dh_dados, aes(x = especificacao, y = pvalor, fill = resultado)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#C05A28", linewidth = 0.8) +
  annotate("text", x = 1.5, y = 0.065,
           label = "Limiar de rejeição H0 (p = 0,05)",
           size = 3.2, colour = "#C05A28") +
  geom_text(aes(label = paste0("p = ", round(pvalor, 3), "\nχ²(5) = ", round(stat, 3))),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Válido" = "#1D9E75")) +
  scale_y_continuous(limits = c(0, 0.45)) +
  labs(
    title    = "Teste Difference-in-Hansen (modelo sem variáveis macro)",
    subtitle = "H0: instrumentos adicionais do Sys-GMM são exógenos. p > 0,05 → Sys-GMM válido",
    x        = NULL, y = "p-valor", fill = ""
  ) +
  tema_grafico

arq_g12 <- salvar_grafico(g12, "grafico_diff_in_hansen.png", largura = 9, altura = 5)

graficos_gerados <- data.frame(
  grafico = c(
    "NPL médio - SSM12",
    "VC_AI médio - SSM12",
    "VC_AI vs NPL",
    "VC_AI por país - SSM12",
    "NPL por país - SSM12",
    "Credit growth vs NPL",
    "GDP growth vs NPL",
    "Unemployment vs NPL",
    "AIC/BIC por modelo FE (com macro)",
    "Coeficiente VC_AI: sem vs com macro",
    "Persistência NPL — todos os estimadores",
    "Diff-in-Hansen (sem macro)"
  ),
  arquivo = basename(c(arq_g1, arq_g2, arq_g3, arq_g4, arq_g5, arq_g6,
                       arq_g7, arq_g8, arq_g9, arq_g10, arq_g11, arq_g12)),
  caminho = c(arq_g1, arq_g2, arq_g3, arq_g4, arq_g5, arq_g6,
              arq_g7, arq_g8, arq_g9, arq_g10, arq_g11, arq_g12)
)

addWorksheet(wb, "Graficos")
writeData(wb, "Graficos", graficos_gerados)

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
inserir_imagem_excel(wb, "Graficos", arq_g9,  89, 1)
inserir_imagem_excel(wb, "Graficos", arq_g10, 89, 9)
inserir_imagem_excel(wb, "Graficos", arq_g11, 109, 1)
inserir_imagem_excel(wb, "Graficos", arq_g12, 109, 9)

arq_excel <- file.path(OUT_DIR, "resultados_completos_vc_ai_FE_GMM_macro.xlsx")

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
      addWorksheet(wb2, "Diff_in_Hansen")
      writeData(wb2, "Diff_in_Hansen", diff_hansen)
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
# 16. FINAL
# ==============================================================================
cat("\nSCRIPT CONCLUÍDO\n")
cat("Arquivos gerados em:", OUT_DIR, "\n")
cat("- painel_final_vc_ai_macro.csv\n")
cat("- resultados_fe_estaticos_macro.csv\n")
cat("- resultados_dinamicos_fe_gmm_macro.csv\n")
cat("- diagnosticos_fe_macro.csv\n")
cat("- diagnosticos_gmm_macro.csv\n")
cat("- diff_in_hansen_macro.csv\n")
cat("- tabela_fe_diffgmm_sysgmm_macro.csv\n")
cat("- resultados_completos_vc_ai_FE_GMM_macro.xlsx\n")
cat("- pasta graficos com arquivos .png\n")

################################################################################
