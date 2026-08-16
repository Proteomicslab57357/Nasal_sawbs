# ============================================================
# Per-metabolite linear mixed model + FDR-corrected emmeans
# Design:
#   abundance ~ Status * Role + (1 | FamilyID)
#
# IMPORTANT:
#   This version assumes the input abundance matrix is ALREADY:
#   1) log2-transformed
#   2) NOT scaled
#
# Input files:
#   ExpressionMatrix_Log2_Transformed.csv   # metabolites in rows, samples in columns
#   metadata.csv                            # sample metadata
#
# Output files:
#   01_model_term_statistics.csv
#   02_emmeans_group_means.csv
#   03_emmeans_pairwise_contrasts.csv
#   04_significant_model_terms_FDR_lt_0.05.csv
#   05_significant_posthoc_contrasts_FDR_lt_0.05.csv
#   06_failed_metabolites.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# -----------------------------
# 1. User settings
# -----------------------------
expression_file <- "ExpressionMatrix_Log2_Transformed.csv"
metadata_file   <- "metadata.csv"
output_dir      <- "LMM_emmeans_results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Read input data
# -----------------------------
raw_df <- read.csv(expression_file, check.names = FALSE, stringsAsFactors = FALSE)
meta   <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

# Basic checks
if (!"Metabolite" %in% colnames(raw_df)) {
  stop("The first column of ExpressionMatrix_Log2_Transformed.csv must be named 'Metabolite' and contain metabolite names.")
}

required_meta_cols <- c("Sample", "Status", "Role", "FamilyID")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
if (length(missing_meta_cols) > 0) {
  stop(
    "metadata.csv is missing required columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

# -----------------------------
# 3. Convert metabolite table to long format
#    Input format:
#      rows    = metabolites
#      columns = samples
# -----------------------------
long_df <- raw_df %>%
  pivot_longer(
    cols = -Metabolite,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(meta, by = "Sample")

# Check join success
if (any(is.na(long_df$Status)) || any(is.na(long_df$Role)) || any(is.na(long_df$FamilyID))) {
  missing_samples <- long_df %>%
    filter(is.na(Status) | is.na(Role) | is.na(FamilyID)) %>%
    distinct(Sample) %>%
    pull(Sample)
  
  stop(
    "Some samples in ExpressionMatrix_Log2_Transformed.csv were not matched in metadata.csv: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Make sure abundance is numeric
long_df$Abundance <- as.numeric(long_df$Abundance)

if (any(!is.finite(long_df$Abundance), na.rm = TRUE)) {
  stop("Non-finite abundance values detected.")
}

# -----------------------------
# 4. Use input abundance directly
#    Input is already log2-transformed only
# -----------------------------
long_df <- long_df %>%
  mutate(
    abundance = Abundance
  )

# -----------------------------
# 5. Factor coding
#    Set explicit reference levels if present
# -----------------------------
long_df <- long_df %>%
  mutate(
    Status   = factor(Status),
    Role     = factor(Role),
    FamilyID = factor(FamilyID),
    Group    = if ("Group" %in% colnames(.)) factor(Group) else interaction(Status, Role, sep = ".")
  )

# Prefer healthy as reference if available
if (all(c("healthy", "pneumonia") %in% levels(long_df$Status))) {
  long_df$Status <- relevel(long_df$Status, ref = "healthy")
}

# Prefer mother as reference if available
if (all(c("mother", "infant") %in% levels(long_df$Role))) {
  long_df$Role <- relevel(long_df$Role, ref = "mother")
}

# Check model identifiability
if (nlevels(long_df$Status) < 2) {
  stop("Status has fewer than 2 levels in metadata.csv. Cannot fit the model.")
}
if (nlevels(long_df$Role) < 2) {
  stop("Role has fewer than 2 levels in metadata.csv. Cannot fit the model.")
}
if (nlevels(long_df$FamilyID) < 2) {
  warning("FamilyID has fewer than 2 levels. The random effect may not be estimable.")
}

# -----------------------------
# 6. Helper function to fit one metabolite
# -----------------------------
fit_one_metabolite <- function(df_met) {
  metabolite_name <- unique(df_met$Metabolite)
  
  # remove missing rows just in case
  df_met <- df_met %>%
    filter(
      is.finite(abundance),
      !is.na(Status),
      !is.na(Role),
      !is.na(FamilyID)
    )
  
  # Need enough data
  if (nrow(df_met) < 4) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = "Too few observations")
    ))
  }
  
  # Must have all factor levels represented for this metabolite
  if (nlevels(droplevels(df_met$Status)) < 2) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = "Status has <2 levels")
    ))
  }
  
  if (nlevels(droplevels(df_met$Role)) < 2) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = "Role has <2 levels")
    ))
  }
  
  # Fit model
  fit <- tryCatch(
    lmer(abundance ~ Status * Role + (1 | FamilyID), data = df_met, REML = FALSE),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = paste("Model fit error:", fit$message))
    ))
  }
  
  # Type III ANOVA from lmerTest
  anova_tab <- tryCatch(
    anova(fit, type = 3),
    error = function(e) e
  )
  
  if (inherits(anova_tab, "error")) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = paste("ANOVA extraction error:", anova_tab$message))
    ))
  }
  
  anova_df <- as.data.frame(anova_tab)
  anova_df$Term <- rownames(anova_df)
  rownames(anova_df) <- NULL
  
  # Keep only fixed-effect terms of interest
  p_col     <- intersect(c("Pr(>F)", "Pr(>Chisq)"), colnames(anova_df))
  f_col     <- intersect(c("F value", "Chisq"), colnames(anova_df))
  numdf_col <- intersect(c("NumDF", "Df"), colnames(anova_df))
  dendf_col <- intersect(c("DenDF"), colnames(anova_df))
  
  if (length(p_col) == 0) {
    return(list(
      term_stats = tibble(),
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = "No p-value column found in ANOVA table")
    ))
  }
  
  term_stats <- anova_df %>%
    filter(Term %in% c("Status", "Role", "Status:Role")) %>%
    transmute(
      Metabolite = metabolite_name,
      Term = Term,
      NumDF = if (length(numdf_col) > 0) .data[[numdf_col[1]]] else NA_real_,
      DenDF = if (length(dendf_col) > 0) .data[[dendf_col[1]]] else NA_real_,
      Statistic = if (length(f_col) > 0) .data[[f_col[1]]] else NA_real_,
      P_value = .data[[p_col[1]]]
    )
  
  # EMMs for Status x Role combinations
  emm_obj <- tryCatch(
    emmeans(fit, ~ Status * Role),
    error = function(e) e
  )
  
  if (inherits(emm_obj, "error")) {
    return(list(
      term_stats = term_stats,
      emm_means = tibble(),
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = paste("emmeans error:", emm_obj$message))
    ))
  }
  
  emm_means <- tryCatch(
    summary(emm_obj, infer = TRUE) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(
        Metabolite = metabolite_name,
        Mean_OriginalScale = 2^emmean,
        Lower_CL_OriginalScale = 2^lower.CL,
        Upper_CL_OriginalScale = 2^upper.CL
      ) %>%
      rename(
        EMM = emmean,
        SE = SE,
        df = df,
        Lower_CL = lower.CL,
        Upper_CL = upper.CL
      ) %>%
      select(
        Metabolite, Status, Role,
        EMM, SE, df, Lower_CL, Upper_CL,
        Mean_OriginalScale, Lower_CL_OriginalScale, Upper_CL_OriginalScale
      ),
    error = function(e) tibble(
      Metabolite = metabolite_name,
      Error = paste("EMM summary error:", e$message)
    )
  )
  
  # All pairwise contrasts among the four Status x Role combinations
  contrast_obj <- tryCatch(
    contrast(emm_obj, method = "pairwise", adjust = "none"),
    error = function(e) e
  )
  
  if (inherits(contrast_obj, "error")) {
    return(list(
      term_stats = term_stats,
      emm_means = emm_means,
      contrasts = tibble(),
      failed = tibble(Metabolite = metabolite_name, Reason = paste("contrast error:", contrast_obj$message))
    ))
  }
  
  contrast_df <- tryCatch(
    summary(contrast_obj, infer = TRUE) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(
        Metabolite = metabolite_name,
        FoldChange = 2^estimate,
        FoldChange_Lower_CL = 2^lower.CL,
        FoldChange_Upper_CL = 2^upper.CL,
        Log2FoldChange = estimate
      ) %>%
      rename(
        Contrast = contrast,
        Estimate = estimate,
        SE = SE,
        df = df,
        t_ratio = t.ratio,
        P_value = p.value,
        Lower_CL = lower.CL,
        Upper_CL = upper.CL
      ) %>%
      select(
        Metabolite, Contrast,
        Estimate, Log2FoldChange, FoldChange,
        SE, df, t_ratio, P_value,
        Lower_CL, Upper_CL,
        FoldChange_Lower_CL, FoldChange_Upper_CL
      ),
    error = function(e) tibble(
      Metabolite = metabolite_name,
      Error = paste("Contrast summary error:", e$message)
    )
  )
  
  list(
    term_stats = term_stats,
    emm_means = emm_means,
    contrasts = contrast_df,
    failed = tibble()
  )
}

# -----------------------------
# 7. Run model for every metabolite
# -----------------------------
metabolite_list <- unique(long_df$Metabolite)

results_list <- lapply(metabolite_list, function(met) {
  df_met <- long_df %>% filter(Metabolite == met)
  fit_one_metabolite(df_met)
})

# -----------------------------
# 8. Combine outputs
# -----------------------------
term_stats_all <- bind_rows(lapply(results_list, `[[`, "term_stats"))
emm_means_all  <- bind_rows(lapply(results_list, `[[`, "emm_means"))
contrasts_all  <- bind_rows(lapply(results_list, `[[`, "contrasts"))
failed_all     <- bind_rows(lapply(results_list, `[[`, "failed"))

# -----------------------------
# 9. FDR correction
# -----------------------------
# 9a. FDR for model terms across metabolites
if (nrow(term_stats_all) > 0) {
  term_stats_all <- term_stats_all %>%
    group_by(Term) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    ungroup()
}

# 9b. FDR for pairwise post-hoc contrasts across metabolites within each contrast
if (nrow(contrasts_all) > 0 && "Contrast" %in% colnames(contrasts_all)) {
  contrasts_all <- contrasts_all %>%
    group_by(Contrast) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    ungroup()
}

# -----------------------------
# 10. Write outputs
# -----------------------------
write.csv(
  term_stats_all,
  file = file.path(output_dir, "01_model_term_statistics.csv"),
  row.names = FALSE
)

write.csv(
  emm_means_all,
  file = file.path(output_dir, "02_emmeans_group_means.csv"),
  row.names = FALSE
)

write.csv(
  contrasts_all,
  file = file.path(output_dir, "03_emmeans_pairwise_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  term_stats_all %>% filter(FDR < 0.05),
  file = file.path(output_dir, "04_significant_model_terms_FDR_lt_0.05.csv"),
  row.names = FALSE
)

write.csv(
  contrasts_all %>% filter(FDR < 0.05),
  file = file.path(output_dir, "05_significant_posthoc_contrasts_FDR_lt_0.05.csv"),
  row.names = FALSE
)

write.csv(
  failed_all,
  file = file.path(output_dir, "06_failed_metabolites.csv"),
  row.names = FALSE
)

# -----------------------------
# 11. Console summary
# -----------------------------
cat("\n============================================================\n")
cat("Analysis completed.\n")
cat("Output directory:", output_dir, "\n")
cat("Input data assumption: already log2-transformed and NOT scaled\n")
cat("Input file:", expression_file, "\n")
cat("Total metabolites attempted:", length(metabolite_list), "\n")
cat(
  "Successfully analyzed (term stats available):",
  dplyr::n_distinct(term_stats_all$Metabolite),
  "\n"
)
cat("Failed metabolites:", nrow(failed_all), "\n")

if (nrow(term_stats_all) > 0) {
  cat("\nSignificant model terms at FDR < 0.05:\n")
  print(
    term_stats_all %>%
      filter(FDR < 0.05) %>%
      count(Term, name = "N_significant")
  )
}

if (nrow(contrasts_all) > 0) {
  cat("\nSignificant post-hoc contrasts at FDR < 0.05:\n")
  print(
    contrasts_all %>%
      filter(FDR < 0.05) %>%
      count(Contrast, name = "N_significant")
  )
}
cat("============================================================\n")