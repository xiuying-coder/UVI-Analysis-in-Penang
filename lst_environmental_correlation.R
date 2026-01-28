# ============================================================
# Linear correlation analysis between LST and
# UVI, Albedo, UV_DRS, NDVI, and NDBI
#
# Input:
#   data/raw/Penang_ClearSky_MultiVar_LST_UV_2015_2024.csv
#
# Output:
#   1. Correlation statistics table (CSV)
#   2. Correlation heatmap (PNG)
#   3. Scatter plots with regression lines by land-use type (PNG)
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "reshape2", "broom")

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
}
lapply(required_pkgs, library, character.only = TRUE)


# ------------------------------------------------------------
# 2. Define project structure and paths
# ------------------------------------------------------------
# The working directory should be set to the project root
# Example:
# setwd("LST_UV_Project")

project_dir <- getwd()

data_dir    <- file.path(project_dir, "data", "raw")
results_dir <- file.path(project_dir, "results")
fig_dir     <- file.path(results_dir, "figures")

if (!dir.exists(results_dir)) dir.create(results_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# ------------------------------------------------------------
# 3. Read input data
# ------------------------------------------------------------
file_path <- file.path(
  data_dir,
  "Penang_ClearSky_MultiVar_LST_UV_2015_2024.csv"
)

data <- read.csv(file_path, stringsAsFactors = FALSE)

data <- data %>%
  filter(
    !is.na(lst),
    !is.na(uvi),
    !is.na(albedo),
    !is.na(uv_drs),
    !is.na(landuse_type)
  )

cat("Data loaded successfully:", nrow(data), "rows retained\n")


# ------------------------------------------------------------
# 4. Function to calculate correlation statistics
# ------------------------------------------------------------
get_corr_stats <- function(df, xvar, yvar = "lst") {
  
  if (unique(df$landuse_type) == "Water Bodies" &&
      xvar %in% c("ndvi", "ndbi")) {
    return(data.frame(
      variable = xvar,
      r = NA,
      p = NA,
      R2 = NA
    ))
  }
  
  if (sum(!is.na(df[[xvar]])) < 5 ||
      length(unique(df[[xvar]][!is.na(df[[xvar]])])) < 5) {
    return(data.frame(
      variable = xvar,
      r = NA,
      p = NA,
      R2 = NA
    ))
  }
  
  test  <- cor.test(df[[xvar]], df[[yvar]], use = "complete.obs")
  model <- lm(df[[yvar]] ~ df[[xvar]])
  
  data.frame(
    variable = xvar,
    r  = as.numeric(test$estimate),
    p  = test$p.value,
    R2 = summary(model)$r.squared
  )
}


# ------------------------------------------------------------
# 5. Correlation analysis by land-use type
# ------------------------------------------------------------
vars <- c("uvi", "albedo", "uv_drs", "ndvi", "ndbi")

corr_results <- data %>%
  group_by(landuse_type) %>%
  do({
    res_list <- lapply(vars, function(v) get_corr_stats(., v))
    do.call(rbind, res_list)
  }) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )


# ------------------------------------------------------------
# 6. Save correlation results
# ------------------------------------------------------------
out_csv <- file.path(
  results_dir,
  "LST_Linear_Correlation_Results.csv"
)

corr_results_export <- corr_results %>%
  mutate(
    r  = ifelse(is.na(r),  "-", sprintf("%.2f", r)),
    p  = ifelse(is.na(p),  "-", format.pval(p, digits = 2, eps = 0.001)),
    R2 = ifelse(is.na(R2), "-", sprintf("%.2f", R2))
  )

write.csv(corr_results_export, out_csv, row.names = FALSE)

cat("Correlation table saved to:", out_csv, "\n")


# ------------------------------------------------------------
# 7. Correlation heatmap
# ------------------------------------------------------------
corr_melt <- corr_results %>%
  filter(!is.na(r)) %>%
  melt(
    id.vars = c("landuse_type", "variable"),
    measure.vars = "r",
    value.name = "corr_value"
  )

p_heatmap <- ggplot(corr_melt,
                    aes(x = variable, y = landuse_type, fill = corr_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", corr_value)), size = 4) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Pearson r"
  ) +
  labs(
    title = "Correlation between LST and Environmental Variables",
    x = "Variable",
    y = "Land Use Type"
  ) +
  theme_minimal(base_size = 14)

heatmap_file <- file.path(fig_dir, "LST_Correlation_Heatmap.png")
ggsave(heatmap_file, p_heatmap, width = 10, height = 6, dpi = 300)

cat("Heatmap saved to:", heatmap_file, "\n")


# ------------------------------------------------------------
# 8. Scatter plots with regression lines
# ------------------------------------------------------------
for (var in vars) {
  
  plot_data <- data
  
  if (var %in% c("ndvi", "ndbi")) {
    plot_data <- plot_data %>%
      filter(landuse_type != "Water Bodies")
  }
  
  p_scatter <- ggplot(plot_data, aes(x = .data[[var]], y = lst)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~landuse_type, scales = "free", ncol = 2) +
    labs(
      title = paste("LST vs", toupper(var), "by Land Use Type"),
      x = toupper(var),
      y = "LST"
    ) +
    theme_bw(base_size = 14)
  
  file_name <- file.path(
    fig_dir,
    paste0("LST_vs_", toupper(var), "_by_LandUse.png")
  )
  
  ggsave(file_name, p_scatter, width = 10, height = 8, dpi = 300)
  cat("Scatter plot saved:", file_name, "\n")
}


cat("Analysis completed successfully\n")
