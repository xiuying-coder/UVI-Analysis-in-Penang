# ==========================================================
# LST Linear Contribution Analysis
# Standardized regression with nonlinear Albedo
# NDVI and NDBI included
#
# Reproducible version using relative paths
# ==========================================================

# ----------------------------------------------------------
# 1. Load required packages
# ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)

# ----------------------------------------------------------
# 2. Define project structure and paths
# ----------------------------------------------------------

project_dir <- getwd()

data_dir    <- file.path(project_dir, "data", "raw")
results_dir <- file.path(project_dir, "results")
fig_dir     <- file.path(results_dir, "figures")

if (!dir.exists(results_dir)) dir.create(results_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# ----------------------------------------------------------
# 3. Read input data
# ----------------------------------------------------------
input_file <- file.path(
  data_dir,
  "Penang_ClearSky_MultiVar_LST_UV_2015_2024.csv"
)

data <- read.csv(input_file, stringsAsFactors = FALSE)

# ----------------------------------------------------------
# 4. Data preprocessing
# ----------------------------------------------------------
data <- data %>%
  mutate(UV_DRS_scaled = uv_drs / 100) %>%
  filter(
    landuse_type != "Water Bodies",
    !is.na(lst),
    !is.na(uvi),
    !is.na(albedo),
    !is.na(ndvi),
    !is.na(ndbi),
    !is.na(UV_DRS_scaled),
    !is.na(landuse_type)
  )

# ----------------------------------------------------------
# 5. Standardization function
# ----------------------------------------------------------
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


# ----------------------------------------------------------
# 6. Output files
# ----------------------------------------------------------
formula_txt <- file.path(
  results_dir,
  "LST_Model_Formulas_AlbedoNonlinear_NDVI_NDBI.txt"
)

if (file.exists(formula_txt)) file.remove(formula_txt)
file.create(formula_txt)

contrib_results <- data.frame()

# ----------------------------------------------------------
# 7. Modeling by land-use type
# ----------------------------------------------------------
for (lu in unique(data$landuse_type)) {
  
  subset <- data %>% filter(landuse_type == lu)
  
  if (nrow(subset) > 100) {
    
    subset_std <- subset %>%
      mutate(across(
        c(lst, uvi, albedo, ndvi, ndbi, UV_DRS_scaled),
        standardize
      )) %>%
      mutate(albedo2 = albedo^2)
    
    model <- lm(
      lst ~ uvi + albedo + albedo2 + UV_DRS_scaled + ndvi + ndbi,
      data = subset_std
    )
    
    sm <- summary(model)
    betas <- sm$coefficients[-1, 1]
    
    formula_text <- paste0(
      "Land Use Type: ", lu, "\n",
      "Standardized model:\n",
      "LST* = ",
      round(betas["uvi"], 3), " * UVI + ",
      round(betas["albedo"], 3), " * Albedo + ",
      round(betas["albedo2"], 3), " * Albedo^2 + ",
      round(betas["UV_DRS_scaled"], 3), " * UV_DRS + ",
      round(betas["ndvi"], 3), " * NDVI + ",
      round(betas["ndbi"], 3), " * NDBI\n",
      "R2 = ", round(sm$r.squared, 3), "\n",
      "---------------------------------------------\n"
    )
    
    write(formula_text, file = formula_txt, append = TRUE)
    
    contrib_df <- data.frame(
      LandUse = lu,
      Variable = c("UVI", "Albedo", "UV_DRS", "NDVI", "NDBI"),
      Beta = c(
        betas["uvi"],
        betas["albedo"] + betas["albedo2"],
        betas["UV_DRS_scaled"],
        betas["ndvi"],
        betas["ndbi"]
      )
    )
    
    abs_beta <- abs(contrib_df$Beta)
    contrib_df$Contribution <- 100 * abs_beta / sum(abs_beta)
    contrib_df$R2 <- sm$r.squared
    
    contrib_results <- rbind(contrib_results, contrib_df)
  }
}

# ----------------------------------------------------------
# 8. Save contribution table
# ----------------------------------------------------------
csv_output <- file.path(
  results_dir,
  "LST_Contribution_By_LandUse_AlbedoNonlinear_NDVI_NDBI.csv"
)

write.csv(contrib_results, csv_output, row.names = FALSE)


# ----------------------------------------------------------
# 9. Define global variable order based on Urbanized Areas
# ----------------------------------------------------------
urban_order <- contrib_results %>%
  filter(LandUse == "Urbanized Areas") %>%
  arrange(desc(Contribution)) %>%
  pull(Variable)

contrib_results$Variable <- factor(
  contrib_results$Variable,
  levels = urban_order
)

# ----------------------------------------------------------
# 10. Plot relative contributions
# ----------------------------------------------------------
p <- ggplot(
  contrib_results,
  aes(x = Variable, y = Contribution, fill = Variable)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = sprintf("%.1f%%", Contribution)),
    vjust = -0.3,
    size = 4.5
  ) +
  facet_wrap(~ LandUse, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 16) +
  labs(
    title = "Relative Contribution of Factors to LST (NDVI and NDBI)",
    x = "Factor (ordered by Urbanized Areas)",
    y = "Relative Contribution (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  expand_limits(y = max(contrib_results$Contribution) * 1.15)

ggsave(
  filename = file.path(
    fig_dir,
    "LST_Relative_Contribution_By_LandUse_NDVI_NDBI_UrbanOrdered.png"
  ),
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Linear contribution analysis completed successfully\n")
