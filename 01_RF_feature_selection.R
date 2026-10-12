# Random Forest feature selection for island and riverine ARG profiles
# Input: ARG_subtype_abundance.xlsx (sheets: Water and Sheet1)
# The first two columns contain metadata; remaining numeric columns contain
# ARG subtype abundances (copies/cell).

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(randomForest)
})

input_file <- "ARG_subtype_abundance.xlsx"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

water <- read_excel(input_file, sheet = "Water")
references <- read_excel(input_file, sheet = "Sheet1")
group_column <- names(references)[1]

island <- water %>% select(-c(1, 2)) %>% select(where(is.numeric))
river <- references %>%
  filter(.data[[group_column]] == "Yangtze River") %>%
  select(-c(1, 2)) %>% select(where(is.numeric))

island[is.na(island)] <- 0
river[is.na(river)] <- 0
common_features <- intersect(names(island), names(river))
island <- island[, common_features, drop = FALSE]
river <- river[, common_features, drop = FALSE]

cat(sprintf("Island samples: %d; riverine samples: %d\n", nrow(island), nrow(river)))

# Screen ARG subtypes for statistically significant enrichment in island water.
# The pseudocount, Wilcoxon parameters, and filtering criteria are unchanged
# from the original analysis script.
screening <- bind_rows(lapply(common_features, function(feature) {
  island_values <- island[[feature]]
  river_values <- river[[feature]]
  island_mean <- mean(island_values)
  river_mean <- mean(river_values)
  island_prevalence <- mean(island_values > 0)
  river_prevalence <- mean(river_values > 0)
  p_value <- tryCatch(
    wilcox.test(island_values, river_values, exact = FALSE)$p.value,
    error = function(e) 1
  )
  data.frame(
    Feature = feature,
    mean_Chongming = island_mean,
    mean_Yangtze = river_mean,
    prev_Chongming = island_prevalence,
    prev_Yangtze = river_prevalence,
    fold_change = (island_mean + 1e-10) / (river_mean + 1e-10),
    p_value = p_value
  )
})) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))

candidates <- screening %>% filter(
  mean_Chongming > mean_Yangtze,
  fold_change >= 1.5,
  prev_Chongming > prev_Yangtze,
  p_adj < 0.05
)

if (nrow(candidates) == 0) stop("No ARG subtypes passed the predefined enrichment criteria.")

# Train the final 1,000-tree RF model and rank retained subtypes by Gini importance.
candidate_features <- candidates$Feature
training_x <- bind_rows(
  island[, candidate_features, drop = FALSE],
  river[, candidate_features, drop = FALSE]
)
training_x[is.na(training_x)] <- 0
training_y <- as.factor(c(
  rep("Chongming", nrow(island)), rep("Yangtze", nrow(river))
))
training_x <- training_x[, colSums(training_x) > 0, drop = FALSE]

set.seed(123)
rf_model <- randomForest(
  x = training_x, y = training_y, ntree = 1000, importance = TRUE
)
importance <- data.frame(
  Feature = rownames(rf_model$importance),
  Gini = rf_model$importance[, "MeanDecreaseGini"]
) %>% arrange(desc(Gini))

ranked_candidates <- candidates %>%
  inner_join(importance, by = "Feature") %>%
  arrange(desc(Gini))
top20 <- head(ranked_candidates, 20)

write.csv(screening, file.path(output_dir, "ARG_subtype_screening.csv"), row.names = FALSE)
write.csv(ranked_candidates, file.path(output_dir, "RF_Gini_ranked_candidates.csv"), row.names = FALSE)
write.csv(top20, file.path(output_dir, "RF_top20_ARG_subtypes.csv"), row.names = FALSE)
saveRDS(rf_model, file.path(output_dir, "RF_final_model.rds"))

cat(sprintf("Enriched candidates: %d; selected features: %d\n", nrow(candidates), nrow(top20)))
print(top20[, c("Feature", "Gini", "p_adj", "fold_change")])
