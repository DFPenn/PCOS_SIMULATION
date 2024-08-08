library(readxl)
library(pROC)
library(dplyr)
library(boot)
library(writexl)  # For writing to Excel

# Read data
data <- read_excel("E:/table/IDI_NRI_article/PCOS/PCOSdata154622.xlsx")

# Define the grouping function
menstrual_group <- function(cycle) {
  if (cycle >= 26 & cycle <= 36) {
    return("Eumeno")
  } else if (cycle >= 37 & cycle < 90) {
    return("Oligo")
  } else {
    return("Ameno")
  }
}

# Create the menstrual group
data <- data %>%
  mutate(menstrual_group = sapply(menstrual_cycle2, menstrual_group))
print(table(data$menstrual_group))

# Define function to fit logistic regression and calculate OR, AUC
Fit_logistic_regression <- function(data, group1, group2, outcome) {
  
  # Filter the data to the specified groups
  filtered_data <- data %>%
    filter(menstrual_group %in% c(group1, group2)) %>%
    mutate(X = ifelse(menstrual_group == group1, 0, 1))
  
  # Baseline model
  baseline_model <- glm(as.formula(paste(outcome, "~ Age + BMI + WC + FAI")), data = filtered_data, family = binomial)
  baseline_pred <- predict(baseline_model, filtered_data, type = "response")
  baseline_roc <- roc(filtered_data[[outcome]], baseline_pred)
  baseline_auc <- auc(baseline_roc)
  
  # New model
  new_model <- glm(as.formula(paste(outcome, "~ Age + BMI + WC + FAI + X")), data = filtered_data, family = binomial)
  new_pred <- predict(new_model, filtered_data, type = "response")
  new_roc <- roc(filtered_data[[outcome]], new_pred)
  new_auc <- auc(new_roc)
  
  # DeLong Test
  delong_test <- roc.test(baseline_roc, new_roc, method = "delong")
  
  # Calculate IDI and NRI
  roc.p <- pROC::roc(filtered_data[[outcome]], baseline_pred)
  cutoffB <- roc.p$thresholds[which.max(roc.p$sensitivities + roc.p$specificities)]
  
  c1 <- cut(baseline_pred, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  c2 <- cut(new_pred, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  
  categ <- improveProb(x1 = as.numeric(c11) * (1 / (length(levels(c11)))), 
                       x2 = as.numeric(c22) * (1 / (length(levels(c22)))), 
                       y = filtered_data[[outcome]])
  
  conti <- improveProb(x1 = baseline_pred, x2 = new_pred, y = filtered_data[[outcome]])
  
  idi <- round(conti$idi, 4)
  categ_nri <- round(categ$nri, 4)
  conti_nri <- round(conti$nri, 4)
  
  return(c(baseline_auc, new_auc, idi, categ_nri, conti_nri, delong_test$p.value))
}

# Define bootstrap function
bootstrap_fun <- function(data, indices, group1, group2, outcome) {
  d <- data[indices, ]  
  Fit_logistic_regression(d, group1, group2, outcome)
}

# Execute Bootstrap analysis
repeat_analysis <- function(data, group1, group2, outcome, n_iter = 1000) {
  # Execute Bootstrap
  results <- boot(data = data, statistic = function(data, indices) bootstrap_fun(data, indices, group1, group2, outcome), R = n_iter)
  
  mean_results <- apply(results$t, 2, mean)
  sd_results <- apply(results$t, 2, sd)
  
  mean_baseline_auc <- mean_results[1]
  mean_new_auc <- mean_results[2]
  mean_idi <- mean_results[3]
  mean_categ_nri <- mean_results[4]
  mean_conti_nri <- mean_results[5]
  delong_significant <- sum(results$t[, 6] < 0.05)
  
  # Save each bootstrap result to an Excel file
  results_df <- as.data.frame(results$t) %>%
    mutate(group1 = group1, group2 = group2)
  
  #colnames(results_df) <- c("Baseline_AUC", "New_AUC", "IDI", "Categorical_NRI", "Continuous_NRI", "DeLong_p_value", "group1", "group2")
  
  # Output to Excel file
  #write_xlsx(results_df, paste0("E:/table/IDI_NRI_article/PCOS/Bootstrap_", outcome, "_", group1, "_vs_", group2, "_results.xlsx"))
  
  return(list(
    mean_baseline_auc = mean_baseline_auc,
    sd_baseline_auc = sd_results[1],
    mean_new_auc = mean_new_auc,
    sd_new_auc = sd_results[2],
    mean_idi = mean_idi,
    sd_idi = sd_results[3],
    mean_categ_nri = mean_categ_nri,
    sd_categ_nri = sd_results[4],
    mean_conti_nri = mean_conti_nri,
    sd_conti_nri = sd_results[5],
    delong_significant = delong_significant
  ))
}

# Define combinations
combinations <- list(c("Eumeno", "Oligo"), c("Eumeno", "Ameno"), c("Oligo", "Ameno"))

# Define outcome list
outcomes <- c("IR", "IFGIGT", "dyslipidemia")

# Loop through each combination and outcome
results <- list()
for (outcome in outcomes) {
  for (comb in combinations) {
    group1 <- comb[1]
    group2 <- comb[2]
    res <- repeat_analysis(data, group1, group2, outcome)
    results[[paste(outcome, group1, "vs", group2)]] <- res
  }
}

# Print results
for (result in names(results)) {
  cat(paste("Combination:", result, "\n"))
  cat(paste("Mean Baseline AUC:", results[[result]]$mean_baseline_auc, "\n"))
  cat(paste("SD Baseline AUC:", results[[result]]$sd_baseline_auc, "\n"))
  cat(paste("Mean New Model AUC:", results[[result]]$mean_new_auc, "\n"))
  cat(paste("SD New Model AUC:", results[[result]]$sd_new_auc, "\n"))
  cat(paste("Mean IDI:", results[[result]]$mean_idi, "\n"))
  cat(paste("SD IDI:", results[[result]]$sd_idi, "\n"))
  cat(paste("Mean Categorical NRI:", results[[result]]$mean_categ_nri, "\n"))
  cat(paste("SD Categorical NRI:", results[[result]]$sd_categ_nri, "\n"))
  cat(paste("Mean Continuous NRI:", results[[result]]$mean_conti_nri, "\n"))
  cat(paste("SD Continuous NRI:", results[[result]]$sd_conti_nri, "\n"))
  cat(paste("Proportion of significant DeLong tests:", results[[result]]$delong_significant, "\n"))
  cat("\n")
}











#make this example reproducible
set.seed(10)

#load boot library
library(boot)

#define dataset
x <- c(12, 14, 14, 15, 18, 21, 25, 29, 32, 35)

#define function to calculate mean
meanFunc <- function(x,i){mean(x[i])}

#calculate standard error using 100 bootstrapped samples
boot(x, meanFunc, 100)



# Make this example reproducible
set.seed(10)

# Load boot library
library(boot)

# Define dataset
x <- c(12, 14, 14, 15, 18, 21, 25, 29, 32, 35)

# Define function to calculate mean
meanFunc <- function(x, i){mean(x[i])}

# Calculate standard error using 100 bootstrapped samples
boot_results <- boot(x, meanFunc, 100)

# Extract the standard error
B_error <- sd(boot_results$t)

# Print the standard error
print(B_error)


