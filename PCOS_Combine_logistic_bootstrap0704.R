library(readxl)
library(pROC)
library(dplyr)
library(boot)

# Read data
data <- read_excel("E:/table/IDI_NRI_article/PCOS/PCOSdata154622.xlsx")

# Define the grouping function
menstrual_group <- function(cycle) {
  if (cycle >= 26 & cycle <= 35) {
    return("Eumeno")
  } else if (cycle >= 36 & cycle < 90) {
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
Fit_logistic_regression <- function(data, indices, group1, group2, y_var, additional_vars) {
  # Subset the data based on bootstrap indices
  d <- data[indices, ]
  
  # Filter the data to the specified groups
  filtered_data <- d %>%
    filter(menstrual_group %in% c(group1, group2)) %>%
    mutate(X = ifelse(menstrual_group == group1, 0, 1))
  
  # The formula for logistic regression
  formula_str <- paste(y_var, "~ X")
  if (length(additional_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(additional_vars, collapse = " + "))
  }
  formula <- as.formula(formula_str)
  
  # Fit the logistic regression model
  model <- glm(formula, data = filtered_data, family = binomial)
  
  # Calculate the Odds Ratio
  OR <- exp(coef(model)["X"])
  
  # P-value for X
  p_value <- summary(model)$coefficients["X", "Pr(>|z|)"]
  
  # Calculate the AUC
  roc_obj <- roc(filtered_data[[y_var]], fitted(model))
  AUC <- auc(roc_obj)
  
  return(c(OR, p_value, AUC))
}

# Define groups, response variables, and additional variables sets
groups <- list(
  c("Eumeno", "Oligo"),
  c("Eumeno", "Ameno"),
  c("Oligo", "Ameno")
)

response_vars <- c("IR", "IFGIGT", "dyslipidemia")

additional_vars_sets <- list(
  list(),  # m1: no additional variables
  c("Age", "BMI", "WC"),  # m2: Age, BMI, WC
  c("Age", "BMI", "WC", "FAI")  # m3: Age, BMI, WC, FAI
)

results <- list()

# Loop through groups, response variables, and additional variables sets
for (group in groups) {
  group1 <- group[1]
  group2 <- group[2]
  
  for (y_var in response_vars) {
    for (i in 1:length(additional_vars_sets)) {
      additional_vars <- additional_vars_sets[[i]]
      model_id <- paste0("m", i)
      
      # Perform bootstrap
      boot_result <- boot(data = data, statistic = function(data, indices) {
        Fit_logistic_regression(data, indices, group1, group2, y_var, additional_vars)
      }, R = 100)
      
      OR_mean <- mean(boot_result$t[, 1])
      p_value_sum <- sum(boot_result$t[, 2] < 0.05)  # count  p-values < 0.05
      AUC_mean <- mean(boot_result$t[, 3])
      
      result_name <- paste(model_id, group1, group2, y_var, sep = "_")
      
      results[[result_name]] <- c(OR = OR_mean, p_value = p_value_sum, AUC = AUC_mean)
    }
  }
}

results_df <- as.data.frame(t(as.data.frame(results)))
print(results_df)
write.csv(results_df, "E:/table/IDI_NRI_article/PCOS/ORAUCP_Bootstrap_results.csv", row.names = TRUE)

