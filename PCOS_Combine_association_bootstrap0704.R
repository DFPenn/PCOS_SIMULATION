library(readxl)
library(pROC)
library(dplyr)

# Read 
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

#create the menstrual group
data <- data %>%
  mutate(menstrual_group = sapply(menstrual_cycle2, menstrual_group))
print(table(data$menstrual_group))

# Define function to fit logistic regression and calculate OR, AUC
Fit_logistic_regression <- function(data, group1, group2, y_var, additional_vars) {
  
  # Filter the data to the specified groups
  filtered_data <- data %>%
    filter(menstrual_group %in% c(group1, group2)) %>%
    mutate(X = ifelse(menstrual_group == group1, 0, 1))
  
  # the formula for logistic regression
  formula_str <- paste(y_var, "~ X")
  if (length(additional_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(additional_vars, collapse = " + "))
  }
  formula <- as.formula(formula_str)
  
  # Fit the logistic regression model
  model <- glm(formula, data = filtered_data, family = binomial)
  
  # Calculate the Odds Ratio
  OR <- exp(coef(model)["X"])
  
  # p-value for X
  p_value <- summary(model)$coefficients["X", "Pr(>|z|)"]
  
  # Calculate the AUC
  roc_obj <- roc(filtered_data[[y_var]], fitted(model))
  AUC <- auc(roc_obj)
  
  return(list(OR = OR, p_value = p_value, AUC = AUC))
}

# Define groups, response variables, and additional variables sets
groups <- list(
  c("Eumeno", "Oligo"),
  c("Eumeno", "Ameno"),
  c("Oligo", "Ameno")
)

response_vars <- c("IR", "IFGIGT","dyslipidemia")
# "dyslipidemia"


additional_vars_sets <- list(
  list(),  # m1: no additional variables
  c("Age", "BMI", "WC"),  # m2: Age, BMI, WC
  c("Age", "BMI", "WC", "FAI")  # m3: Age, BMI, WC, FAI
)

# Set the number of bootstrap samples
n_bootstrap <- 100

results <- list()

# Loop through groups, response variables, and additional variables sets
for (group in groups) {
  group1 <- group[1]
  group2 <- group[2]
  
  for (y_var in response_vars) {
    for (i in 1:length(additional_vars_sets)) {
      
      additional_vars <- additional_vars_sets[[i]]
      model_id <- paste0("m", i)
      
      # Initialize vectors to store bootstrap results
      ORs <- numeric(n_bootstrap)
      p_values <- numeric(n_bootstrap)
      AUCs <- numeric(n_bootstrap)
      
      for (b in 1:n_bootstrap) {
        # Bootstrap sampling
        set.seed(b)
        bootstrap_sample <- data[sample(nrow(data), replace = TRUE), ]
        result <- Fit_logistic_regression(bootstrap_sample, group1, group2, y_var, additional_vars)
        
        ORs[b] <- result$OR
        p_values[b] <- result$p_value
        AUCs[b] <- result$AUC
      }
      
      # Calculate the mean of the bootstrap results
      mean_OR <- mean(ORs)
      mean_p_value <- mean(p_values)
      mean_AUC <- mean(AUCs)
      
      result_name <- paste(model_id, group1, group2, y_var, sep = "_")
      
      results[[result_name]] <- c(OR = mean_OR, p_value = mean_p_value, AUC = mean_AUC)
    }
  }
}

results_df <- as.data.frame(t(as.data.frame(results)))
print(results_df)
