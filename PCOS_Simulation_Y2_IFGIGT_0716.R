library(MASS)
library(pROC)
library(Hmisc)
library(openxlsx)

# Function to generate data
generate_data <- function(NN, beta0, betaX, gamaZ, seed, Z_type = "Z1") {
  
  set.seed(seed)
  
  # Simulated age
  X1 <- rnorm(n = NN, mean = 27, sd = 4.6)
  X1 <- floor(X1)
  X1 <- pmax(X1, min(X1[X1 >= 0]))
  
  # Simulated BMI
  X2 <- rnorm(n = NN, mean = 23, sd = 7)
  X2 <- floor(X2)
  X2 <- pmax(X2, min(X2[X2 >= 0]))
  
  # Simulated WC
  X3 <- rnorm(n = NN, mean = 77, sd = 10)
  X3 <- floor(X3)
  X3 <- pmax(X3, min(X3[X3 >= 0]))
  
  # Simulated FAI (对数正态分布)
  X4 <- rlnorm(n = NN, meanlog = 1.6, sdlog = 0.5)
  X4 <- floor(X4)
  X4 <- pmax(X4, min(X4[X4 >= 0]))
  
  # Simulated Z based on Z_type
  if (Z_type == "Z1") {
    Z <- rbinom(n = NN, size = 1, prob = 0.6)
  } else if (Z_type == "Z2") {
    Z <- rpois(n = NN, lambda = 2)
  } else if (Z_type == "Z3") {
    Z <- sample(1:4, size = NN, replace = TRUE, prob = c(0.15, 0.35, 0.35, 0.15))
  } else {
    stop("Invalid Z_type")
  }
  
  X <- cbind(1, X1, X2, X3, X4, Z)
  y <- rbinom(NN, size = 1, prob = exp(X %*% c(beta0, betaX, gamaZ)) / (1 + exp(X %*% c(beta0, betaX, gamaZ))))
  data.frame(y, X1, X2, X3, X4, Z)
}

# Function to calculate AUC, IDI, NRI, and Delong test
calculate_metrics <- function(dat) {
  m0 <- glm(y ~ X1 + X2 + X3 + X4, family = binomial, data = dat)
  m1 <- glm(y ~ X1 + X2 + X3 + X4 + Z, family = binomial, data = dat)
  
  yhat0 <- m0$fitted.values
  yhat1 <- m1$fitted.values
  
  roc_model0 <- roc(dat$y, yhat0)
  roc_model1 <- roc(dat$y, yhat1)
  
  AUC0 <- roc_model0$auc
  AUC1 <- roc_model1$auc
  
  delong_test_result <- roc.test(roc_model0, roc_model1)
  Delong_p_value <- delong_test_result$p.value
  
  roc.p <- pROC::roc(dat$y, yhat0)
  cutoffB <- roc.p$thresholds[which.max(roc.p$sensitivities + roc.p$specificities)]
  
  c1 <- cut(yhat0, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  c2 <- cut(yhat1, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  
  categ <- improveProb(x1 = as.numeric(c11) * (1 / (length(levels(c11)))),
                       x2 = as.numeric(c22) * (1 / (length(levels(c22)))), y = dat[, 1])
  
  conti <- improveProb(x1 = yhat0, x2 = yhat1, y = dat[, 1])
  
  IDI_model <- round(conti$idi, 4)
  NRI_categ_model <- round(categ$nri, 4)
  NRI_conti_model <- round(conti$nri, 4)
  
  list(AUC0 = AUC0, AUC1 = AUC1, Delong_p_value = Delong_p_value,
       IDI_model = IDI_model, NRI_categ_model = NRI_categ_model, NRI_conti_model = NRI_conti_model)
}

# Main simulation
simulate_metrics <- function(iterations, NN, beta0, betaX, gamaZ, Z_type) {
  results <- data.frame(
    iteration = 1:iterations,
    AUC0 = numeric(iterations),
    AUC1 = numeric(iterations),
    IDI = numeric(iterations),
    NRI_categ = numeric(iterations),
    NRI_conti = numeric(iterations),
    Delong_p_value = numeric(iterations)
  )
  
  AUC0_list <- numeric(iterations)
  AUC1_list <- numeric(iterations)
  IDI_list <- numeric(iterations)
  NRI_categ_list <- numeric(iterations)
  NRI_conti_list <- numeric(iterations)
  Delong_significant_count <- 0
  
  
  for (i in 1:iterations) {
    seed <- i
    dat <- generate_data(NN, beta0, betaX, gamaZ, seed, Z_type)
    metrics <- calculate_metrics(dat)
    
    results$AUC0[i] <- metrics$AUC0
    results$AUC1[i] <- metrics$AUC1
    results$IDI[i] <- metrics$IDI_model
    results$NRI_categ[i] <- metrics$NRI_categ_model
    results$NRI_conti[i] <- metrics$NRI_conti_model
    results$Delong_p_value[i] <- metrics$Delong_p_value
    
    AUC0_list[i] <- metrics$AUC0
    AUC1_list[i] <- metrics$AUC1
    IDI_list[i] <- metrics$IDI_model
    NRI_categ_list[i] <- metrics$NRI_categ_model
    NRI_conti_list[i] <- metrics$NRI_conti_model
    
    if (metrics$Delong_p_value < 0.05) {
      Delong_significant_count <- Delong_significant_count + 1
    }
  }
  
  overall_stats <- list(
    AUC0_mean = mean(AUC0_list),
    AUC0_sd = sd(AUC0_list),
    
    AUC1_mean = mean(AUC1_list),
    AUC1_sd = sd(AUC1_list),
    
    IDI_mean = mean(IDI_list),
    IDI_sd = sd(IDI_list),
    
    NRI_categ_mean = mean(NRI_categ_list),
    NRI_categ_sd = sd(NRI_categ_list),
    
    NRI_conti_mean = mean(NRI_conti_list),
    NRI_conti_sd = sd(NRI_conti_list),
    Delong_significant_count = Delong_significant_count
  )
  
  list(results = results, overall_stats = overall_stats)
  
}

# Set parameters IFGIGT
NN <- 3000
beta0 <- -5   #50% prevalence
betaX <- c(-0.02, 0.15, 0.02, 0.02)
gamaZ <- 0.64
iterations <- 100

# Simulate metrics for different Z types
Z_types <- c("Z1", "Z2", "Z3")
results_list <- list()
overall_stats_list <- list()

# Create a workbook
wb <- createWorkbook()

for (Z_type in Z_types) {
  cat(paste("Z variable:", Z_type, "\n"))
  simulation_result <- simulate_metrics(iterations, NN, beta0, betaX, gamaZ, Z_type)
  results_list[[Z_type]] <- simulation_result$results
  overall_stats_list[[Z_type]] <- simulation_result$overall_stats
  
  # Add results to workbook
  addWorksheet(wb, paste0(Z_type, "_details"))
  writeData(wb, sheet = paste0(Z_type, "_details"), simulation_result$results)
  
  # Add overall stats to workbook
  addWorksheet(wb, paste0(Z_type, "_summary"))
  writeData(wb, sheet = paste0(Z_type, "_summary"), as.data.frame(simulation_result$overall_stats))
}

# Save workbook
saveWorkbook(wb, "E:/table/simulation_IFGIGT_Z_results0716.xlsx", overwrite = TRUE)

overall_stats_list
