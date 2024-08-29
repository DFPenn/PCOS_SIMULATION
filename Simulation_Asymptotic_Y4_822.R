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


# Function to calculate AUC, IDI, NRI, and Delong test, including theoretical variance calculation
calculate_metrics <- function(dat) {
  m0 <- glm(y ~ X1 + X2 + X3 + X4, family = binomial, data = dat)
  m1 <- glm(y ~ X1 + X2 + X3 + X4 + Z, family = binomial, data = dat)
  
  yhat0 <- m0$fitted.values
  yhat1 <- m1$fitted.values
  
  roc_model0 <- roc(dat$y, yhat0)
  roc_model1 <- roc(dat$y, yhat1)
  
  AUC0 <- roc_model0$auc
  AUC1 <- roc_model1$auc
  Delta_AUC <- AUC1 - AUC0
  
  delong_test_result <- roc.test(roc_model0, roc_model1)
  Delong_p_value <- delong_test_result$p.value
  
  # Calculate the theoretical variance of IDI and NRI
  delta_p <- yhat1 - yhat0
  delta_p_case <- delta_p[dat$y == 1]
  delta_p_con <- delta_p[dat$y == 0]
  
  Asy_var_IDI <- var(delta_p_case) / length(delta_p_case) + var(delta_p_con) / length(delta_p_con)
  Asy_var_IDI_sd <- sqrt(Asy_var_IDI)
  
  Asy_var_NRI <- 4 * (var(ifelse(delta_p_case > 0, 1, 0)) / length(delta_p_case) + 
                        var(ifelse(delta_p_con > 0, 1, 0)) / length(delta_p_con))
  Asy_var_NRI_sd <- sqrt(Asy_var_NRI)
  
  # Calculate the theoretical variance of AUC
  Asy_var_AUC0_sd <- sqrt(var(roc_model0, method = "delong"))
  Asy_var_AUC1_sd <- sqrt(var(roc_model1, method = "delong"))
  Asy_deltaAUC_sd <- sqrt(Asy_var_AUC1_sd^2 + Asy_var_AUC0_sd^2)
  
  # NRI and IDI calculations with the function
  categ <- improveProb(x1 = yhat0, x2 = yhat1, y = dat$y)
  conti <- improveProb(x1 = yhat0, x2 = yhat1, y = dat$y)
  
  IDI_model <- mean(delta_p_case) - mean(delta_p_con)
  NRI_categ_model <- round(categ$nri, 4)
  NRI_conti_model <- round(conti$nri, 4)
  
  list(AUC0 = AUC0, Asy_var_AUC0_sd = Asy_var_AUC0_sd, AUC1 = AUC1, Asy_var_AUC1_sd = Asy_var_AUC1_sd, 
       Asy_deltaAUC_sd = Asy_deltaAUC_sd, Delta_AUC = Delta_AUC, Delong_p_value = Delong_p_value,
       IDI_model = IDI_model, Asy_var_IDI_sd = Asy_var_IDI_sd, Asy_var_NRI_sd = Asy_var_NRI_sd,
       NRI_categ_model = NRI_categ_model, NRI_conti_model = NRI_conti_model)
}

# Main simulation
simulate_metrics <- function(iterations, NN, beta0, betaX, gamaZ, Z_type) {
  results <- data.frame(
    iteration = 1:iterations,
    AUC0 = numeric(iterations),
    AUC1 = numeric(iterations),
    DeltaAUC = numeric(iterations),
    Asy_var_AUC0_sd = numeric(iterations),
    Asy_var_AUC1_sd = numeric(iterations),
    Asy_deltaAUC_sd = numeric(iterations),
    IDI = numeric(iterations),
    Asy_var_IDI_sd = numeric(iterations),
    Asy_var_NRI_sd = numeric(iterations),
    NRI_categ = numeric(iterations),
    NRI_conti = numeric(iterations),
    Delong_p_value = numeric(iterations)
  )
  
  AUC0_list <- numeric(iterations)
  AUC1_list <- numeric(iterations)
  DeltaAUC_list <- numeric(iterations)
  Asy_var_AUC0_sd_list <- numeric(iterations)
  Asy_var_AUC1_sd_list <- numeric(iterations)
  Asy_deltaAUC_sd_list <- numeric(iterations)
  IDI_list <- numeric(iterations)
  Asy_var_IDI_sd_list <- numeric(iterations)  
  Asy_var_NRI_sd_list <- numeric(iterations)
  NRI_categ_list <- numeric(iterations)
  NRI_conti_list <- numeric(iterations)
  Delong_significant_count <- 0
  
  for (i in 1:iterations) {
    
    seed <- i+2024
    dat <- generate_data(NN, beta0, betaX, gamaZ, seed, Z_type)
    metrics <- calculate_metrics(dat)
    
    results[i,] <- c(i, metrics$AUC0, metrics$AUC1, metrics$Delta_AUC, metrics$Asy_var_AUC0_sd,
                     metrics$Asy_var_AUC1_sd, metrics$Asy_deltaAUC_sd, metrics$IDI_model, 
                     metrics$Asy_var_IDI_sd, metrics$Asy_var_NRI_sd, metrics$NRI_categ_model, 
                     metrics$NRI_conti_model, metrics$Delong_p_value)
    
    AUC0_list[i] <- metrics$AUC0
    AUC1_list[i] <- metrics$AUC1
    DeltaAUC_list[i] <- metrics$Delta_AUC
    Asy_var_AUC0_sd_list[i] <- metrics$Asy_var_AUC0_sd
    Asy_var_AUC1_sd_list[i] <- metrics$Asy_var_AUC1_sd
    Asy_deltaAUC_sd_list[i] <- metrics$Asy_deltaAUC_sd
    IDI_list[i] <- metrics$IDI_model
    Asy_var_IDI_sd_list[i] <- metrics$Asy_var_IDI_sd
    Asy_var_NRI_sd_list[i] <- metrics$Asy_var_NRI_sd
    NRI_categ_list[i] <- metrics$NRI_categ_model
    NRI_conti_list[i] <- metrics$NRI_conti_model
    
    if (metrics$Delong_p_value < 0.05) {
      Delong_significant_count <- Delong_significant_count + 1
    }
  }
  
  overall_stats <- data.frame(
    AUC0_mean = mean(AUC0_list),
    AUC0_sd = sd(AUC0_list),
    Asy_var_AUC0_sd_mean = mean(Asy_var_AUC0_sd_list), 
    AUC1_mean = mean(AUC1_list),
    AUC1_sd = sd(AUC1_list),
    Asy_var_AUC1_sd_mean = mean(Asy_var_AUC1_sd_list), 
    Asy_deltaAUC_sd_mean = mean(Asy_deltaAUC_sd_list),
    DeltaAUC_sd = sd(DeltaAUC_list),
    IDI_mean = mean(IDI_list),
    IDI_sd = sd(IDI_list),
    Asy_var_IDI_sd_mean = mean(Asy_var_IDI_sd_list),  
    NRI_categ_mean = mean(NRI_categ_list),
    NRI_categ_sd = sd(NRI_categ_list),
    NRI_conti_mean = mean(NRI_conti_list),
    NRI_conti_sd = sd(NRI_conti_list),
    Asy_var_NRI_sd_mean = mean(Asy_var_NRI_sd_list),
    Delong_significant_count = Delong_significant_count
  )
  
  list(results = results, overall_stats = overall_stats)
}



# Set parameters IFGIGT
NN <- 3000
beta0 <- 0.97
betaX <- c(-0.01, -0.09, 0.01, -0.02)
gamaZ <- 0.2
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
saveWorkbook(wb, "E:/table/822_Simulation_Y4dyslipidemia_asymptotic_IDINRIAUC01.xlsx", overwrite = TRUE)


