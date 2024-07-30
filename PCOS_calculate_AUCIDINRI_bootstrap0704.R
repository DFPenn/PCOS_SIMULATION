library(readxl)
library(dplyr)
library(pROC)
library(PredictABEL)
library(Hmisc)
library(boot)  # For bootstrap

# 读取数据
data <- read_excel("E:/table/IDI_NRI_article/PCOS/PCOSdata154622.xlsx")

rep <- 100

# 定义分组函数并赋值
menstrual_group <- function(cycle) {
  if (cycle >= 26 & cycle <= 35) {
    return(0)  # Eumeno
  } else if (cycle >= 36 & cycle < 90) {
    return(1)  # Oligo
  } else {
    return(2)  # Ameno
  }
}

data$menstrual_cycle_group <- sapply(data$menstrual_cycle2, menstrual_group)
data$menstrual_cycle_group <- factor(data$menstrual_cycle_group, levels = c(0, 1, 2), labels = c("Eumeno", "Oligo", "Ameno"))

print(table(data$menstrual_cycle_group))

# 定义计算AUC, IDI, NRI的函数
analyze_outcome <- function(data, outcome) {
  # 基线模型
  baseline_model <- glm(as.formula(paste(outcome, "~ Age + BMI + WC + FAI")), data = data, family = binomial)
  baseline_pred <- predict(baseline_model, data, type = "response")
  baseline_roc <- roc(data[[outcome]], baseline_pred)
  baseline_auc <- auc(baseline_roc)
  
  # 新模型
  new_model <- glm(as.formula(paste(outcome, "~ Age + BMI + WC + FAI + menstrual_cycle_group")), data = data, family = binomial)
  new_pred <- predict(new_model, data, type = "response")
  new_roc <- roc(data[[outcome]], new_pred)
  new_auc <- auc(new_roc)
  
  # DeLong Test
  delong_test <- roc.test(baseline_roc, new_roc, method = "delong")
  
  # 计算IDI和NRI
  roc.p <- pROC::roc(data[[outcome]], baseline_pred)  # Find the optimal cutoffB
  cutoffB <- roc.p$thresholds[which.max(roc.p$sensitivities + roc.p$specificities)]
  
  c1 <- cut(baseline_pred, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  c2 <- cut(new_pred, breaks = c(0, cutoffB, 1), include.lowest = TRUE, right = FALSE)
  
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  
  categ <- improveProb(x1 = as.numeric(c11) * (1 / (length(levels(c11)))), 
                       x2 = as.numeric(c22) * (1 / (length(levels(c22)))), 
                       y = data[[outcome]])
  
  conti <- improveProb(x1 = baseline_pred, x2 = new_pred, y = data[[outcome]])
  
  idi <- round(conti$idi, 4)
  categ_nri <- round(categ$nri, 4)
  conti_nri <- round(conti$nri, 4)
  
  return(c(baseline_auc, new_auc, idi, categ_nri, conti_nri, delong_test$p.value))
}

# 定义bootstrap函数
bootstrap_fun <- function(data, indices, outcome) {
  d <- data[indices, ]  
  analyze_outcome(d, outcome)
}

# 执行Bootstrap分析的函数
repeat_analysis <- function(data, outcome, n_iter = rep) {
  # 执行Bootstrap
  results <- boot(data = data, statistic = function(data, indices) bootstrap_fun(data, indices, outcome), R = n_iter)
  
  mean_baseline_auc <- mean(results$t[, 1])
  mean_new_auc <- mean(results$t[, 2])
  mean_idi <- mean(results$t[, 3])
  mean_categ_nri <- mean(results$t[, 4])
  mean_conti_nri <- mean(results$t[, 5])
  delong_significant <- mean(results$t[, 6] < 0.05)
  
  return(list(
    mean_baseline_auc = mean_baseline_auc,
    mean_new_auc = mean_new_auc,
    mean_idi = mean_idi,
    mean_categ_nri = mean_categ_nri,
    mean_conti_nri = mean_conti_nri,
    delong_significant = delong_significant
  ))
}

# 定义结局变量列表
outcomes <- c("IR", "IFGIGT", "dyslipidemia")

# 遍历每个结局变量并分析
results <- lapply(outcomes, function(outcome) {
  repeat_analysis(data, outcome)
})

# 打印结果
for (i in 1:length(outcomes)) {
  cat(paste("Outcome:", outcomes[i], "\n"))
  cat(paste("Mean Baseline AUC:", results[[i]]$mean_baseline_auc, "\n"))
  cat(paste("Mean New Model AUC:", results[[i]]$mean_new_auc, "\n"))
  cat(paste("Mean IDI:", results[[i]]$mean_idi, "\n"))
  cat(paste("Mean Categorical NRI:", results[[i]]$mean_categ_nri, "\n" ))
  cat(paste("Mean Continuous NRI:", results[[i]]$mean_conti_nri, "\n" ))
  cat(paste("Proportion of significant DeLong tests:", results[[i]]$delong_significant, "\n"))
  cat("\n")
}
