library(dplyr)
library(log4r)
library(qs)
library(CoxBoost)

path <- "./step_pre_plot_survival_and_timeROC.r"
source(path) 
method_name <- "CoxBoost"
cindex_results <- data.frame()
cindex_result <- data.frame()

results <- list()

set.seed(1) 
pen <- optimCoxBoostPenalty(train_data[, "time"],
  train_data[, "event"],
  as.matrix(train_data[, -c(1, 2)]),
  trace = TRUE,
  start.penalty = 500,
  parallel = TRUE
)
cv_res <- cv.CoxBoost(train_data[, "time"],
  train_data[, "event"],
  as.matrix(train_data[, -c(1, 2)]),
  maxstepno = 500,
  K = 5,
  type = "verweij",
  penalty = pen$penalty,
  multicore = 1
)
cox_fit <- CoxBoost(train_data[, "time"],
  train_data[, "event"],
  as.matrix(train_data[, -c(1, 2)]),
  stepno = cv_res$optimal.step,
  penalty = pen$penalty
)
risk_socre <- as.numeric(predict(cox_fit,
  newdata = train_data[, -c(1, 2)],
  newtime = train_data[, "time"],
  newstatus = train_data[, "event"],
  type = "lp"
))
plot_survival_and_timeROC(
  train_data = train_data,
  pred_prob = risk_socre,
  store_name = paste0(store_path, "train", "_coxboost")
)
result <- calculate_roc_cindex(
  train_data = train_data,
  pred_prob = risk_socre,
  store_name = paste0(store_path, "train", "_coxboost")
)
cindex_result <- data.frame(
  Method = method_name,
  Seed = i,
  cindex = result$coxph_index
) %>% rename(!!train_data_name := cindex)
#Verify data
rs <- lapply(names(verify_data_list), function(name) {
  verify_data <- verify_data_list[[name]]
  risk_socre <- as.numeric(predict(cox_fit,
    newdata = verify_data[, -c(1, 2)],
    newtime = verify_data[, "time"],
    newstatus = verify_data[, "event"],
    type = "lp"
  ))
  plot_survival_and_timeROC(
    train_data = verify_data,
    pred_prob = risk_socre,
    store_name = paste0(store_path, name, "_coxboost_seed", i)
  )
  result <- calculate_roc_cindex(
    train_data = verify_data,
    pred_prob = risk_socre,
    store_name = paste0(store_path, name, "_coxboost_seed", i)
  )
  data.frame(
    Method = method_name,
    Seed = i,
    cindex = result$coxph_index
  ) %>% rename(!!name := cindex)
})
for (result in rs) {
  cindex_result <- left_join(cindex_result, result, by = c("Method", "Seed"))
}
cindex_results <- bind_rows(cindex_results, cindex_result)
# choose the xgboost genes
sig_gene_multi_cox <- c()
step_logplik <- predict(cox_fit,
  newdata = as.matrix(train_data[, -c(1, 2)]),
  newtime = train_data[, "time"],
  newstatus = train_data[, "event"],
  at.step = 0:cv_res$optimal.step,
  type = "logplik"
)
sig_gene_multi_cox <- cox_fit$xnames[cox_fit$coefficients[which.max(step_logplik), ] != 0]
if (length(sig_gene_multi_cox) == 0) {
  print(paste0("xgboost"," has no selected_genes."))
}
sig_gene_multi_cox <- c(sig_gene_multi_cox, "time", "event")
results[["xgboost"]] <- sig_gene_multi_cox


qsave(results, paste0(store_path, "coxboost_selected_genes.qs"))
write.csv(cindex_results, paste0(store_path, "coxboost_cindex_results.csv"))
