rm(list = ls())
library(qs)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(timeROC)
library(survival)
library(survminer)

path <- "./step_pre_plot_survival_and_timeROC.r"
source(path)  #import the common function

ori_data <- train_data

val_data_list <- verify_data_list
pre_var <- colnames(ori_data)[-c(1:2)]
estimate_data <- ori_data[, c("time", "event", pre_var)]
val_dd_list <- lapply(val_data_list, function(x) {
  x[, c("time", "event", pre_var)]
})

rf_nodesize <- 5
seed <- 1

cindex_results <- data.frame()
cindex_result <- data.frame()

##################################
#### Method RSF ####
##################################
method_name <- "RSF"
set.seed(seed)
fit <- rfsrc(Surv(time, event) ~ .,
  data = estimate_data,
  ntree = 1000, nodesize = rf_nodesize,
  splitrule = "logrank",
  importance = TRUE,
  proximity = TRUE,
  forest = TRUE,
  seed = seed
)

RS <- as.numeric(predict(fit, newdata = estimate_data)$predicted)
store_name <- paste0(store_path, "training", "_RSF_")
result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
cindex_result <- data.frame(
  Method = method_name,
  cindex = result$coxph_index
) %>% rename(!!train_data_name := cindex)
plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
})

rs <- lapply(names(rs), function(name) {
  train_data <- rs[[name]]
  store_name <- paste0(store_path, name, "_rsf")
  write.csv(train_data, paste(store_name, "_predict.csv"))
  plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
  result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
  data.frame(
    Method = method_name,
    cindex = result$coxph_index
  ) %>% rename(!!name := cindex)
})

for (result in rs) {
  cindex_result <- left_join(cindex_result, result, by = c("Method"))
}
cindex_results <- bind_rows(cindex_results, cindex_result)

##################################
#### Method Enet ####
##################################
method_name <- "Enet"
x1 <- as.matrix(estimate_data[, pre_var])
x2 <- as.matrix(Surv(estimate_data$time, estimate_data$event))

for (alpha in seq(0, 1, 0.1)) {
  set.seed(seed)
  fit <- cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 30, gamma = 0)
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))
  })
  RS <- as.numeric(predict(fit, type = "link", newx = as.matrix(x1), s = fit$lambda.min))
  store_name <- paste0(store_path, "training", "_enet_", alpha)
  train_result <- calculate_roc_cindex(train_data, RS, store_name = store_name)
  plot_survival_and_timeROC(train_data, RS, store_name = store_name)

  cindex_result <- data.frame(
    Method = paste0(method_name, "[α=", alpha, "]"),
    cindex = train_result$coxph_index
  ) %>% rename(!!train_data_name := cindex)

  rs <- lapply(names(rs), function(name) {
    train_data <- rs[[name]]
    store_name <- paste0(store_path, name, "_enet_", alpha)
    result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
    plot_survival_and_timeROC(
      train_data,
      train_data$RS,
      store_name = store_name
    )
    data.frame(
      Method = paste0(method_name, "[α=", alpha, "]"),
      cindex = result$coxph_index
    ) %>% rename(!!name := cindex)
  })

  for (result in rs) {
    cindex_result <- left_join(cindex_result, result, by = c("Method"))
  }
  cindex_results <- bind_rows(cindex_results, cindex_result)
}

##################################
#### Method StepCox ####
##################################
method_name <- "StepCox "
for (direction_name in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time, event) ~ ., estimate_data), direction = direction_name)
  rs <- lapply(val_dd_list, function(x) {
    cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
  })

  RS <- as.numeric(predict(fit, type = "risk", newdata = estimate_data))
  store_name <- paste0(store_path, "training", "_StepCox_", direction_name)
  result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
  cindex_result <- data.frame(
    Method = paste0(method_name, "[", direction_name, "]"),
    cindex = result$coxph_index
  ) %>% rename(!!train_data_name := cindex)
  plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

  rs <- lapply(names(rs), function(name) {
    train_data <- rs[[name]]
    store_name <- paste0(store_path, name, "_stepcox_", direction_name)
    write.csv(train_data, paste(store_name, "_predict.csv"))
    result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
    plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
    data.frame(
      Method = paste0(method_name, "[", direction_name, "]"),
      cindex = result$coxph_index
    ) %>% rename(!!name := cindex)
  })

  for (result in rs) {
    cindex_result <- left_join(cindex_result, result, by = c("Method"))
  }
  cindex_results <- bind_rows(cindex_results, cindex_result)
}
##################################
#### Method StepCox+survivalsvm ####
##################################
method_name <- "StepCox+survivalsvm"
for (direction_name in c("both", "backward")) {
  fit <- step(coxph(Surv(time, event) ~ ., estimate_data), direction = direction_name)
  rid <- names(coef(fit))
  estimate_data2 <- estimate_data[, c("time", "event", rid)]
  val_dd_list2 <- lapply(val_data_list, function(x) {
    x[, c("time", "event", rid)]
  })

  fit <- survivalsvm(Surv(time, event) ~ ., data = estimate_data2, gamma.mu = 1)

  RS <- as.numeric(predict(fit, estimate_data)$predicted)
  store_name <- paste0(store_path, "training", "_StepCoxSurv_", direction_name)
  result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
  cindex_result <- data.frame(
    Method = paste0(method_name, "[", direction_name, "]"),
    cindex = result$coxph_index
  ) %>% rename(!!train_data_name := cindex)
  plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

  rs <- lapply(val_dd_list2, function(x) {
    cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
  })

  rs <- lapply(names(rs), function(name) {
    train_data <- rs[[name]]
    store_name <- paste0(store_path, name, "_stepcox_surv_", direction_name)
    write.csv(train_data, paste(store_name, "_predict.csv"))
    plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
    result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
    data.frame(
      Method = paste0(method_name, "[", direction_name, "]"),
      cindex = result$coxph_index
    ) %>% rename(!!name := cindex)
  })
  for (result in rs) {
    cindex_result <- left_join(cindex_result, result, by = c("Method"))
  }
  cindex_results <- bind_rows(cindex_results, cindex_result)
}

##################################
#### Method superpc####
##################################
method_name <- "superpc"
data <- list(
  x = t(estimate_data[, -c(1, 2)]),
  y = estimate_data$time,
  censoring.status = estimate_data$event,
  featurenames = colnames(estimate_data)[-c(1, 2)]
)
set.seed(seed)
fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5) # default
cv.fit <- superpc.cv(fit, data,
  n.threshold = 20, # default
  n.fold = 10,
  n.components = 3,
  min.features = 5,
  max.features = nrow(data$x),
  compute.fullcv = TRUE,
  compute.preval = TRUE
)

ff <- superpc.predict(fit, data,
  data,
  threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])],
  n.components = 1
)
RS <- as.numeric(ff$v.pred)

store_name <- paste0(store_path, "training", "_superpc")
result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
cindex_result <- data.frame(
  Method = method_name,
  cindex = result$coxph_index
) %>% rename(!!train_data_name := cindex)
plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

rs <- lapply(val_dd_list, function(w) {
  test <- list(x = t(w[, -c(1, 2)]), y = w$time, censoring.status = w$event, featurenames = colnames(w)[-c(1, 2)])
  ff <- superpc.predict(fit, data,
    test,
    threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])],
    n.components = 1
  )
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  rr2
})

rs <- lapply(names(rs), function(name) {
  train_data <- rs[[name]]
  store_name <- paste0(store_path, name, "_superpc")
  write.csv(train_data, paste(store_name, "_predict.csv"))
  plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
  result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
  data.frame(
    Method = method_name,
    cindex = result$coxph_index
  ) %>% rename(!!name := cindex)
})
for (result in rs) {
  cindex_result <- left_join(cindex_result, result, by = c("Method"))
}
cindex_results <- bind_rows(cindex_results, cindex_result)
##################################
#### Method GBM ####
##################################
method_name <- "GBM"
set.seed(seed)
fit <- gbm(
  formula = Surv(time, event) ~ ., data = estimate_data, distribution = "coxph",
  n.trees = 10000,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  cv.folds = 10, n.cores = 6
)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(
  formula = Surv(time, event) ~ ., data = estimate_data, distribution = "coxph",
  n.trees = best,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  cv.folds = 10, n.cores = 8
)

RS <- as.numeric(predict(fit, estimate_data, n.trees = best, type = "link"))
store_name <- paste0(store_path, "training", "_GBM")
result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
cindex_result <- data.frame(
  Method = method_name,
  cindex = result$coxph_index
) %>% rename(!!train_data_name := cindex)
plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
})

rs <- lapply(names(rs), function(name) {
  train_data <- rs[[name]]
  store_name <- paste0(store_path, name, "_gbm")
  write.csv(train_data, paste(store_name, "_predict.csv"))
  plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
  result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
  data.frame(
    Method = method_name,
    cindex = result$coxph_index
  ) %>% rename(!!name := cindex)
})
for (result in rs) {
  cindex_result <- left_join(cindex_result, result, by = c("Method"))
}
cindex_results <- bind_rows(cindex_results, cindex_result)
##################################
#### Method survivalsvm ####
##################################
method_name <- "survivalsvm"
fit <- survivalsvm(Surv(time, event) ~ ., data = estimate_data, gamma.mu = 1)

RS <- as.numeric(predict(fit, estimate_data)$predicted)
store_name <- paste0(store_path, "training", "_survivalsvm")
result <- calculate_roc_cindex(estimate_data, RS, store_name = store_name)
cindex_result <- data.frame(
  Method = method_name,
  cindex = result$coxph_index
) %>% rename(!!train_data_name := cindex)
plot_survival_and_timeROC(estimate_data, RS, store_name = store_name)

rs <- lapply(val_dd_list, function(x) {
  cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
})

rs <- lapply(names(rs), function(name) {
  train_data <- rs[[name]]
  store_name <- paste0(store_path, name, "_survivalsvm")
  write.csv(train_data, paste(store_name, "_predict.csv"))
  result <- calculate_roc_cindex(train_data, train_data$RS, store_name = store_name)
  plot_survival_and_timeROC(train_data, train_data$RS, store_name = store_name)
  data.frame(
    Method = method_name,
    cindex = result$coxph_index
  ) %>% rename(!!name := cindex)
})
for (result in rs) {
  cindex_result <- left_join(cindex_result, result, by = c("Method"))
}
cindex_results <- bind_rows(cindex_results, cindex_result)

#save the final result
write.csv(cindex_result, "final_cindex_results.csv")
