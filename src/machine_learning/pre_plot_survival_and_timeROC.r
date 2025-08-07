plot_survival_and_timeROC <- function(train_data, pred_prob, store_name) {
  # grouping: divide the risk scores into high- and low-risk groups based on the median.
  pred_prob <- as.numeric(pred_prob)
  risk_score_plot <- ifelse(pred_prob > median(pred_prob), "High", "Low")
  # build the plot database
  risk_data_plot <- data.frame(
    time = train_data$time,
    event = train_data$event,
    risk_score_plot = factor(
      risk_score_plot,
      levels = c("High", "Low"),
      labels = c("High", "Low")
    )
  )

  # Survival Analysis Model
  fit <- survfit(Surv(time, as.numeric(event)) ~ risk_score_plot, data = risk_data_plot)

  # plot the Kaplan–Meier curve
  lasso_km <- ggsurvplot(
    fit,
    data = risk_data_plot,
    legend.labs = c("High", "Low"),
    pval = TRUE,
    pval.size = 8,
    risk.table = TRUE,
    size = 0.6, 
    risk.table.height = 0.25,
    surv.median.line = "hv",
    legend.title = "RiskScore",
    title = "Overall survival",
    ylab = "Survival probability",
    xlab = "Time (Days)",
    censor.shape = 124,
    censor.size = 2,
    conf.int = FALSE,
    break.x.by = 720,
    ggtheme = theme_classic() +  
      theme(
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18),   
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18)
      )
  )

  pdf(paste0(store_name, "_km.pdf"), width = 8, height = 12)
  print(lasso_km)
  dev.off()

  # TimeROC
  times <- c(365, 2 * 365, 3 * 365, 4 * 365, 5 * 365)
  time_roc <- timeROC(
    T = train_data$time,
    delta = train_data$event,
    marker = pred_prob,
    cause = 1,
    times = times,
    iid = TRUE
  )

  # timeROC
  pdf(paste0(store_name, "_timeroc.pdf"), width = 8, height = 8)
  par(cex.lab = 1.8)
  par(cex.axis = 1.8)
  plot(time_roc, time = 365, col = "#1f77b4", title = FALSE, lwd = 2)
  plot(time_roc, time = 2 * 365, add = TRUE, lwd = 1.5, col = "#ff7f0e")
  plot(time_roc, time = 3 * 365, add = TRUE, lwd = 1.5, col = "#2ca02c")
  plot(time_roc, time = 4 * 365, add = TRUE, lwd = 1.5, col = "#d62728")
  plot(time_roc, time = 5 * 365, add = TRUE, lwd = 1.5, col = "#9467bd")
  legend("bottomright",
         legend = paste0(
           c("1-year", "2-year", "3-year", "4-year", "5-year"),
           " AUC=", round(time_roc$AUC, 3)
         ),
         col = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"),
         lwd = 2,
         cex = 1.8)
  dev.off()
}

calculate_roc_cindex <- function(train_data, pred_prob, store_name, is_looger = FALSE, logger = NULL) {
  roc_obj <- roc(train_data$event, pred_prob)

  concordance_cindex <- concordance.index(
    x = pred_prob,
    surv.time = train_data$time,
    surv.event = train_data$event
  )

  cal_data <- data.frame(time = train_data$time, event = train_data$event, pred_prob)

  coxph_cindex <- as.numeric(summary(coxph(Surv(time, event) ~ pred_prob, cal_data))$concordance[1])

  write.csv(
    cal_data,
    paste0(store_name, "_predict_riskscore.csv")
  )
  print(store_name)
  print(paste0(
    "coxph-cindex_value:",
    coxph_cindex
  ))
  print(paste0("concordance-index_value:", concordance_cindex$c.index))
  print(paste0("auc_value:",  as.numeric(auc(roc_obj))))
  print(store_name)
  if (is_looger) {
    info(logger, "coxph-index:", coxph_cindex)
    info(logger, "concordance:", as.numeric(concordance_cindex$c.index))
    info(logger, "auc:", as.numeric(auc(roc_obj)))
  }
  result <- data.frame(
    coxph_index = coxph_cindex,
    concordance = as.numeric(concordance_cindex$c.index),
    auc = as.numeric(auc(roc_obj))
  )
  write.csv(
    result,
    paste0(store_name, "_calculated_cindex_resulst.csv")
  )
  return(result)
}
