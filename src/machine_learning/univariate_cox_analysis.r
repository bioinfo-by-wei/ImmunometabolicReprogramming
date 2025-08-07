var_genes <- colnames(bulk_data[, 3:ncol(bulk_data)])
univ_formulas <- sapply(
  var_genes,
  function(x) as.formula(paste("Surv(time, event)~", x))
)
univ_models <- lapply(
  univ_formulas,
  function(x) {
    coxph(x, data = bulk_data)
  }
)

univ_results <- lapply(
  univ_models,
  function(x) {
    x <- summary(x)
    p_value <- signif(x$wald["pvalue"], digits = 3)
    wald_test <- signif(x$wald["test"], digits = 3)
    beta <- signif(x$coef[1], digits = 3) # coeficient beta
    hr <- signif(x$coef[2], digits = 3) # exp(beta)
    hr_confint_lower <- signif(x$conf.int[, "lower .95"], 3)
    hr_confint_upper <- signif(x$conf.int[, "upper .95"], 3)
    hr <- paste0(
      hr, " (",
      hr_confint_lower, "-", hr_confint_upper, ")"
    )
    res <- c(beta, hr, wald_test, p_value)
    names(res) <- c("beta", "HR (95% CI for HR)", "wald_test", "p_value")
    res
  }
)
univ_results <- t(as.data.frame(univ_results, check.names = FALSE))
univ_results <- as.data.frame(univ_results)
univ_results$p_value <- as.numeric(univ_results$p_value)
results[[name]] <- univ_results
